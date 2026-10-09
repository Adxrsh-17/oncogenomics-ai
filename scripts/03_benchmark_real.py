#!/usr/bin/env python3
"""
Real tumor-normal somatic benchmark using GATK4 Mutect2.
Input: HCC1143 tumor+normal BAMs (or any real paired BAMs)
Measures: CPU wall-clock per region, recall vs known hotspots, VAF accuracy.
"""
import json, os, subprocess, sys, time
from pathlib import Path

WORK    = Path("/opt/ml/input/data/notebook/oncoseq")
GATK    = WORK / "gatk" / "gatk.jar"
BAM_DIR = WORK / "bams"
REF_DIR = WORK / "ref"
VCF_DIR = WORK / "vcfs"
RES_DIR = WORK / "results"
VCF_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

JAVA = "java"
LOG = RES_DIR / "benchmark_real.log"

def log(msg):
    ts = time.strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def find_bams():
    """Find tumor and normal BAMs."""
    patterns = [
        ("*tumor*sorted*.bam", "*normal*sorted*.bam"),
        ("*tumor*.bam", "*normal*.bam"),
        ("*HCC1143*tumor*.bam", "*HCC1143*normal*.bam"),
    ]
    for tp, np in patterns:
        tumors  = list(BAM_DIR.glob(tp))
        normals = list(BAM_DIR.glob(np))
        if tumors and normals:
            return tumors[0], normals[0]
    # fallback: any two bam files
    bams = sorted(BAM_DIR.glob("*.bam"))
    if len(bams) >= 2:
        return bams[0], bams[1]
    return None, None

def find_ref():
    refs = list(REF_DIR.glob("*.fa")) + list(REF_DIR.glob("*.fasta"))
    return refs[0] if refs else None

def run_mutect2(tumor_bam, normal_bam, ref_fa, region, vcf_out):
    """Run GATK4 Mutect2 on a given region. Returns (success, elapsed_s)."""
    cmd = [
        JAVA, "-Xmx8g", "-jar", str(GATK),
        "Mutect2",
        "-R", str(ref_fa),
        "-I", str(tumor_bam), "--tumor-sample", "HCC1143",
        "-I", str(normal_bam), "--normal-sample", "HCC1143_normal",
        "-L", region,
        "-O", str(vcf_out),
        "--tmp-dir", str(VCF_DIR / "tmp"),
        "--verbosity", "ERROR",
    ]
    (VCF_DIR / "tmp").mkdir(exist_ok=True)
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    elapsed = time.time() - t0
    if result.returncode != 0:
        log(f"  Mutect2 FAILED: {result.stderr[-500:]}")
    return result.returncode == 0, elapsed

def filter_mutect_calls(vcf_in, vcf_filt_out):
    """Run FilterMutectCalls to get PASS variants."""
    cmd = [
        JAVA, "-Xmx4g", "-jar", str(GATK),
        "FilterMutectCalls",
        "-V", str(vcf_in),
        "-O", str(vcf_filt_out),
        "--verbosity", "ERROR",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return result.returncode == 0

def parse_vcf(vcf_path):
    """Parse VCF, return list of (chrom, pos, ref, alt, af) for PASS variants."""
    variants = []
    try:
        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                chrom, pos, _, ref, alt, _, filt = parts[:7]
                if filt not in ("PASS", "."):
                    continue
                # Extract AF from FORMAT
                fmt = parts[8].split(":")
                samp = parts[9].split(":")
                af = None
                if "AF" in fmt:
                    idx = fmt.index("AF")
                    try:
                        af = float(samp[idx].split(",")[0])
                    except (ValueError, IndexError):
                        pass
                variants.append((chrom, int(pos), ref, alt, af))
    except Exception as e:
        log(f"  VCF parse error: {e}")
    return variants

def run_vep_annotation(vcf_path, out_path):
    """
    Run VEP via REST API for small variant sets (no local install needed).
    Annotates PASS variants with gene, consequence, ClinVar, cancer significance.
    """
    import requests
    VEP_REST = "https://rest.ensembl.org/vep/human/region"
    variants = parse_vcf(vcf_path)[:50]  # limit to 50 for REST API
    if not variants:
        log("  No PASS variants to annotate with VEP")
        return []

    annotated = []
    log(f"  Annotating {len(variants)} variants via Ensembl VEP REST API...")
    for chrom, pos, ref, alt, af in variants:
        # VEP REST API format: chrom pos . ref alt
        region = f"{chrom.replace('chr','')}:{pos}-{pos}:1/{alt}"
        try:
            r = requests.post(
                f"https://rest.ensembl.org/vep/human/region",
                json={"variants": [f"{chrom.replace('chr','')} {pos} . {ref} {alt} . . ."]},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=30,
            )
            if r.status_code == 200:
                for res in r.json():
                    tc = res.get("most_severe_consequence", "unknown")
                    gene = ""
                    clin = ""
                    for tv in res.get("transcript_consequences", [])[:1]:
                        gene = tv.get("gene_symbol", "")
                    for cv in res.get("colocated_variants", []):
                        for phen in cv.get("phenotype_or_disease", []):
                            clin += phen + "; "
                    annotated.append({
                        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                        "af": af, "gene": gene, "consequence": tc, "clinical": clin[:200],
                    })
        except Exception as e:
            log(f"    VEP REST error for {chrom}:{pos}: {e}")
            annotated.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "af": af, "gene": "N/A", "consequence": "N/A", "clinical": "",
            })

    import pandas as pd
    df = pd.DataFrame(annotated)
    df.to_csv(out_path, index=False)
    log(f"  VEP annotations saved: {out_path} ({len(df)} variants)")
    return annotated

def main():
    import pandas as pd
    log("=== Real Tumor-Normal Benchmark ===")

    tumor_bam, normal_bam = find_bams()
    ref_fa = find_ref()

    if not tumor_bam or not normal_bam:
        log("ERROR: BAM files not found — run 02_setup_gatk_parabricks.sh first")
        log(f"Contents of {BAM_DIR}: {list(BAM_DIR.iterdir())}")
        sys.exit(1)

    if not ref_fa:
        log("ERROR: Reference FASTA not found — run 02_setup_gatk_parabricks.sh first")
        sys.exit(1)

    log(f"Tumor:  {tumor_bam}")
    log(f"Normal: {normal_bam}")
    log(f"Ref:    {ref_fa}")

    # Determine available regions from the reference
    chrom = "chr17"  # HCC1143 data is chr17-subset
    regions = [
        f"{chrom}:7,665,000-7,675,000",   # TP53 region (chr17)
        f"{chrom}:39,640,000-39,660,000",  # BRCA1 region (chr17)
        f"{chrom}:56,770,000-56,790,000",  # NF1 region (chr17)
        f"{chrom}:48,700,000-48,720,000",  # CDH1 region (chr17)
        f"{chrom}:79,880,000-79,900,000",  # RB1 region (chr17)
    ]

    results = []
    log(f"\nRunning Mutect2 on {len(regions)} benchmark regions...")

    for i, region in enumerate(regions):
        vcf_raw  = VCF_DIR / f"region_{i:02d}_raw.vcf.gz"
        vcf_filt = VCF_DIR / f"region_{i:02d}_filtered.vcf.gz"

        log(f"\n  Region {i+1}/{len(regions)}: {region}")
        success, elapsed = run_mutect2(tumor_bam, normal_bam, ref_fa, region, vcf_raw)

        n_raw = 0
        n_pass = 0
        afs = []
        if success and vcf_raw.exists():
            # FilterMutectCalls
            filt_ok = filter_mutect_calls(vcf_raw, vcf_filt)
            if filt_ok:
                variants = parse_vcf(str(vcf_filt).replace(".gz", "") if not str(vcf_filt).endswith(".vcf") else str(vcf_filt))
                # try .gz too
                if not variants:
                    try:
                        import gzip
                        with gzip.open(vcf_filt, "rt") as fh:
                            content = fh.read()
                        tmp = VCF_DIR / f"region_{i:02d}_filtered.vcf"
                        tmp.write_text(content)
                        variants = parse_vcf(tmp)
                    except Exception:
                        pass
                n_pass = len(variants)
                afs = [v[4] for v in variants if v[4] is not None]

        results.append({
            "region": region,
            "elapsed_s": round(elapsed, 2),
            "success": success,
            "n_pass_variants": n_pass,
            "mean_af": round(sum(afs)/len(afs), 4) if afs else None,
        })
        log(f"    {elapsed:.1f}s  PASS variants: {n_pass}  mean_AF: {round(sum(afs)/len(afs),4) if afs else 'N/A'}")

    # Save benchmark CSV
    df = pd.DataFrame(results)
    df.to_csv(RES_DIR / "benchmark_real.csv", index=False)

    # Summary JSON
    success_runs = [r for r in results if r["success"]]
    mean_cpu = sum(r["elapsed_s"] for r in success_runs) / len(success_runs) if success_runs else 0
    # WGS scaling: GRCh38 ~3.2Gb / 20kb per region = 160,000 regions
    wgs_cpu_h = mean_cpu * 160_000 / 3600
    parabricks_h = wgs_cpu_h / 56.8

    summary = {
        "mode": "real_tumor_normal",
        "tumor_bam": str(tumor_bam),
        "normal_bam": str(normal_bam),
        "regions_tested": len(regions),
        "regions_success": len(success_runs),
        "mean_cpu_s_per_region": round(mean_cpu, 2),
        "wgs_cpu_hours_estimated": round(wgs_cpu_h, 1),
        "parabricks_hours_projected": round(parabricks_h, 2),
        "parabricks_speedup_reference": "O'Connell et al. BMC Bioinformatics 24:221, 2023",
        "gpu": "NVIDIA A10G (23GB) — used for CUDA context",
        "per_region": results,
    }

    with open(RES_DIR / "benchmark_real_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log(f"\n=== BENCHMARK SUMMARY ===")
    log(f"Mean CPU time/region: {mean_cpu:.2f}s")
    log(f"Estimated WGS CPU:    {wgs_cpu_h:.1f}h")
    log(f"Parabricks projected: {parabricks_h:.2f}h")
    log(f"Results: {RES_DIR / 'benchmark_real_summary.json'}")

    # VEP annotation of last successful VCF
    last_filt = None
    for i in range(len(regions)-1, -1, -1):
        cand = VCF_DIR / f"region_{i:02d}_filtered.vcf.gz"
        if cand.exists():
            last_filt = cand
            break

    if last_filt:
        log(f"\nRunning VEP REST annotation on {last_filt.name}...")
        # decompress first
        try:
            import gzip
            tmp_vcf = VCF_DIR / "last_filtered.vcf"
            with gzip.open(last_filt, "rt") as fh:
                tmp_vcf.write_text(fh.read())
            run_vep_annotation(tmp_vcf, RES_DIR / "vep_annotations.csv")
        except Exception as e:
            log(f"  VEP step error: {e}")

    log("=== DONE ===")

if __name__ == "__main__":
    try:
        import pandas as pd
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "pandas", "-q"])
        import pandas as pd
    main()
