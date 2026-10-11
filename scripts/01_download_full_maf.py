#!/usr/bin/env python3
"""
Full TCGA-LUAD MAF download — targeting all 557 patients.
Uses GDC bulk API. Runs on /opt/ml/input/data/notebook/oncoseq/mafs/
"""
import json, os, requests, subprocess, sys, time
from pathlib import Path

LOCAL_RESULTS = Path(__file__).resolve().parent.parent / "results"
LOCAL_WORKDIR = LOCAL_RESULTS / "mafs"
SAGEMAKER_WORKDIR = Path("/opt/ml/input/data/notebook/oncoseq/mafs")
WORKDIR = SAGEMAKER_WORKDIR if Path("/opt/ml").exists() else LOCAL_WORKDIR
WORKDIR.mkdir(parents=True, exist_ok=True)
LOG = WORKDIR / "download.log"

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA  = "https://api.gdc.cancer.gov/data"

def log(msg):
    ts = time.strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def get_all_luad_maf_ids():
    """Fetch all open-access TCGA-LUAD MAF file UUIDs via GDC API."""
    params = {
        "filters": json.dumps({
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-LUAD"}},
                {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}},
                {"op": "=", "content": {"field": "access", "value": "open"}},
                {"op": "=", "content": {"field": "data_format", "value": "MAF"}},
            ]
        }),
        "fields": "file_id,file_name,cases.case_id,cases.submitter_id",
        "format": "json",
        "size": "1000",
    }
    r = requests.get(GDC_FILES, params=params, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    log(f"GDC returned {len(hits)} MAF file entries")
    return hits

def download_bulk(file_ids, outdir):
    """Download a batch of file UUIDs as a tarball."""
    ids_str = "\n".join(file_ids)
    cmd = [
        "curl", "-s", "--retry", "5", "--retry-delay", "5",
        "--max-time", "600",
        "-H", "Content-Type: application/json",
        "-X", "POST",
        "-d", json.dumps({"ids": file_ids}),
        "-o", str(outdir / "batch.tar.gz"),
        GDC_DATA,
    ]
    result = subprocess.run(cmd, capture_output=True, timeout=700)
    return result.returncode == 0

def extract_and_merge(mafs_dir):
    """Extract all .maf.gz files and merge into one file."""
    import gzip, io
    import pandas as pd

    log("Extracting MAF files...")
    all_tarballs = list(mafs_dir.glob("batch_*.tar.gz"))
    all_maf_files = []

    # Extract each tarball
    for tb in all_tarballs:
        extract_dir = mafs_dir / f"extracted_{tb.stem}"
        extract_dir.mkdir(exist_ok=True)
        subprocess.run(["tar", "-xzf", str(tb), "-C", str(extract_dir)], check=True)
        all_maf_files.extend(extract_dir.rglob("*.maf.gz"))
        all_maf_files.extend(extract_dir.rglob("*.maf"))

    log(f"Found {len(all_maf_files)} MAF files to merge")

    dfs = []
    for maf_path in all_maf_files:
        try:
            opener = gzip.open if str(maf_path).endswith(".gz") else open
            with opener(maf_path, "rt") as fh:
                lines = []
                for line in fh:
                    if not line.startswith("#"):
                        lines.append(line)
            df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False)
            dfs.append(df)
        except Exception as e:
            log(f"  Skip {maf_path.name}: {e}")

    if not dfs:
        log("ERROR: no MAF dataframes loaded")
        return None

    merged = pd.concat(dfs, ignore_index=True)
    out = mafs_dir / "tcga_luad_full_merged.maf"
    merged.to_csv(out, sep="\t", index=False)
    log(f"Merged MAF: {len(merged)} rows, {merged['Tumor_Sample_Barcode'].nunique()} patients → {out}")
    return out

def main():
    log("=== TCGA-LUAD Full MAF Download ===")

    # Check if already downloaded
    merged_path = WORKDIR / "tcga_luad_full_merged.maf"
    if merged_path.exists():
        import pandas as pd
        df = pd.read_csv(merged_path, sep="\t", low_memory=False, nrows=5)
        log(f"Already merged: {merged_path}")
        return

    log("Step 1: Get file IDs from GDC...")
    hits = get_all_luad_maf_ids()
    file_ids = [h["file_id"] for h in hits]
    log(f"Total MAF files: {len(file_ids)}")

    # Save manifest
    with open(WORKDIR / "manifest.json", "w") as f:
        json.dump(hits, f, indent=2)

    # Download in batches of 50
    BATCH = 50
    batches = [file_ids[i:i+BATCH] for i in range(0, len(file_ids), BATCH)]
    log(f"Downloading {len(batches)} batches of ≤{BATCH} files each...")

    for idx, batch in enumerate(batches):
        batch_out = WORKDIR / f"batch_{idx:03d}.tar.gz"
        if batch_out.exists():
            log(f"  Batch {idx+1}/{len(batches)}: already downloaded, skip")
            continue
        log(f"  Batch {idx+1}/{len(batches)}: {len(batch)} files...")
        ids_str = "\n".join(batch)
        cmd = (
            f'curl -s --retry 5 --retry-delay 5 --max-time 600 '
            f'-H "Content-Type: application/json" -X POST '
            f'-d \'{{"ids": {json.dumps(batch)}}}\' '
            f'-o {batch_out} '
            f'{GDC_DATA}'
        )
        ret = subprocess.run(cmd, shell=True, timeout=700)
        if ret.returncode == 0 and batch_out.stat().st_size > 1000:
            log(f"    OK ({batch_out.stat().st_size/1e6:.1f} MB)")
        else:
            log(f"    FAILED (code {ret.returncode})")

    log("Step 3: Extract and merge...")
    extract_and_merge(WORKDIR)
    log("=== DONE ===")

if __name__ == "__main__":
    try:
        import pandas as pd
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "pandas", "-q"], check=True)
        import pandas as pd
    main()
