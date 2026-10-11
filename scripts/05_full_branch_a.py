#!/usr/bin/env python3
"""
Full Branch A analysis on 566-patient TCGA-LUAD MAF.
Produces updated patient counts, pathway analysis, clinical report
matching the planned 557+ patient cohort.
"""
import json, sys
from pathlib import Path
import pandas as pd
import numpy as np

import os
LOCAL_MAF = Path(__file__).resolve().parent.parent / "results" / "tcga_luad_merged.maf"
DEFAULT_MAF = str(LOCAL_MAF) if LOCAL_MAF.exists() else "/opt/ml/input/data/notebook/oncoseq/mafs/tcga_luad_full_merged.maf"
LOCAL_OUT = Path(__file__).resolve().parent.parent / "results" / "sagemaker"
DEFAULT_OUT = str(LOCAL_OUT) if not Path("/opt/ml").exists() else "/opt/ml/input/data/notebook/oncoseq/results"

MAF_FULL = Path(os.environ.get("ONCOSEQ_MAF", DEFAULT_MAF))
OUT_DIR  = Path(os.environ.get("ONCOSEQ_OUT", DEFAULT_OUT))
OUT_DIR.mkdir(parents=True, exist_ok=True)

GNOMAD_MAX_AF_THRESH = 0.001

PATHWAY_GENES = {
    "RAS/MAPK": ["KRAS","NRAS","HRAS","BRAF","RAF1","MAP2K1","MAP2K2","NF1"],
    "RTK/Growth": ["EGFR","ERBB2","ERBB3","MET","FGFR1","FGFR2","FGFR3","RET","ALK","ROS1"],
    "PI3K/AKT": ["PIK3CA","PIK3R1","PTEN","AKT1","AKT2","TSC1","TSC2","MTOR"],
    "Cell Cycle": ["TP53","RB1","CDKN2A","CDK4","CDK6","CCND1","MDM2","SMAD4"],
    "Chromatin": ["SMARCA4","ARID1A","ARID2","KDM6A","KMT2C","KMT2D","SETD2","EP300"],
    "DNA Repair": ["BRCA1","BRCA2","ATM","ATR","CHEK1","CHEK2","PALB2","FANCM"],
    "Wnt": ["APC","CTNNB1","RNF43","AXIN1","AXIN2","TCF7L2"],
}

DRUG_MAP = {
    "G12C":  ("KRAS",  "Sotorasib",    "FDA-approved", "CodeBreaK 100 — 37.1% RR"),
    "L858R": ("EGFR",  "Osimertinib",  "FDA-approved", "FLAURA — 80% RR"),
    "V600E": ("BRAF",  "Dabrafenib+Trametinib", "FDA-approved", "GEOMETRY — 44% RR"),
}

def identify_targets(maf):
    """Find KRAS G12C and EGFR L858R patients (exact gene + exact protein change)."""
    hgvs = maf["HGVSp_Short"].fillna("")
    gene = maf["Hugo_Symbol"]
    kras_mask = (gene == "KRAS") & (hgvs == "p.G12C")
    egfr_mask = (gene == "EGFR") & (hgvs == "p.L858R")
    kras_pts = set(maf.loc[kras_mask, "Tumor_Sample_Barcode"].dropna().unique())
    egfr_pts = set(maf.loc[egfr_mask, "Tumor_Sample_Barcode"].dropna().unique())
    return kras_pts, egfr_pts

def pathway_analysis(maf, patients):
    """Count patients with ≥1 mutation in each pathway."""
    maf_filt = maf[maf["Tumor_Sample_Barcode"].isin(patients)]
    gene_col = "Hugo_Symbol"
    results = []
    for pw, genes in PATHWAY_GENES.items():
        pw_maf = maf_filt[maf_filt[gene_col].isin(genes)]
        pts = pw_maf["Tumor_Sample_Barcode"].nunique()
        top_gene = pw_maf[gene_col].value_counts().idxmax() if len(pw_maf) > 0 else "N/A"
        top_n = pw_maf[gene_col].value_counts().max() if len(pw_maf) > 0 else 0
        results.append({
            "Pathway": pw, "Patients": pts,
            "Pct": round(100 * pts / len(patients), 1) if patients else 0,
            "Top_Gene": top_gene, "Top_Gene_N": int(top_n),
        })
    return sorted(results, key=lambda x: -x["Patients"])

def mutation_burden(maf):
    """Mutations per patient."""
    return maf.groupby("Tumor_Sample_Barcode").size()

def main():
    print(f"Loading MAF: {MAF_FULL}")
    maf = pd.read_csv(MAF_FULL, sep="\t", comment="#", low_memory=False)
    print(f"  {len(maf)} mutations, {maf['Tumor_Sample_Barcode'].nunique()} patients")

    # Gap 4: Explicit gnomAD filter
    if "gnomAD_AF" in maf.columns:
        af = pd.to_numeric(maf["gnomAD_AF"], errors="coerce").fillna(0.0)
        pre_len = len(maf)
        maf = maf[af <= GNOMAD_MAX_AF_THRESH].copy()
        removed = pre_len - len(maf)
        print(f"  Filtered {removed} variants with gnomAD_AF > {GNOMAD_MAX_AF_THRESH}")

    all_pts = set(maf["Tumor_Sample_Barcode"].dropna().unique())
    N = len(all_pts)
    print(f"  Total patients: {N}")

    # Mutation burden
    burden = mutation_burden(maf)
    median_burden = int(burden.median())
    print(f"  Median mutation burden: {median_burden}")

    # Target identification
    kras_pts, egfr_pts = identify_targets(maf)
    union_pts = kras_pts | egfr_pts
    print(f"  KRAS G12C: {len(kras_pts)} patients ({100*len(kras_pts)/N:.1f}%)")
    print(f"  EGFR L858R: {len(egfr_pts)} patients ({100*len(egfr_pts)/N:.1f}%)")
    print(f"  Target cohort: {len(union_pts)} patients")

    # Pathway analysis
    pw = pathway_analysis(maf, all_pts)
    pw_df = pd.DataFrame(pw)
    pw_df.to_csv(OUT_DIR / "pathway_analysis_full.csv", index=False)
    print(f"\nPathway analysis:")
    for row in pw:
        print(f"  {row['Pathway']:15s}: {row['Patients']:3d} pts ({row['Pct']:.1f}%) — top: {row['Top_Gene']} ({row['Top_Gene_N']})")

    # Clinical report for target patients
    target_maf = maf[maf["Tumor_Sample_Barcode"].isin(union_pts)].copy()
    target_maf["variant"] = target_maf["HGVSp_Short"].fillna("") if "HGVSp_Short" in target_maf.columns else ""
    target_maf["is_kras_g12c"] = target_maf["Tumor_Sample_Barcode"].isin(kras_pts)
    target_maf["is_egfr_l858r"] = target_maf["Tumor_Sample_Barcode"].isin(egfr_pts)

    reports = []
    for pt in sorted(union_pts):
        b = int(burden.get(pt, 0))
        is_kras = pt in kras_pts
        is_egfr = pt in egfr_pts
        drug = " + ".join(d for d, f in (("Sotorasib", is_kras), ("Osimertinib", is_egfr)) if f)
        variant = " + ".join(v for v, f in (("KRAS p.G12C", is_kras), ("EGFR p.L858R", is_egfr)) if f)
        trial = " + ".join(t for t, f in (("CodeBreaK 100 (37.1% RR)", is_kras), ("FLAURA (80% RR)", is_egfr)) if f)
        reports.append({
            "patient": pt, "variant": variant, "drug": drug,
            "trial": trial, "mutation_burden": b,
            "high_burden_flag": b > 500,
        })

    rep_df = pd.DataFrame(reports).sort_values("variant")
    rep_df.to_csv(OUT_DIR / "patient_clinical_report_full.csv", index=False)
    print(f"\nClinical report: {len(rep_df)} actionable patients")

    # Summary JSON
    summary = {
        "cohort": "TCGA-LUAD full open-access",
        "total_patients": N,
        "total_mutations": len(maf),
        "median_mutation_burden": median_burden,
        "kras_g12c_patients": len(kras_pts),
        "kras_g12c_pct": round(100*len(kras_pts)/N, 1),
        "egfr_l858r_patients": len(egfr_pts),
        "egfr_l858r_pct": round(100*len(egfr_pts)/N, 1),
        "target_cohort_union": len(union_pts),
        "pathway_summary": pw,
    }
    with open(OUT_DIR / "oncoseq_summary_full.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nSummary saved: {OUT_DIR / 'oncoseq_summary_full.json'}")
    print("=== DONE ===")

if __name__ == "__main__":
    main()
