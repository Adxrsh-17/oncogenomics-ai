#!/usr/bin/env python3
"""
Structured ClinVar Pathogenicity Parsing & External Model Validation.
Parses ClinVar annotations into standard ACMG 5-tier classifications:
  - Pathogenic (P)
  - Likely Pathogenic (LP)
  - Variant of Uncertain Significance (VUS)
  - Likely Benign (LB)
  - Benign (B)

Evaluates the trained ML prioritizer scores against ClinVar P/LP vs B/LB
as an independent, external benchmark (addressing Gap 7 and Stage 9 audit gap).
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
MAF_FILE = RESULTS_DIR / "tcga_luad_merged.maf"
SCORES_FILE = RESULTS_DIR / "ml_prioritizer_scores.csv"
OUT_JSON = RESULTS_DIR / "clinvar_validation_metrics.json"

def parse_clinvar_string(raw_sig):
    """
    Parse a raw ClinVar significance string (from MAF CLIN_SIG or VEP)
    into a standardized classification (P, LP, VUS, LB, B, Other).
    """
    if pd.isna(raw_sig) or not str(raw_sig).strip():
        return "Unannotated"
    s = str(raw_sig).lower().replace(" ", "_").replace("/", ";")
    
    # Check for conflicting interpretations
    if "conflicting" in s:
        return "Conflicting"
    
    # Priority matching
    if "pathogenic" in s and "likely_pathogenic" not in s:
        return "P"
    if "likely_pathogenic" in s:
        return "LP"
    if "uncertain_significance" in s or "vus" in s:
        return "VUS"
    if "likely_benign" in s:
        return "LB"
    if "benign" in s:
        return "B"
    return "Other"

def parse_vep_colocated_clinvar(colocated_variants):
    """
    Extract structured ClinVar classification from Ensembl VEP REST colocated_variants.
    """
    if not colocated_variants or not isinstance(colocated_variants, list):
        return "Unannotated"
    for cv in colocated_variants:
        # Check clin_sig list or clinvar_somatic_classification
        sigs = cv.get("clin_sig", [])
        if sigs:
            for sig in sigs:
                parsed = parse_clinvar_string(sig)
                if parsed in ("P", "LP", "VUS", "LB", "B"):
                    return parsed
        somatic = cv.get("clinvar_somatic_classification", "")
        if somatic:
            parsed = parse_clinvar_string(somatic)
            if parsed in ("P", "LP", "VUS", "LB", "B"):
                return parsed
    return "Unannotated"

def run_clinvar_validation():
    print("=" * 60)
    print("ClinVar Structured Parsing & External ML Validation")
    print("=" * 60)

    if not MAF_FILE.exists():
        print(f"ERROR: MAF file {MAF_FILE} not found.")
        return

    # Load MAF CLIN_SIG
    print("Loading MAF CLIN_SIG column...")
    maf = pd.read_csv(
        MAF_FILE,
        sep="\t",
        comment="#",
        low_memory=False,
        usecols=["Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short", "CLIN_SIG"]
    )

    maf["clinvar_tier"] = maf["CLIN_SIG"].apply(parse_clinvar_string)
    dist = maf["clinvar_tier"].value_counts().to_dict()
    print("ClinVar Classification Distribution across MAF:")
    for k, v in sorted(dist.items(), key=lambda x: -x[1]):
        print(f"  {k:15s}: {v:5d} variants")

    # Load ML Scores
    if not SCORES_FILE.exists():
        print(f"Scores file {SCORES_FILE} not found; skipping external metric evaluation.")
        return

    print("\nMerging with ML prioritizer scores...")
    scores = pd.read_csv(SCORES_FILE)
    
    # Merge on patient + gene + HGVSp_Short
    merged = pd.merge(
        scores,
        maf[["Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short", "clinvar_tier"]].drop_duplicates(),
        on=["Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short"],
        how="inner"
    )

    print(f"Merged scored variants with ClinVar: {len(merged)}")
    
    # Score distribution per tier
    tier_stats = {}
    for tier, grp in merged.groupby("clinvar_tier"):
        tier_stats[tier] = {
            "count": int(len(grp)),
            "rf_mean": float(round(grp["rf_score"].mean(), 4)),
            "rf_median": float(round(grp["rf_score"].median(), 4)),
            "xgb_mean": float(round(grp["xgb_score"].mean(), 4)),
            "xgb_median": float(round(grp["xgb_score"].median(), 4)),
            "ensemble_mean": float(round(grp["ensemble_score"].mean(), 4)),
        }

    print("\nModel Scores across ClinVar Tiers:")
    for t, s in tier_stats.items():
        print(f"  Tier {t:6s} (n={s['count']:3d}) -> Ensemble Mean: {s['ensemble_mean']:.4f}, RF Mean: {s['rf_mean']:.4f}")

    # Binary discrimination: Pathogenic (P or LP) vs Non-Pathogenic (B or LB)
    eval_set = merged[merged["clinvar_tier"].isin(["P", "LP", "B", "LB"])].copy()
    eval_set["is_pathogenic"] = eval_set["clinvar_tier"].isin(["P", "LP"]).astype(int)

    metrics = {
        "clinvar_maf_distribution": dist,
        "scored_tier_stats": tier_stats,
    }

    if eval_set["is_pathogenic"].nunique() == 2 and len(eval_set) >= 10:
        y_true = eval_set["is_pathogenic"]
        rf_auc = float(roc_auc_score(y_true, eval_set["rf_score"]))
        rf_prauc = float(average_precision_score(y_true, eval_set["rf_score"]))
        xgb_auc = float(roc_auc_score(y_true, eval_set["xgb_score"]))
        xgb_prauc = float(average_precision_score(y_true, eval_set["xgb_score"]))
        ens_auc = float(roc_auc_score(y_true, eval_set["ensemble_score"]))
        ens_prauc = float(average_precision_score(y_true, eval_set["ensemble_score"]))

        print("\nExternal Validation on ClinVar (P/LP vs B/LB):")
        print(f"  Evaluation set size: {len(eval_set)} (Pathogenic: {int(y_true.sum())}, Benign: {int(len(y_true) - y_true.sum())})")
        print(f"  Random baseline PR-AUC: {y_true.mean():.4f}")
        print(f"  Random Forest ROC-AUC : {rf_auc:.4f} | PR-AUC: {rf_prauc:.4f}")
        print(f"  XGBoost       ROC-AUC : {xgb_auc:.4f} | PR-AUC: {xgb_prauc:.4f}")
        print(f"  Ensemble      ROC-AUC : {ens_auc:.4f} | PR-AUC: {ens_prauc:.4f}")

        metrics["external_validation"] = {
            "sample_size": int(len(eval_set)),
            "pathogenic_count": int(y_true.sum()),
            "benign_count": int(len(y_true) - y_true.sum()),
            "baseline_prevalence": float(round(y_true.mean(), 4)),
            "rf_roc_auc": round(rf_auc, 4),
            "rf_pr_auc": round(rf_prauc, 4),
            "xgb_roc_auc": round(xgb_auc, 4),
            "xgb_pr_auc": round(xgb_prauc, 4),
            "ensemble_roc_auc": round(ens_auc, 4),
            "ensemble_pr_auc": round(ens_prauc, 4),
        }
    else:
        print("\nNote: Insufficient B/LB variants in scored subset for ROC-AUC calculation.")

    with open(OUT_JSON, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nSaved metrics: {OUT_JSON}")
    print("=" * 60)

if __name__ == "__main__":
    run_clinvar_validation()
