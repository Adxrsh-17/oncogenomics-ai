#!/usr/bin/env python3
"""
ML Variant Prioritizer — Random Forest + XGBoost variant ranking.
Features: VAF, depth, functional consequence, pathway relevance,
          clonality estimate, mutation type, cancer driver score.
Input:  results/patient_clinical_report.csv + results/pathway_analysis.csv + MAF
Output: results/ml_prioritizer_scores.csv + model performance report
"""
import json, sys, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

WORK      = Path(__file__).resolve().parent.parent / "results"
MAF_FILE  = WORK / "tcga_luad_merged.maf"
REPORT    = WORK / "patient_clinical_report.csv"
OUT_SCORES = WORK / "ml_prioritizer_scores.csv"
OUT_MODEL  = WORK / "ml_model_report.json"

# ---------- known cancer driver genes for scoring ----------
TIER1_DRIVERS = {
    "KRAS", "EGFR", "TP53", "STK11", "KEAP1", "RB1", "CDKN2A",
    "SMARCA4", "NF1", "BRAF", "MET", "ERBB2", "PIK3CA", "PTEN",
}
TIER2_DRIVERS = {
    "ATM", "BRCA2", "RET", "ALK", "ROS1", "NTRK1", "NTRK2", "NTRK3",
    "FGFR1", "FGFR2", "FGFR3", "CDH1", "ARID1A", "KDM6A", "SMAD4",
}

PATHWAY_MAP = {
    "KRAS": "RAS/MAPK", "NRAS": "RAS/MAPK", "BRAF": "RAS/MAPK", "NF1": "RAS/MAPK",
    "EGFR": "RTK/Growth", "MET": "RTK/Growth", "ERBB2": "RTK/Growth",
    "PIK3CA": "PI3K/AKT", "PTEN": "PI3K/AKT", "AKT1": "PI3K/AKT",
    "TP53": "Cell Cycle", "RB1": "Cell Cycle", "CDKN2A": "Cell Cycle",
    "BRCA1": "DNA Repair", "BRCA2": "DNA Repair", "ATM": "DNA Repair",
    "SMARCA4": "Chromatin", "ARID1A": "Chromatin", "KDM6A": "Chromatin",
    "APC": "Wnt", "CTNNB1": "Wnt",
    "STK11": "Other", "KEAP1": "Other",
}

ACTIONABLE_VARIANTS = {
    "KRAS_G12C": ("Sotorasib", "FDA-approved", "A"),
    "EGFR_L858R": ("Osimertinib", "FDA-approved", "A"),
    "EGFR_ex19del": ("Osimertinib", "FDA-approved", "A"),
    "BRAF_V600E": ("Dabrafenib+Trametinib", "FDA-approved", "A"),
    "MET_ex14skip": ("Capmatinib", "FDA-approved", "A"),
    "RET_fusion": ("Selpercatinib", "FDA-approved", "A"),
    "ALK_fusion": ("Alectinib", "FDA-approved", "A"),
    "ERBB2_amp": ("Trastuzumab", "FDA-approved B", "B"),
    "PIK3CA_mut": ("Alpelisib (trial)", "Phase II", "B"),
    "BRCA2_mut": ("Olaparib (trial)", "Phase II", "B"),
}

def consequence_score(variant_classification):
    """Map variant classification to numerical severity (0–10)."""
    scores = {
        "Nonsense_Mutation": 9, "Frame_Shift_Del": 9, "Frame_Shift_Ins": 9,
        "Splice_Site": 8, "Translation_Start_Site": 8,
        "Missense_Mutation": 7, "In_Frame_Del": 6, "In_Frame_Ins": 6,
        "Nonstop_Mutation": 7,
        "Silent": 1, "3'UTR": 1, "5'UTR": 1, "Intron": 0,
        "RNA": 0, "lincRNA": 0, "IGR": 0,
    }
    return scores.get(variant_classification, 3)

def driver_tier(gene):
    if gene in TIER1_DRIVERS: return 2
    if gene in TIER2_DRIVERS: return 1
    return 0

def build_feature_matrix(maf_df):
    """Engineer features from MAF dataframe for ML prioritization."""
    df = maf_df.copy()

    # VAF
    if "tumor_f" in df.columns:
        df["vaf"] = pd.to_numeric(df["tumor_f"], errors="coerce")
    elif "t_alt_count" in df.columns and "t_depth" in df.columns:
        depth_vals = pd.to_numeric(df["t_depth"], errors="coerce")
        depth_vals = depth_vals.where(depth_vals != 0, np.nan)
        df["vaf"] = pd.to_numeric(df["t_alt_count"], errors="coerce") / depth_vals
    else:
        df["vaf"] = 0.15  # median imputation

    # Depth
    if "t_depth" in df.columns:
        df["depth"] = pd.to_numeric(df["t_depth"], errors="coerce").fillna(30)
    else:
        df["depth"] = 30

    # Consequence score
    df["consequence_score"] = df["Variant_Classification"].apply(consequence_score)

    # Driver gene score
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else "Gene"
    df["driver_tier_legacy"] = df[gene_col].apply(driver_tier)
    
    # Load OncoKB for curated driver_tier
    df["driver_tier"] = 0
    oncokb_path = Path(__file__).resolve().parent.parent / "data" / "oncokb_cancer_gene_list.txt"
    if oncokb_path.exists():
        oncokb_df = pd.read_csv(oncokb_path, sep="\t")
        oncokb_genes = set(oncokb_df[oncokb_df["OncoKB Annotated"] == "Yes"]["Hugo Symbol"])
        df["driver_tier"] = df[gene_col].apply(lambda x: 2 if x in oncokb_genes else 0)
    else:
        df["driver_tier"] = df["driver_tier_legacy"]

    # Pathway membership
    df["pathway"] = df[gene_col].map(PATHWAY_MAP).fillna("Other")

    # Clonality estimate (VAF > 0.35 → likely clonal in diploid)
    df["clonality_score"] = (df["vaf"].fillna(0) > 0.35).astype(float)

    # Mutation type
    df["is_snv"] = (df["Variant_Type"] == "SNP").astype(float) if "Variant_Type" in df.columns else 1.0

    # Known hotspot (regex)
    hotspot_set = {"G12C", "G12V", "G12D", "L858R", "V600E", "R175H", "R248W", "R248Q", "R273H"}
    if "HGVSp_Short" in df.columns:
        df["is_hotspot_regex"] = df["HGVSp_Short"].fillna("").apply(
            lambda x: float(any(h in str(x) for h in hotspot_set))
        )
    else:
        df["is_hotspot_regex"] = 0.0

    # New curated source (CancerHotspots)
    df["is_hotspot"] = 0.0
    hotspots_path = Path(__file__).resolve().parent.parent / "data" / "hotspots.csv"
    if hotspots_path.exists() and "HGVSp_Short" in df.columns:
        import re
        hotspots_df = pd.read_csv(hotspots_path)
        hotspot_locs = set(zip(hotspots_df["Hugo_Symbol"], hotspots_df["Amino_Acid_Position"].astype(str)))
        
        def check_curated_hotspot(row):
            hgvsp = str(row.get("HGVSp_Short", ""))
            gene = str(row.get(gene_col, ""))
            if not str(hgvsp).startswith("p."): return 0.0
            match = re.search(r"p\.[A-Za-z]+(\d+)", str(hgvsp))
            if match:
                pos = match.group(1)
                return float((gene, pos) in hotspot_locs)
            return 0.0
            
        df["is_hotspot"] = df.apply(check_curated_hotspot, axis=1)

    # Population frequency proxy (ExAC/gnomAD — use 0 if not available as germline filter)
    if "ExAC_AF" in df.columns:
        df["pop_af"] = pd.to_numeric(df["ExAC_AF"], errors="coerce").fillna(0)
    else:
        df["pop_af"] = 0.0

    features = [
        "vaf", "depth", 
        "clonality_score", "is_snv", "is_hotspot", "pop_af",
    ]
    return df, features, gene_col

def train_rf_xgb(X, y_labels, feature_names, groups):
    """Train Random Forest and XGBoost on binary labels (driver=1, passenger=0)."""
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import GroupKFold, cross_val_score
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import roc_auc_score
    import numpy as np

    # Replace NaN
    X = np.nan_to_num(X, nan=0.0, posinf=1.0, neginf=0.0)

    print(f"Training on {X.shape[0]} variants, {X.shape[1]} features")
    print(f"  Driver variants (y=1): {y_labels.sum()} / {len(y_labels)}")

    rf = RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42, n_jobs=6)
    xgb = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.05, random_state=42)

    cv = GroupKFold(n_splits=5)

    if y_labels.sum() < 5:
        print("  Warning: too few positive labels for CV — using full dataset")
        rf.fit(X, y_labels)
        xgb.fit(X, y_labels)
        rf_auc = 0.0
        xgb_auc = 0.0
    else:
        rf_scores  = cross_val_score(rf,  X, y_labels, groups=groups, cv=cv, scoring="roc_auc")
        xgb_scores = cross_val_score(xgb, X, y_labels, groups=groups, cv=cv, scoring="roc_auc")
        rf_auc  = rf_scores.mean()
        xgb_auc = xgb_scores.mean()
        print(f"  RF  CV AUC: {rf_auc:.3f} ± {rf_scores.std():.3f}")
        print(f"  XGB CV AUC: {xgb_auc:.3f} ± {xgb_scores.std():.3f}")
        rf.fit(X, y_labels)
        xgb.fit(X, y_labels)

    fi_rf  = dict(zip(feature_names, rf.feature_importances_))
    fi_xgb = dict(zip(feature_names, xgb.feature_importances_))

    return rf, xgb, rf_auc, xgb_auc, fi_rf, fi_xgb

def main():
    print("=== ML Variant Prioritizer ===")

    # Load MAF
    if MAF_FILE.exists():
        print(f"Loading MAF: {MAF_FILE}")
        maf = pd.read_csv(MAF_FILE, sep="\t", comment="#", low_memory=False)
    else:
        # Try the full MAF if available
        alt_maf = Path("/opt/ml/input/data/notebook/oncoseq/mafs/tcga_luad_full_merged.maf")
        if alt_maf.exists():
            print(f"Loading full MAF: {alt_maf}")
            maf = pd.read_csv(alt_maf, sep="\t", comment="#", low_memory=False)
        else:
            print("ERROR: MAF file not found")
            sys.exit(1)

    print(f"MAF: {len(maf)} rows, {maf['Tumor_Sample_Barcode'].nunique() if 'Tumor_Sample_Barcode' in maf.columns else '?'} patients")

    # Build feature matrix
    df, features, gene_col = build_feature_matrix(maf)

    # Labels: tier-1 driver gene missense/nonsense = positive class
    y = ((df["driver_tier"] >= 1) & (df["consequence_score"] >= 7)).astype(int).values
    X = df[features].values
    groups = df[gene_col].values

    # Train models
    rf, xgb_model, rf_auc, xgb_auc, fi_rf, fi_xgb = train_rf_xgb(X, y, features, groups)

    # Score all variants
    rf_prob  = rf.predict_proba(np.nan_to_num(X, nan=0.0))[:, 1]
    xgb_prob = xgb_model.predict_proba(np.nan_to_num(X, nan=0.0))[:, 1]
    df["rf_score"]  = rf_prob
    df["xgb_score"] = xgb_prob
    df["ensemble_score"] = (rf_prob + xgb_prob) / 2

    # Actionability annotation
    def get_action(row):
        gene = str(row.get(gene_col, ""))
        hgvsp = str(row.get("HGVSp_Short", ""))
        for key, (drug, approval, tier) in ACTIONABLE_VARIANTS.items():
            g, mut = key.split("_", 1) if "_" in key else (key, "")
            if gene == g and (not mut or mut in hgvsp):
                return drug, approval, tier
        return "", "", ""

    df[["drug", "approval", "tier"]] = df.apply(
        lambda r: pd.Series(get_action(r)), axis=1
    )

    # Save scored variants (top 5000 by ensemble score)
    out_cols = [gene_col, "Tumor_Sample_Barcode", "HGVSp_Short", "Variant_Classification",
                "pathway"] + features + ["rf_score", "xgb_score", "ensemble_score", "drug", "approval", "tier"]
    out_cols = [c for c in out_cols if c in df.columns]

    top_variants = df.sort_values("ensemble_score", ascending=False).head(5000)
    top_variants[out_cols].to_csv(OUT_SCORES, index=False)
    print(f"\nScored variants saved: {OUT_SCORES}")

    # Top 20 for display
    print("\nTop 20 prioritized variants:")
    disp_cols = [gene_col, "HGVSp_Short", "Variant_Classification", "ensemble_score", "drug"]
    disp_cols = [c for c in disp_cols if c in top_variants.columns]
    print(top_variants[disp_cols].head(20).to_string(index=False))

    # Model report
    report = {
        "rf_cv_auc": round(rf_auc, 4),
        "xgb_cv_auc": round(xgb_auc, 4),
        "n_variants": len(df),
        "n_positive_class": int(y.sum()),
        "features": features,
        "rf_feature_importance": {k: round(float(v), 4) for k, v in fi_rf.items()},
        "xgb_feature_importance": {k: round(float(v), 4) for k, v in fi_xgb.items()},
        "actionable_variants_found": int((df["drug"] != "").sum()),
    }
    with open(OUT_MODEL, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nModel report: {OUT_MODEL}")
    print(f"RF AUC: {rf_auc:.3f}  XGB AUC: {xgb_auc:.3f}")
    print(f"Actionable variants: {report['actionable_variants_found']}")
    
    # Report agreement between regex and curated hotspot
    agreement = (df["is_hotspot"] == df["is_hotspot_regex"]).mean()
    print(f"Hotspot Agreement (Regex vs Curated): {agreement:.2%}")
    print("=== DONE ===")

if __name__ == "__main__":
    # Install sklearn if needed
    try:
        from sklearn.ensemble import RandomForestClassifier
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "scikit-learn", "xgboost", "-q"], check=True)
    main()
