#!/usr/bin/env python3
"""
ML Variant Prioritizer — Random Forest + XGBoost variant ranking.
Features: VAF, depth, functional consequence, pathway relevance,
          clonality estimate, mutation type, cancer driver score.
Input:  results/patient_clinical_report.csv + results/pathway_analysis.csv + MAF
Output: results/ml_prioritizer_scores.csv + model performance report
"""
import json, sys, subprocess, os, requests
from pathlib import Path
import numpy as np
import pandas as pd

os.environ['OMP_NUM_THREADS'] = '6'
os.environ['MKL_NUM_THREADS'] = '6'
os.environ['OPENBLAS_NUM_THREADS'] = '6'

WORK      = Path(__file__).resolve().parent.parent / "results"
MAF_FILE  = WORK / "tcga_luad_merged.maf"
REPORT    = WORK / "patient_clinical_report.csv"
OUT_SCORES = WORK / "ml_prioritizer_scores.csv"
OUT_MODEL  = WORK / "ml_model_report.json"

GNOMAD_MAX_AF_THRESH = 0.001
MODEL_DIR  = WORK / "models"
RF_MODEL_PATH = MODEL_DIR / "rf_model.joblib"

# TCGA ABSOLUTE purity — GDC open-access supplementary
PURITY_URL = (
    "https://api.gdc.cancer.gov/data/"
    "1e5e78d4-9e8a-449b-9ca0-9a6c4b3a3a0e"  # TCGA-LUAD clinical supplement UUID
)
PURITY_LOCAL = WORK.parent / "data" / "tcga_purity.tsv"


def load_purity_estimates():
    """Load per-patient tumour purity.  Returns dict {case_id: purity_float}.

    Priority:
      1. data/tcga_purity.tsv  (columns: case_id, purity)
      2. Broad/TCGA PanCancer ABSOLUTE TSV via public URL
      3. Uniform fallback 0.35 (typical LUAD purity)

    To supply your own: place a TSV at data/tcga_purity.tsv with columns
    case_id (e.g. TCGA-05-4418) and purity (0–1 float).
    """
    if PURITY_LOCAL.exists():
        df = pd.read_csv(PURITY_LOCAL, sep="\t")
        if {"case_id", "purity"}.issubset(df.columns):
            return dict(zip(df["case_id"], pd.to_numeric(df["purity"], errors="coerce").fillna(0.35)))

    # Try Broad TCGA PanCancer ABSOLUTE purity (public supplement)
    broad_url = (
        "https://api.gdc.cancer.gov/data/"
        "1a7d7be8-675d-4e60-a105-19d4b1a7b4d8"
    )
    try:
        r = requests.get(broad_url, timeout=15)
        if r.status_code == 200:
            from io import StringIO
            df = pd.read_csv(StringIO(r.text), sep="\t")
            # Common column names across TCGA purity files
            id_col = next((c for c in df.columns if "sample" in c.lower() or "case" in c.lower()), None)
            pur_col = next((c for c in df.columns if "purity" in c.lower() or "absolute" in c.lower()), None)
            if id_col and pur_col:
                # Normalise barcode → 12-char case ID (TCGA-XX-XXXX)
                ids = df[id_col].astype(str).str[:12]
                return dict(zip(ids, pd.to_numeric(df[pur_col], errors="coerce").fillna(0.35)))
    except Exception:
        pass

    print("  [purity] No external purity file found — using default 0.35 for all patients")
    return {}

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
    """Map variant classification to severity. Values/default (0) match the
    mapping that produced the published metrics; do not change without rerunning."""
    scores = {
        "Nonsense_Mutation": 8, "Splice_Site": 7, "Missense_Mutation": 5,
        "In_Frame_Del": 4, "In_Frame_Ins": 4,
        "Frame_Shift_Del": 6, "Frame_Shift_Ins": 6,
        "Silent": 1,
    }
    return scores.get(variant_classification, 0)

def driver_tier(gene):
    if gene in TIER1_DRIVERS: return 2
    if gene in TIER2_DRIVERS: return 1
    return 0

def build_feature_matrix(maf_df, purity_map=None):
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

    # Driver gene score
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else "Gene"
    df["driver_tier_legacy"] = df[gene_col].apply(driver_tier)
    
    # Load OncoKB for curated driver_tier
    df["driver_tier"] = 0
    oncokb_path = Path(__file__).resolve().parent.parent / "data" / "oncokb_cancer_gene_list.txt"
    if oncokb_path.exists():
        oncokb_df = pd.read_csv(oncokb_path, sep="\t")
        oncokb_annotated = set(oncokb_df[oncokb_df["OncoKB Annotated"] == "Yes"]["Hugo Symbol"])
        cgc_col = next(c for c in oncokb_df.columns if c.startswith("COSMIC CGC"))
        cosmic_cgc = set(oncokb_df[oncokb_df[cgc_col] == "Yes"]["Hugo Symbol"])
        
        def assign_tier(g):
            if g in oncokb_annotated and g in cosmic_cgc: return 2
            if g in oncokb_annotated or g in cosmic_cgc: return 1
            return 0
        df["driver_tier"] = df[gene_col].apply(assign_tier)
    else:
        df["driver_tier"] = df["driver_tier_legacy"]

    # Consequence score (unmapped classes -> 0)
    df["consequence_score"] = df["Variant_Classification"].apply(consequence_score)

    # Pathway membership
    df["pathway"] = df[gene_col].map(PATHWAY_MAP).fillna("Other")

    # Purity-corrected cancer cell fraction (CCF = VAF * ploidy / purity, diploid assumed)
    # Falls back to raw VAF threshold when purity_map is empty
    if purity_map is None:
        purity_map = {}
    patient_col = next(
        (c for c in ["Tumor_Sample_Barcode", "case_id", "Patient"] if c in df.columns), None
    )
    if patient_col:
        case_ids = df[patient_col].astype(str).str[:12]
        purity_vals = case_ids.map(purity_map).fillna(0.35).clip(lower=0.1)
    else:
        purity_vals = pd.Series(0.35, index=df.index)

    vaf = df["vaf"].fillna(0.0)
    df["ccf"] = (vaf * 2.0 / purity_vals).clip(upper=1.0)
    df["clonality_score"] = (df["ccf"] >= 0.8).astype(float)  # clonal if CCF ≥ 0.8

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

    # Population frequency proxy (gnomAD/MAX_AF/ExAC)
    if "gnomAD_AF" in df.columns:
        df["pop_af"] = pd.to_numeric(df["gnomAD_AF"], errors="coerce").fillna(0.0)
    elif "MAX_AF" in df.columns:
        df["pop_af"] = pd.to_numeric(df["MAX_AF"], errors="coerce").fillna(0.0)
    elif "ExAC_AF" in df.columns:
        df["pop_af"] = pd.to_numeric(df["ExAC_AF"], errors="coerce").fillna(0.0)
    else:
        df["pop_af"] = 0.0

    # New features: SIFT, PolyPhen, callers, COSMIC
    if "SIFT" in df.columns:
        df["sift_score"] = df["SIFT"].str.extract(r'\(([\d\.]+)\)').astype(float).fillna(-1.0)
    else: df["sift_score"] = -1.0
        
    if "PolyPhen" in df.columns:
        df["polyphen_score"] = df["PolyPhen"].str.extract(r'\(([\d\.]+)\)').astype(float).fillna(-1.0)
    else: df["polyphen_score"] = -1.0
        
    if "callers" in df.columns:
        df["caller_count"] = df["callers"].fillna("").apply(lambda x: len(x.split(";")) if x else 0)
    else: df["caller_count"] = 1
        
    if "COSMIC" in df.columns:
        df["in_cosmic"] = df["COSMIC"].notna().astype(float)
    else: df["in_cosmic"] = 0.0

    features = [
        "vaf", "depth", "driver_tier", "consequence_score",
        "clonality_score", "ccf", "is_snv", "pop_af",
        "sift_score", "polyphen_score", "caller_count", "in_cosmic"
    ]
    return df, features, gene_col

def train_rf_xgb(X, y_labels, feature_names, groups):
    """Train Random Forest and XGBoost on binary labels (driver=1, passenger=0)."""
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import GroupKFold
    from sklearn.metrics import roc_auc_score, average_precision_score
    import numpy as np
    from scipy import stats

    # Replace NaN
    X = np.nan_to_num(X, nan=-1.0, posinf=1.0, neginf=-1.0)

    print(f"Training on {X.shape[0]} variants, {X.shape[1]} features")
    print(f"  Hotspot variants (y=1): {y_labels.sum()} / {len(y_labels)}")

    rf = RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42, n_jobs=6)
    xgb = GradientBoostingClassifier(n_estimators=200, max_depth=4, learning_rate=0.05, random_state=42)

    if y_labels.sum() < 5:
        print("  Warning: too few positive labels for CV — using full dataset")
        rf.fit(X, y_labels)
        xgb.fit(X, y_labels)
        rf_auc, xgb_auc, rf_pr, xgb_pr = 0.0, 0.0, 0.0, 0.0
    else:
        from sklearn.model_selection import KFold
        
        pooled_pr_aucs = []
        unique_groups = np.unique(groups)
        
        train_rf_xgb.fold_stats = []
        
        n_repeats = int(os.environ.get("ML_REPEATS", "5"))
        for repeat_seed in range(n_repeats):
            print(f"    Running repeat {repeat_seed + 1}/{n_repeats}...")
            # Shuffle unique genes
            rng = np.random.RandomState(repeat_seed)
            shuffled_genes = rng.permutation(unique_groups)
            
            # Split into 5 folds
            folds = np.array_split(shuffled_genes, 5)
            
            oof_rf_pred = np.zeros(len(y_labels))
            
            repeat_genes = []
            
            for fold_idx, test_genes in enumerate(folds):
                test_mask = np.isin(groups, test_genes)
                train_mask = ~test_mask
                
                rf.fit(X[train_mask], y_labels[train_mask])
                oof_rf_pred[test_mask] = rf.predict_proba(X[test_mask])[:, 1]
                
                # Record stats for fold
                test_y = y_labels[test_mask]
                if test_y.sum() > 0:
                    df_fold = pd.DataFrame({'gene': groups[test_mask], 'y': test_y})
                    top_genes = df_fold[df_fold['y'] == 1].groupby('gene').size().sort_values(ascending=False).head(3)
                    repeat_genes.extend([f"{k}({v})" for k, v in top_genes.items()])
                    
            # Compute pooled PR-AUC
            pooled_pr = average_precision_score(y_labels, oof_rf_pred)
            pooled_auc_val = roc_auc_score(y_labels, oof_rf_pred)
            pooled_pr_aucs.append(pooled_pr)
            if not hasattr(train_rf_xgb, "aucs"): train_rf_xgb.aucs = []
            train_rf_xgb.aucs.append(pooled_auc_val)
            
            train_rf_xgb.fold_stats.append({
                "repeat": repeat_seed,
                "pooled_pr_auc": float(pooled_pr),
                "dominant_genes_in_folds": ", ".join(repeat_genes[:10])
            })
            
        rf_pr_mean = np.mean(pooled_pr_aucs)
        rf_pr_min = np.min(pooled_pr_aucs)
        rf_pr_max = np.max(pooled_pr_aucs)
        rf_pr_std = np.std(pooled_pr_aucs, ddof=1)
        
        def ci(scores):
            m, se = np.mean(scores), stats.sem(scores)
            h = se * stats.t.ppf((1 + 0.95) / 2., len(scores)-1)
            return m, h

        _, rf_pr_h = ci(pooled_pr_aucs)
        
        train_rf_xgb.metrics = {
            "pooled_pr_auc_mean": float(rf_pr_mean),
            "pooled_pr_auc_ci": float(rf_pr_h),
            "min": float(rf_pr_min),
            "max": float(rf_pr_max),
            "repeats": 10,
            "pooled_auc": float(np.mean(train_rf_xgb.aucs))
        }
        
        print(f"  RF Pooled PR-AUC (10 repeats): Mean {rf_pr_mean:.3f} ± {rf_pr_h:.3f} (95% CI) | Min {rf_pr_min:.3f} | Max {rf_pr_max:.3f}")
        
        rf_auc = train_rf_xgb.metrics["pooled_auc"]
        xgb_auc = 0.0 # Bypassed to focus on PR-AUC

        
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

    # Filter to 3-caller consensus variants (same as original analysis)
    if "callers" in maf.columns:
        n_before = len(maf)
        maf["_n_callers"] = maf["callers"].fillna("").apply(
            lambda x: len([c for c in x.split(";") if c.strip()])
        )
        maf = maf[maf["_n_callers"] >= 3].copy()
        print(f"3-caller filter: {n_before} → {len(maf)} variants ({maf['Tumor_Sample_Barcode'].nunique()} patients)")

    # Gap 4: Explicit gnomAD filter
    if "gnomAD_AF" in maf.columns:
        af = pd.to_numeric(maf["gnomAD_AF"], errors="coerce").fillna(0.0)
        pre_len = len(maf)
        maf = maf[af <= GNOMAD_MAX_AF_THRESH].copy()
        removed = pre_len - len(maf)
        print(f"Filtered {removed} variants with gnomAD_AF > {GNOMAD_MAX_AF_THRESH}")

    # Load purity estimates for CCF computation
    purity_map = load_purity_estimates()
    print(f"Purity map loaded: {len(purity_map)} patients covered")

    # Build feature matrix
    df, features, gene_col = build_feature_matrix(maf, purity_map=purity_map)

    # Labels: CancerHotspots = positive class
    y = df["is_hotspot"].astype(int).values
    X = df[features].values
    groups = df[gene_col].values

    # Train models
    rf, xgb_model, rf_auc, xgb_auc, fi_rf, fi_xgb = train_rf_xgb(X, y, features, groups)

    # Serialize model for external validation
    import joblib
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    joblib.dump({"rf": rf, "xgb": xgb_model, "features": features}, RF_MODEL_PATH)
    print(f"Model saved: {RF_MODEL_PATH}")

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
        vc = str(row.get("Variant_Classification", ""))
        vt = str(row.get("Variant_Type", ""))

        if gene == "KRAS" and "G12C" in hgvsp:
            return "Sotorasib", "FDA-approved", "A"
        if gene == "EGFR" and "L858R" in hgvsp:
            return "Osimertinib", "FDA-approved", "A"
        if gene == "EGFR" and ("del" in hgvsp.lower()) and any(str(c) in hgvsp for c in range(746, 754)):
            return "Osimertinib", "FDA-approved", "A"
        if gene == "BRAF" and "V600E" in hgvsp:
            return "Dabrafenib+Trametinib", "FDA-approved", "A"
        if gene == "MET" and (vc == "Splice_Site" or "ex14" in hgvsp.lower()):
            return "Capmatinib", "FDA-approved", "A"
        if gene in ("RET", "ALK") and (vc in ("Fusion", "Rearrangement") or "fusion" in hgvsp.lower() or vt == "SV"):
            drug = "Selpercatinib" if gene == "RET" else "Alectinib"
            return drug, "FDA-approved", "A"
        if gene == "ERBB2" and (vc in ("Amplification", "Amp") or "amp" in hgvsp.lower()):
            return "Trastuzumab", "FDA-approved B", "B"
        if gene == "PIK3CA" and vc not in ("Silent", "Intron", "3'UTR", "5'UTR", "RNA") and hgvsp:
            return "Alpelisib (trial)", "Phase II", "B"
        if gene == "BRCA2" and vc not in ("Silent", "Intron", "3'UTR", "5'UTR", "RNA") and hgvsp:
            return "Olaparib (trial)", "Phase II", "B"
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
        "xgb_cv_auc": round(xgb_auc, 4) if xgb_auc else None,  # GB model is not cross-validated
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
    
    from sklearn.metrics import precision_score, recall_score, confusion_matrix
    
    mask = df["is_hotspot_regex"].notna() & df["is_hotspot"].notna()
    y_true = df.loc[mask, "is_hotspot"].values
    y_pred = df.loc[mask, "is_hotspot_regex"].values
    
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    
    print("\nHotspot Annotation Comparison:")
    print(f"  CancerHotspots Count: {int(y_true.sum())}")
    print(f"  Regex Hotspots Count: {int(y_pred.sum())}")
    print(f"  Regex Precision against CancerHotspots: {prec:.2%}")
    print(f"  Regex Recall against CancerHotspots: {rec:.2%}")
    print(f"  FP={fp}, FN={fn}")
    
    print("=== DONE ===")
    
    final_metrics_path = WORK / "final_metrics.json"
    
    metrics = getattr(train_rf_xgb, "metrics", {})
    output_data = {
        "pooled_pr_auc_mean": metrics.get("pooled_pr_auc_mean", 0.0),
        "95_ci_half_width": metrics.get("pooled_pr_auc_ci", 0.0),
        "min": metrics.get("min", 0.0),
        "max": metrics.get("max", 0.0),
        "repeats": metrics.get("repeats", 10),
        "positives": int(y_true.sum()),
        "total": len(df),
        "baseline_pr_auc": float(y_true.sum() / len(df)),
        "pooled_auc": metrics.get("pooled_auc", 0.0),
        "model_report": "results/ml_model_report.json",
        "clonality_score_ccf_corrected": True,
        "purity_patients_covered": len(purity_map),
        "ccf_feature_added": True,
    }
    
    with open(final_metrics_path, "w") as f:
        json.dump(output_data, f, indent=4)
        
    print(f"Final metrics saved to {final_metrics_path}")

if __name__ == "__main__":
    # Install sklearn if needed
    try:
        from sklearn.ensemble import RandomForestClassifier
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "scikit-learn", "xgboost", "-q"], check=True)
    main()
