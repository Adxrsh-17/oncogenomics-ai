import sys
import importlib
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score

def run_ablation():
    sys.path.append(str(Path(__file__).resolve().parent))
    ml_module = importlib.import_module("04_ml_prioritizer")
    
    WORK = Path(__file__).resolve().parent.parent / "results"
    maf_path = WORK / "tcga_luad_merged.maf"
    maf = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False)
    
    df, features, gene_col = ml_module.build_feature_matrix(maf)
    y = df["is_hotspot"].astype(int).values
    groups = df[gene_col].values
    
    # Fill NaN
    df[features] = df[features].fillna(-1.0)
    X_full = df[features].values
    
    feature_sets = {
        "Full Features": features,
        "Without COSMIC": [f for f in features if f != "in_cosmic"],
        "Without COSMIC & driver_tier": [f for f in features if f not in ("in_cosmic", "driver_tier")],
        "Only SIFT/PolyPhen/Consequence": ["sift_score", "polyphen_score", "consequence_score"]
    }
    
    baseline_pr = y.sum() / len(y)
    print(f"Total variants: {len(y)}, Positives: {y.sum()}")
    print(f"PR-AUC Baseline (Positive Rate): {baseline_pr:.4f}\n")
    
    rf = RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42, n_jobs=6)
    cv = GroupKFold(n_splits=5)
    
    for name, f_list in feature_sets.items():
        X = df[f_list].values
        aucs, prs = [], []
        
        # Per fold stats
        print(f"--- {name} ---")
        for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups)):
            if y[train_idx].sum() == 0 or y[test_idx].sum() == 0: continue
            
            rf.fit(X[train_idx], y[train_idx])
            pred = rf.predict_proba(X[test_idx])[:, 1]
            
            auc = roc_auc_score(y[test_idx], pred)
            pr = average_precision_score(y[test_idx], pred)
            aucs.append(auc)
            prs.append(pr)
            
            # Find dominating genes in this fold
            test_df = df.iloc[test_idx].copy()
            test_df["pred"] = pred
            test_df["y"] = y[test_idx]
            pos_genes = test_df[test_df["y"] == 1].groupby(gene_col).size().sort_values(ascending=False).head(3)
            genes_str = ", ".join([f"{k}({v})" for k,v in pos_genes.items()])
            
            print(f"  Fold {fold}: AUC={auc:.3f}, PR={pr:.3f} | Positives: {y[test_idx].sum()} | Top Genes: {genes_str}")
            
        print(f"  Overall {name}: CV AUC = {np.mean(aucs):.3f} ± {np.std(aucs):.3f}, PR-AUC = {np.mean(prs):.3f} ± {np.std(prs):.3f}\n")

if __name__ == "__main__":
    run_ablation()
