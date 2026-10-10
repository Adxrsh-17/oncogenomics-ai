#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

WORK = Path(__file__).resolve().parent.parent / "results"
MAF_FILE = WORK / "tcga_luad_merged.maf"

def main():
    import sys
    import importlib.util
    sys.path.append(str(Path(__file__).resolve().parent))
    spec = importlib.util.spec_from_file_location("ml_prioritizer", str(Path(__file__).resolve().parent / "04_ml_prioritizer.py"))
    ml_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ml_mod)
    build_feature_matrix = ml_mod.build_feature_matrix
    
    maf = pd.read_csv(MAF_FILE, sep="\t", comment="#", low_memory=False)
    af = pd.to_numeric(maf["gnomAD_AF"], errors="coerce").fillna(0.0)
    maf = maf[af <= 0.001].copy()
    
    df, features, gene_col = build_feature_matrix(maf)
    y = df["is_hotspot"].astype(int).values
    groups = df[gene_col].values
    
    cv = GroupKFold(n_splits=5)
    
    def evaluate(X, name):
        from joblib import Parallel, delayed
        unique_groups = np.unique(groups)
        
        def run_repeat(repeat_seed):
            rng = np.random.RandomState(repeat_seed)
            shuffled_genes = rng.permutation(unique_groups)
            folds = np.array_split(shuffled_genes, 5)
            
            oof_pred = np.zeros(len(y))
            rf = RandomForestClassifier(n_estimators=200, max_depth=6, random_state=42, n_jobs=1)
            
            for test_genes in folds:
                test_mask = np.isin(groups, test_genes)
                train_mask = ~test_mask
                if y[train_mask].sum() == 0 or y[test_mask].sum() == 0:
                    continue
                rf.fit(X[train_mask], y[train_mask])
                oof_pred[test_mask] = rf.predict_proba(X[test_mask])[:, 1]
                
            return average_precision_score(y, oof_pred)
            
        pooled_pr_aucs = Parallel(n_jobs=6, backend="threading")(delayed(run_repeat)(i) for i in range(10))
            
        m = np.mean(pooled_pr_aucs)
        s = np.std(pooled_pr_aucs, ddof=1)
        print(f"{name:30s} | PR-AUC: {m:.3f} ± {s:.3f} (std) | Min: {np.min(pooled_pr_aucs):.3f} | Max: {np.max(pooled_pr_aucs):.3f}")

    print(f"Random Baseline PR-AUC: {y.mean():.3f}\n")
    
    # 1. Missense only
    X_missense = (df["Variant_Classification"] == "Missense_Mutation").astype(float).values.reshape(-1, 1)
    evaluate(X_missense, "Missense only")
    
    # 2. Driver tier only
    X_tier = df[["driver_tier"]].values
    evaluate(X_tier, "driver_tier only")
    
    # 3. All except COSMIC and driver_tier
    f_no_prior = [f for f in features if f not in ("COSMIC", "in_cosmic", "driver_tier")]
    evaluate(np.nan_to_num(df[f_no_prior].values, nan=-1.0), "No COSMIC/driver_tier")
    
    # 4. All except COSMIC
    f_no_cosmic = [f for f in features if f not in ("COSMIC", "in_cosmic")]
    evaluate(np.nan_to_num(df[f_no_cosmic].values, nan=-1.0), "No COSMIC")
    
    # 5. Full Config
    evaluate(np.nan_to_num(df[features].values, nan=-1.0), "All features")

if __name__ == "__main__":
    main()
