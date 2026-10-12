#!/usr/bin/env python3
"""
External validation of the ML variant prioritizer on CPTAC-LUAD.

Downloads open-access CPTAC-3 LUAD somatic MAF from GDC, applies the same
feature engineering as 04_ml_prioritizer.py, scores variants with the saved
RF model, and computes PR-AUC vs CancerHotspots labels.

Outputs:
  results/cptac_luad_maf.maf             — downloaded CPTAC-3 MAF
  results/cptac_validation_metrics.json  — PR-AUC, AUC-ROC, comparison to TCGA

Usage:
  python scripts/08_cptac_validation.py
"""
import json
import sys
import subprocess
import os
import gzip
import shutil
from pathlib import Path
from io import StringIO, BytesIO

import requests
import numpy as np
import pandas as pd

# Allow importing from this scripts directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

RESULTS    = Path(__file__).resolve().parent.parent / "results"
DATA_DIR   = Path(__file__).resolve().parent.parent / "data"
MODEL_PATH = RESULTS / "models" / "rf_model.joblib"
CPTAC_MAF  = RESULTS / "cptac_luad_maf.maf"
OUT_METRICS = RESULTS / "cptac_validation_metrics.json"

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA  = "https://api.gdc.cancer.gov/data"


# ---------- Step 1: Download CPTAC-3 LUAD MAF from GDC ----------

def fetch_cptac_maf_ids():
    """Query GDC for open-access CPTAC-3 LUAD somatic MAF file IDs."""
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "CPTAC-3"}},
            {"op": "=", "content": {"field": "cases.primary_site", "value": "Bronchus and lung"}},
            {"op": "=", "content": {"field": "data_type", "value": "Masked Somatic Mutation"}},
            {"op": "=", "content": {"field": "access", "value": "open"}},
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,file_size",
        "format": "json",
        "size": "400",
    }
    r = requests.get(GDC_FILES, params=params, timeout=30)
    r.raise_for_status()
    hits = r.json().get("data", {}).get("hits", [])
    return [h["file_id"] for h in hits]


def download_and_merge_mafs(file_ids):
    """Bulk-download MAF files and merge into a single DataFrame."""
    if not file_ids:
        raise RuntimeError("No CPTAC-3 LUAD MAF files found on GDC. Check project/disease filters.")

    print(f"  Downloading {len(file_ids)} MAF file(s) from GDC ...")
    r = requests.post(
        GDC_DATA,
        data=json.dumps({"ids": file_ids}),
        headers={"Content-Type": "application/json"},
        timeout=120,
        stream=True,
    )
    r.raise_for_status()

    content_type = r.headers.get("Content-Type", "")
    raw = b"".join(r.iter_content(chunk_size=1 << 20))

    dfs = []

    if "zip" in content_type or "tar" in content_type or len(file_ids) > 1:
        # Unpack tar/zip archive produced by bulk download
        import tarfile, zipfile
        try:
            with tarfile.open(fileobj=BytesIO(raw)) as tf:
                for member in tf.getmembers():
                    if member.name.endswith(".maf") or member.name.endswith(".maf.gz"):
                        fobj = tf.extractfile(member)
                        if fobj:
                            content = fobj.read()
                            if member.name.endswith(".gz"):
                                content = gzip.decompress(content)
                            dfs.append(_parse_maf_bytes(content))
        except tarfile.TarError:
            try:
                with zipfile.ZipFile(BytesIO(raw)) as zf:
                    for name in zf.namelist():
                        if name.endswith(".maf") or name.endswith(".maf.gz"):
                            content = zf.read(name)
                            if name.endswith(".gz"):
                                content = gzip.decompress(content)
                            dfs.append(_parse_maf_bytes(content))
            except zipfile.BadZipFile:
                # Single file returned directly
                dfs.append(_parse_maf_bytes(raw))
    else:
        dfs.append(_parse_maf_bytes(raw))

    if not dfs:
        raise RuntimeError("Could not parse any MAF files from GDC response.")

    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(CPTAC_MAF, sep="\t", index=False)
    print(f"  Saved merged MAF: {CPTAC_MAF} ({len(merged)} variants, {merged['Tumor_Sample_Barcode'].nunique()} patients)")
    return merged


def _parse_maf_bytes(raw_bytes):
    text = raw_bytes.decode("utf-8", errors="replace")
    # Skip comment lines starting with #
    lines = [l for l in text.splitlines() if not l.startswith("#")]
    return pd.read_csv(StringIO("\n".join(lines)), sep="\t", low_memory=False)


# ---------- Step 2: Apply CancerHotspots labels ----------

def apply_hotspot_labels(df, data_dir):
    """Mirror the labelling from 04_ml_prioritizer.py."""
    import re
    df["is_hotspot"] = 0.0
    hotspots_path = data_dir / "hotspots.csv"
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in df.columns else "Gene"

    if hotspots_path.exists() and "HGVSp_Short" in df.columns:
        hotspots_df = pd.read_csv(hotspots_path)
        hotspot_locs = set(zip(
            hotspots_df["Hugo_Symbol"],
            hotspots_df["Amino_Acid_Position"].astype(str)
        ))

        def check(row):
            hgvsp = str(row.get("HGVSp_Short", ""))
            gene  = str(row.get(gene_col, ""))
            if not hgvsp.startswith("p."):
                return 0.0
            m = re.search(r"p\.[A-Za-z]+(\d+)", hgvsp)
            return float((gene, m.group(1)) in hotspot_locs) if m else 0.0

        df["is_hotspot"] = df.apply(check, axis=1)
    else:
        # Fallback: regex hotspot set
        hotspot_set = {"G12C", "G12V", "G12D", "L858R", "V600E", "R175H", "R248W", "R248Q", "R273H"}
        if "HGVSp_Short" in df.columns:
            df["is_hotspot"] = df["HGVSp_Short"].fillna("").apply(
                lambda x: float(any(h in str(x) for h in hotspot_set))
            )
    return df


# ---------- Step 3: Score and evaluate ----------

def validate(maf_df, model_path, data_dir):
    import joblib
    from sklearn.metrics import roc_auc_score, average_precision_score

    if not model_path.exists():
        raise FileNotFoundError(
            f"Trained model not found at {model_path}. "
            "Run python scripts/04_ml_prioritizer.py first to train and save the model."
        )

    bundle   = joblib.load(model_path)
    rf       = bundle["rf"]
    features = bundle["features"]

    # Import build_feature_matrix from 04_ml_prioritizer
    # The module already has `if __name__ == "__main__": main()` guard so main() won't run
    import importlib.util, sys as _sys
    _mod_path = Path(__file__).resolve().parent / "04_ml_prioritizer.py"
    spec = importlib.util.spec_from_file_location("ml_prioritizer", _mod_path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    _sys.modules["ml_prioritizer"] = mod
    spec.loader.exec_module(mod)  # type: ignore[union-attr]

    df, _feat_cols, _gene_col = mod.build_feature_matrix(maf_df)

    # Align feature columns — fill any missing with 0
    for f in features:
        if f not in df.columns:
            df[f] = 0.0
    X = df[features].values
    X = np.nan_to_num(X, nan=0.0)

    # Label
    maf_df = apply_hotspot_labels(maf_df, data_dir)
    y = maf_df["is_hotspot"].astype(int).values

    if y.sum() == 0:
        print("  WARNING: No CancerHotspot positives found in CPTAC MAF — check annotation columns.")

    # Predict
    rf_prob = rf.predict_proba(X)[:, 1]

    # Metrics
    pr_auc  = average_precision_score(y, rf_prob) if y.sum() > 0 else float("nan")
    auc_roc = roc_auc_score(y, rf_prob) if y.sum() > 0 else float("nan")
    baseline = float(y.sum()) / len(y) if len(y) > 0 else 0.0
    lift = pr_auc / baseline if baseline > 0 else float("nan")

    print(f"\n  CPTAC-LUAD external validation:")
    print(f"    Variants:    {len(y)}")
    print(f"    Positives:   {int(y.sum())} ({baseline*100:.3f}% baseline)")
    print(f"    PR-AUC:      {pr_auc:.4f}")
    print(f"    AUC-ROC:     {auc_roc:.4f}")
    print(f"    Lift vs rand: {lift:.1f}×")

    return {
        "dataset": "CPTAC-3 LUAD",
        "n_variants": len(y),
        "n_positives": int(y.sum()),
        "baseline_pr_auc": round(baseline, 6),
        "pr_auc_external": round(float(pr_auc), 4),
        "auc_roc_external": round(float(auc_roc), 4),
        "lift_over_random": round(float(lift), 1),
        "tcga_internal_pr_auc": 0.258,
        "note": (
            "CPTAC-3 and TCGA share no patients. "
            "CPTAC MAF uses Illumina WGS, TCGA used WXS — coverage differences may affect recall."
        ),
    }


def main():
    print("=== CPTAC-LUAD External Validation ===")

    # Load or download MAF
    if CPTAC_MAF.exists():
        print(f"Loading cached CPTAC MAF: {CPTAC_MAF}")
        maf = pd.read_csv(CPTAC_MAF, sep="\t", comment="#", low_memory=False)
    else:
        print("Fetching CPTAC-3 LUAD MAF file IDs from GDC ...")
        file_ids = fetch_cptac_maf_ids()
        print(f"Found {len(file_ids)} file(s)")
        maf = download_and_merge_mafs(file_ids)

    print(f"MAF: {len(maf)} variants, {maf['Tumor_Sample_Barcode'].nunique()} patients")

    metrics = validate(maf, MODEL_PATH, DATA_DIR)

    with open(OUT_METRICS, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nMetrics saved: {OUT_METRICS}")

    # Compare to internal
    delta = metrics["pr_auc_external"] - metrics["tcga_internal_pr_auc"]
    direction = "higher" if delta >= 0 else "lower"
    print(f"\nSummary: external PR-AUC {metrics['pr_auc_external']:.3f} vs internal {metrics['tcga_internal_pr_auc']:.3f} "
          f"({abs(delta):.3f} {direction}) — label-overlap confound {'likely reduced' if delta < 0 else 'not inflating result'}.")


if __name__ == "__main__":
    main()
