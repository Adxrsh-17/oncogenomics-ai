#!/usr/bin/env python3
"""
Query GDC for matched blood-normal BAM file IDs for the 29 actionable TCGA-LUAD patients.

Outputs:
  results/gdc_normal_manifest.txt  — TSV: case_id, tumor_bam_id, normal_bam_id, size_gb

IMPORTANT: BAMs are GDC controlled-access (dbGaP phs000178).
To download them you need a GDC token. Once you have one:

  gdc-client download -m results/gdc_normal_manifest.txt -t ~/gdc-user-token.txt -n 4

This script only builds the manifest — no controlled data is fetched here.
"""
import json
import sys
import subprocess
from pathlib import Path

import requests
import pandas as pd

RESULTS = Path(__file__).resolve().parent.parent / "results"
REPORT  = RESULTS / "patient_clinical_report.csv"
OUT     = RESULTS / "gdc_normal_manifest.txt"

GDC_FILES  = "https://api.gdc.cancer.gov/files"
GDC_CASES  = "https://api.gdc.cancer.gov/cases"

# Only actionable patients
def load_target_cases():
    df = pd.read_csv(REPORT)
    # Expect columns: Patient (TCGA barcode), KRAS_G12C, EGFR_L858R, Actionable
    actionable_col = next((c for c in df.columns if "action" in c.lower()), None)
    if actionable_col:
        targets = df[df[actionable_col].notna() & (df[actionable_col] != "")]
    else:
        # Fall back: patients with KRAS G12C or EGFR L858R
        kras_col = next((c for c in df.columns if "kras" in c.lower()), None)
        egfr_col = next((c for c in df.columns if "egfr" in c.lower()), None)
        mask = pd.Series([False] * len(df))
        if kras_col:
            mask = mask | df[kras_col].astype(str).str.contains("G12C|True|1", na=False)
        if egfr_col:
            mask = mask | df[egfr_col].astype(str).str.contains("L858R|True|1", na=False)
        targets = df[mask]

    patient_col = next((c for c in df.columns if "patient" in c.lower() or c == "Patient"), "Patient")
    # Normalise to 12-char TCGA barcode
    return list({str(x)[:12] for x in targets[patient_col].tolist()})


def query_bams_for_case(case_id):
    """Return list of {file_id, file_name, data_category, experimental_strategy, file_size}."""
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.submitter_id", "value": case_id}},
            {"op": "=", "content": {"field": "data_format", "value": "BAM"}},
            {"op": "in", "content": {
                "field": "experimental_strategy",
                "value": ["WGS", "WXS"]
            }},
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,data_category,experimental_strategy,file_size,cases.samples.sample_type",
        "format": "json",
        "size": "20",
    }
    r = requests.get(GDC_FILES, params=params, timeout=30)
    r.raise_for_status()
    return r.json().get("data", {}).get("hits", [])


def classify_sample(hit):
    """Return 'tumor' or 'normal' based on sample_type."""
    try:
        sample_type = hit["cases"][0]["samples"][0]["sample_type"].lower()
    except (KeyError, IndexError):
        # Guess from filename
        name = hit.get("file_name", "").lower()
        sample_type = "normal" if "normal" in name or "_bl" in name else "tumor"
    if "normal" in sample_type or "blood" in sample_type:
        return "normal"
    return "tumor"


def main():
    print("=== GDC Matched Normal BAM Query ===")

    if not REPORT.exists():
        print(f"ERROR: {REPORT} not found. Run 05_full_branch_a.py first.")
        sys.exit(1)

    cases = load_target_cases()
    print(f"Target cases: {len(cases)}")

    rows = []
    missing = []

    for i, case_id in enumerate(cases, 1):
        print(f"  [{i}/{len(cases)}] {case_id} ...", end=" ", flush=True)
        try:
            hits = query_bams_for_case(case_id)
        except Exception as e:
            print(f"ERROR: {e}")
            missing.append(case_id)
            continue

        tumor_hits  = [h for h in hits if classify_sample(h) == "tumor"]
        normal_hits = [h for h in hits if classify_sample(h) == "normal"]

        tumor_id  = tumor_hits[0]["file_id"]  if tumor_hits  else ""
        normal_id = normal_hits[0]["file_id"] if normal_hits else ""
        tumor_size  = round(tumor_hits[0].get("file_size", 0) / 1e9, 2)  if tumor_hits  else 0
        normal_size = round(normal_hits[0].get("file_size", 0) / 1e9, 2) if normal_hits else 0

        status = "complete" if tumor_id and normal_id else ("tumor_only" if tumor_id else "missing")
        print(status)

        rows.append({
            "case_id": case_id,
            "tumor_bam_id":  tumor_id,
            "tumor_size_gb": tumor_size,
            "normal_bam_id": normal_id,
            "normal_size_gb": normal_size,
            "status": status,
        })

    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT, sep="\t", index=False)

    complete = (manifest["status"] == "complete").sum()
    print(f"\nManifest saved: {OUT}")
    print(f"  Complete pairs: {complete}/{len(cases)}")
    print(f"  Missing normals: {(manifest['status'] != 'complete').sum()}")

    total_gb = manifest["tumor_size_gb"].sum() + manifest["normal_size_gb"].sum()
    print(f"  Estimated download: {total_gb:.1f} GB (controlled access)")

    if missing:
        print(f"\nFailed queries ({len(missing)}): {', '.join(missing)}")

    print("""
To download (requires GDC token):
  pip install gdc-client
  gdc-client download -m results/gdc_normal_manifest.txt -t ~/gdc-user-token.txt -n 4

Get your token at: https://portal.gdc.cancer.gov (Login → Download Token)
""")


if __name__ == "__main__":
    try:
        import requests
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "requests", "-q"], check=True)
        import requests
    main()
