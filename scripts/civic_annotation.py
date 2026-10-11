#!/usr/bin/env python3
"""
CIViC Clinical Variant Interpretation Engine
Queries the public CIViC (Clinical Interpretation of Variants in Cancer) GraphQL API
to retrieve evidence-linked clinical assertions, drug associations, clinical trials,
and evidence levels (A-D) for TCGA-LUAD target variants.

Outputs:
  - results/civic_clinical_evidence.json
  - results/civic_clinical_evidence.csv
"""

import json
import os
import sys
import time
import urllib.request
from pathlib import Path
import pandas as pd

CIVIC_GRAPHQL_URL = "https://civicdb.org/api/graphql"
OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Key actionable / prominent variants in TCGA-LUAD
TARGET_VARIANTS = [
    {"gene": "KRAS", "variant": "G12C", "profile": "KRAS G12C"},
    {"gene": "KRAS", "variant": "G12D", "profile": "KRAS G12D"},
    {"gene": "KRAS", "variant": "G12V", "profile": "KRAS G12V"},
    {"gene": "EGFR", "variant": "L858R", "profile": "EGFR L858R"},
    {"gene": "EGFR", "variant": "T790M", "profile": "EGFR T790M"},
    {"gene": "BRAF", "variant": "V600E", "profile": "BRAF V600E"},
    {"gene": "PIK3CA", "variant": "E545K", "profile": "PIK3CA E545K"},
    {"gene": "PIK3CA", "variant": "H1047R", "profile": "PIK3CA H1047R"},
    {"gene": "TP53", "variant": "R175H", "profile": "TP53 R175H"},
    {"gene": "TP53", "variant": "R273H", "profile": "TP53 R273H"},
    {"gene": "ERBB2", "variant": "V777L", "profile": "ERBB2 V777L"},
    {"gene": "MET", "variant": "Exon 14 Skipping", "profile": "MET Exon 14 Skipping"},
]

GRAPHQL_QUERY = """
query GetMolecularProfileEvidence($name: String!) {
  molecularProfiles(name: $name, first: 1) {
    nodes {
      id
      name
      link
      evidenceItems(first: 20) {
        nodes {
          id
          name
          evidenceType
          evidenceLevel
          evidenceDirection
          significance
          therapies {
            name
          }
          disease {
            name
          }
          source {
            citation
            sourceUrl
            publicationYear
          }
        }
      }
    }
  }
}
"""

def query_civic(profile_name):
    """Query CIViC GraphQL API for a molecular profile."""
    payload = json.dumps({
        "query": GRAPHQL_QUERY,
        "variables": {"name": profile_name}
    }).encode("utf-8")

    req = urllib.request.Request(
        CIVIC_GRAPHQL_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "OncoSeq-AI-Pipeline/1.0"
        }
    )

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            nodes = data.get("data", {}).get("molecularProfiles", {}).get("nodes", [])
            if nodes:
                return nodes[0]
    except Exception as e:
        print(f"  [WARN] Failed query for {profile_name}: {e}")
    return None

def main():
    print("=" * 60)
    print("CIViC Evidence Annotation — OncoSeq-AI Clinical Interpretation")
    print(f"Querying CIViC API ({CIVIC_GRAPHQL_URL}) for {len(TARGET_VARIANTS)} variants...")
    print("=" * 60)

    all_evidence = []
    summary_records = []

    for item in TARGET_VARIANTS:
        pname = item["profile"]
        print(f"Fetching CIViC profile: {pname}...")
        profile = query_civic(pname)
        time.sleep(0.3)  # Respectful rate limiting

        if not profile:
            print(f"  -> No direct hit for {pname}")
            continue

        e_items = profile.get("evidenceItems", {}).get("nodes", [])
        print(f"  -> Found {len(e_items)} clinical evidence items")

        for ev in e_items:
            therapies = ", ".join([t["name"] for t in ev.get("therapies", [])]) or "None"
            disease = ev.get("disease", {}).get("name") if ev.get("disease") else "Not specified"
            src = ev.get("source") or {}
            citation = src.get("citation", "")
            src_url = src.get("sourceUrl", "")

            record = {
                "gene": item["gene"],
                "variant": item["variant"],
                "molecular_profile": profile["name"],
                "civic_eid": ev["name"],
                "evidence_type": ev.get("evidenceType"),
                "evidence_level": ev.get("evidenceLevel"),
                "clinical_significance": ev.get("significance"),
                "evidence_direction": ev.get("evidenceDirection"),
                "therapies": therapies,
                "disease": disease,
                "citation": citation,
                "source_url": src_url,
                "civic_url": f"https://civicdb.org{profile.get('link', '')}"
            }
            all_evidence.append(record)

        # High-level summary per variant
        top_therapies = set()
        levels = set()
        for ev in e_items:
            for t in ev.get("therapies", []):
                top_therapies.add(t["name"])
            if ev.get("evidenceLevel"):
                levels.add(ev.get("evidenceLevel"))

        summary_records.append({
            "gene": item["gene"],
            "variant": item["variant"],
            "profile": pname,
            "evidence_item_count": len(e_items),
            "evidence_levels": "/".join(sorted(levels)) if levels else "N/A",
            "approved_therapies": "; ".join(sorted(top_therapies)) if top_therapies else "N/A",
            "civic_profile_id": profile["id"]
        })

    # Save detailed evidence table
    ev_df = pd.DataFrame(all_evidence)
    ev_csv = OUT_DIR / "civic_clinical_evidence.csv"
    ev_df.to_csv(ev_csv, index=False)

    ev_json = OUT_DIR / "civic_clinical_evidence.json"
    with open(ev_json, "w") as f:
        json.dump(all_evidence, f, indent=2)

    sum_df = pd.DataFrame(summary_records)
    sum_csv = OUT_DIR / "civic_actionable_summary.csv"
    sum_df.to_csv(sum_csv, index=False)

    print("\n" + "=" * 60)
    print(f"CIViC Clinical Annotation Complete:")
    print(f"  - Total clinical evidence items retrieved: {len(all_evidence)}")
    print(f"  - Target profiles annotated: {len(summary_records)}")
    print(f"  - Output CSV: {ev_csv}")
    print(f"  - Output JSON: {ev_json}")
    print(f"  - Actionable Summary CSV: {sum_csv}")
    print("=" * 60)
    print("\nSample Actionable Evidence:")
    print(sum_df[["gene", "variant", "evidence_item_count", "evidence_levels", "approved_therapies"]].to_string(index=False))

if __name__ == "__main__":
    main()
