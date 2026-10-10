# Agent Changelog — Annotation Gaps Fixes

## 1. Gap 1: ML Target Leakage (AUC = 1.0)
- **Fix:** Redefined the ML task entirely to prevent leakage. The model now predicts whether a variant is a true hotspot (label `y` = `is_hotspot` from CancerHotspots.org) using biological features. `driver_tier` and `consequence_score` are now used strictly as features, not labels. `GroupKFold` grouped by gene is maintained.
- **Before:** AUC 1.0
- **After:** RF CV AUC 0.987 ± 0.011, PR-AUC 0.268 ± 0.252. The task is highly imbalanced (248 positive examples out of 62k).

## 2. Gap 2: Regex is_hotspot
- **Fix:** Used the curated CancerHotspots DB (v2) for variant matching instead of regex.
- **Metrics against Curated DB:**
  - CancerHotspots Count: 248
  - Regex Count: 63
  - Precision: 90.48% (Regex rarely flags non-hotspots)
  - Recall: 22.98% (Regex missed 191 true hotspots!)

## 3. Gap 3: Hardcoded driver_tier
- **Fix:** Replaced hardcoded driver tiers with the official OncoKB Cancer Gene List, downloaded via their public API (no token required). Kept tier distinction: variants receive Tier 2 if the gene is in both OncoKB and COSMIC CGC, Tier 1 if only one of them, and Tier 0 otherwise.

## 4. Gap 4: Explicit gnomAD filter
- **Fix:** Filtered variants with `gnomAD_AF > 0.001`.
- **Result:** Removed exactly 3 variants. (There were exactly 3 variants strictly greater than the 0.001 threshold in this dataset).

## 5. Gap 5: Debug the hotspot column
- **Finding:** The native `hotspot` column in the MAF contains `'Y'` and `'N'` instead of boolean `TRUE`/`FALSE`. Exactly 253 variants have `'Y'`. Documented this limitation.

## 6. Gap 6: Unused annotations (SIFT/PolyPhen/COSMIC/Callers)
- **Fix:** Successfully parsed and integrated `SIFT`, `PolyPhen`, `caller_count`, and `COSMIC` into the ML features. 
- **Circularity Check:** Checked `COSMIC` for circularity against the hotspot label. While highly correlated (nearly all hotspots are in COSMIC), it is not perfectly circular because 34,384 non-hotspot variants are also in COSMIC (false positives).

## 7. Gap 7: ClinVar
- **Decision:** Documented that ClinVar was excluded from features due to extreme sparsity (0.9% coverage). It is not useful for training the main ML prioritizer.

## 8. Gap 8: Reproducibility
- **Fix:** Pinned all Python dependencies to `requirements.txt` using `pip freeze`, and set fixed `random_state=42` across ML models for determinism.
