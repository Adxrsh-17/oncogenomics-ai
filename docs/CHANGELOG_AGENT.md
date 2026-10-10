# Agent Changelog — Annotation Gaps Fixes

## 1. Gap 1: ML Target Leakage (AUC = 1.0)
- **Fix:** Dropped the calculated label. Introduced a new, strictly independent label by downloading the `CancerHotspots.org` (v2) curated database. The task is now to predict if a variant is a true known hotspot. The resulting baseline RF CV AUC is 0.987 with a PR-AUC of 0.268 (on 248 positive labels vs a random baseline PR-AUC of 0.004).
- **Ablation & Metrics:** 
  - **Full Features:** RF CV AUC 0.987, PR-AUC 0.268. The task is highly imbalanced (248 positive examples out of 62k; Baseline PR: 0.004).
  - **Without COSMIC:** CV AUC 0.983, PR-AUC 0.216 (Confirms performance isn't just relying on COSMIC as a leak).
  - **Without COSMIC & driver_tier:** CV AUC 0.808, PR-AUC 0.023 (Model fails without gene-level priors).
  - **Baselines:** Missense only gives PR-AUC 0.005. `driver_tier` only gives PR-AUC 0.092.
  - **PR-AUC Variance (±0.252):** Explained by `GroupKFold` grouping by gene. Fold 2 has 108 positives (dominated by TP53 with 79), Fold 4 has 88 positives (dominated by KRAS with 60), while Folds 0 and 3 have 9 and 10 positives respectively.
- **Limitation:** The task predicts a database label (CancerHotspots) that is available by lookup. "Hotspot" is not the same as "driver," and CancerHotspots overlaps with TCGA (where this GDC MAF likely originates).

## 2. Gap 2: Regex is_hotspot
- **Fix:** Used the curated CancerHotspots DB (v2) for variant matching instead of regex.
- **Metrics against Curated DB:**
  - CancerHotspots Count: 248
  - Regex Count: 63
  - Precision: 90.48% (Regex rarely flags non-hotspots)
  - Recall: 22.98% (Regex missed 191 true hotspots!)

## 3. Gap 3: Hardcoded driver_tier
- **Fix:** Replaced hardcoded driver tiers with the official OncoKB Cancer Gene List, downloaded via their public API (no token required). Kept tier distinction: variants receive Tier 2 if the gene is in both OncoKB and COSMIC CGC, Tier 1 if only one of them, and Tier 0 otherwise.
- **Licensing/Attribution:** OncoKB is free for academic use (requires a license for commercial). The COSMIC CGC status was extracted purely from the OncoKB file, meaning we did not bypass COSMIC's strict registration wall to download CGC independently.

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
- **Fix:** Pinned only direct Python dependencies to `requirements.txt` (`pandas`, `numpy`, `scikit-learn`, `xgboost`, `scipy`), omitting shared-environment bloat. Set fixed `random_state=42` across ML models for determinism. Conducted a clean-clone re-run to verify metrics.
