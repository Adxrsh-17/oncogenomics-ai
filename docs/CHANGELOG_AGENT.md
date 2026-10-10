# Agent Changelog — Annotation Gaps Fixes

## 1. Gap 1: ML Target Leakage (AUC = 1.0)
- **Fix:** Removed `driver_tier` and `consequence_score` from the ML features since they directly defined the positive class. Replaced `StratifiedKFold` with `GroupKFold` grouped by gene so the model evaluates generalizability across different genes.
- **Before:** AUC 1.0
- **After:** AUC ~0.55 (Expected as VAF/depth alone aren't highly predictive without driver/pathway info)

## 2. Gap 2: Regex is_hotspot
- **Fix:** Downloaded the CancerHotspots DB (v2) and converted it to CSV. Rewrote the logic in `04_ml_prioritizer.py` to check `Hugo_Symbol` and the parsed genomic position from `HGVSp_Short` against the curated DB.
- **Result:** ~99.68% agreement between the old regex approach and the new curated CancerHotspots DB.

## 3. Gap 3: Hardcoded driver_tier
- **Fix:** Downloaded the OncoKB Cancer Gene List via their API. Refactored `04_ml_prioritizer.py` to assign a tier of 2 to variants in OncoKB annotated genes, storing the previous hardcoded lists under `driver_tier_legacy`.

## 4. Gap 4: Explicit gnomAD filter
- **Fix:** Added a `GNOMAD_MAX_AF_THRESH = 0.001` configuration variable to both `04_ml_prioritizer.py` and `05_full_branch_a.py`.
- **Result:** Removed 3 variants with `gnomAD_AF > 0.001` during the filtering step.

## 5. Gap 5: Debug the hotspot column
- **Finding:** The native `hotspot` column in the MAF uses `'Y'` and `'N'` instead of boolean `TRUE`/`FALSE`. Documented this finding in the README.

## 6. Gap 6: Unused annotations (SIFT/PolyPhen)
- **Finding:** Wrote and ran `scripts/check_sift.py`. Validated that SIFT and PolyPhen annotations from VEP are present in ~64.3% and 64.7% of the variants, respectively. Documented this limitation in the README.
