# TCGA-LUAD Pipeline — Context & Progress

## Project Summary
Clinical ML pipeline for somatic variant prioritization in lung adenocarcinoma (LUAD), using TCGA-LUAD as the discovery cohort. End goal: identify and visualize actionable mutations for targeted therapy matching.

---

## Cohort
| Metric | Value |
|---|---|
| Patients | 198 (post-QC) |
| Somatic variants | 62,146 |
| Actionable variants | 29 |
| KRAS G12C → Sotorasib | 21 patients |
| EGFR L858R → Osimertinib | 8 patients |

---

## Pipeline Stages & Status

### ✅ 1. MAF Download (`scripts/01_download_full_maf.py`)
- Downloads open-access TCGA-LUAD masked somatic mutation MAF from GDC bulk API
- Output: `results/tcga_luad_merged.maf` — 618 files, 194,729 variants, 616 patients
- Note: script hardcodes SageMaker path `/opt/ml/...`; local run requires inline Python workaround

### ✅ 2. ML Prioritizer (`scripts/04_ml_prioritizer.py`)
- Random Forest + XGBoost ensemble scoring variants against CancerHotspots v2 labels
- **Purity-corrected CCF** added: `ccf = min(1.0, VAF × 2.0 / purity)`, clonal if CCF ≥ 0.8
- 12-feature model: `vaf, depth, driver_tier, consequence_score, clonality_score, ccf, is_snv, pop_af, sift_score, polyphen_score, caller_count, in_cosmic`
- Model serialized to `results/models/rf_model.joblib` for external validation reuse
- Internal PR-AUC: **0.258** (honest, no label-overlap inflation)
- Hotspots downloaded from CancerHotspots v2 API → `data/hotspots.csv` (1,165 positions)
- Outputs: `results/final_metrics.json`, `results/ml_model_report.json`

### ✅ 3. GDC Normal BAM Manifest (`scripts/07_gdc_query_normals.py`)
- Queries GDC API for matched tumor/normal BAM IDs for the 29 actionable patients
- Output: `results/gdc_normal_manifest.txt` (case_id, tumor_bam_id, normal_bam_id, sizes)
- Controlled-access download requires GDC dbGaP token: `gdc-client download -m manifest.txt -t token.txt`

### ✅ 4. CPTAC External Validation (`scripts/08_cptac_validation.py`)
- Downloads CPTAC-3 LUAD open-access MAF from GDC
- Imports `build_feature_matrix` from `04_ml_prioritizer.py` (no duplication)
- Loads saved RF model, scores variants, computes PR-AUC vs CancerHotspots labels
- Output: `results/cptac_validation_metrics.json`

### ✅ 5. 3D Mutation Visualization (`viz3d/`)
- React + Vite app with 3Dmol.js for interactive protein structure viewing
- **GPU-predicted structures** via NVIDIA BioNeMo OpenFold2 hosted inference (API key secured in `~/.config/bionemo/credentials`)
  - KRAS G12C: `results/kras_g12c.pdb` (171 aa, OpenFold2)
  - EGFR L858R: `results/egfr_l858r.pdb` (169 aa kinase fragment, OpenFold2)
- **Crystal + drug structures** from RCSB PDB:
  - KRAS G12C + Sotorasib: `results/kras_g12c_crystal.pdb` (PDB 6OIM)
  - EGFR L858R + Osimertinib: `results/egfr_l858r_crystal.pdb` (PDB 4ZAU)
- Features: protein selector, GPU vs crystal toggle, cartoon/surface/ball-stick styles, mutation site highlight, drug molecule toggle, auto-rotate
- Dev server: `cd viz3d && npm run dev` → `http://localhost:5173` (also `http://10.24.39.19:5173` on network)

---

## Key File Locations

```
dna/
├── scripts/
│   ├── 04_ml_prioritizer.py     # ML model with purity-corrected CCF
│   ├── 07_gdc_query_normals.py  # GDC matched normal BAM manifest
│   └── 08_cptac_validation.py   # External validation on CPTAC-LUAD
├── results/
│   ├── tcga_luad_merged.maf     # 194,729 somatic variants
│   ├── patient_clinical_report.csv
│   ├── final_metrics.json
│   ├── gdc_normal_manifest.txt  # 29 actionable patients' BAM IDs
│   ├── kras_g12c.pdb            # GPU-predicted (OpenFold2)
│   ├── egfr_l858r.pdb           # GPU-predicted (OpenFold2)
│   ├── kras_g12c_crystal.pdb    # PDB 6OIM (crystal + Sotorasib)
│   ├── egfr_l858r_crystal.pdb   # PDB 4ZAU (crystal + Osimertinib)
│   └── models/
│       └── rf_model.joblib      # Serialized RF model for CPTAC reuse
├── data/
│   └── hotspots.csv             # CancerHotspots v2 (1,165 positions)
└── viz3d/                       # React 3D visualization app
    ├── src/
    │   ├── App.jsx
    │   ├── Viewer3D.jsx
    │   └── pdbData.js           # All 4 PDB structures as JS literals (518KB)
    └── vite.config.js           # Bound to 0.0.0.0:5173 for network access
```

---

## Pending / Next Steps
- **Run `04_ml_prioritizer.py`** to produce updated PR-AUC with CCF feature and serialize model
- **Run `08_cptac_validation.py`** after model is serialized
- **Full fq2bam** (separate GPU branch): requires 2× A10G, Parabricks 4.7.1, GDC BAM access
- **GPU upgrade**: Parabricks `pbrun mutectcaller` on raw BAMs for the 29 actionable patients once GDC token obtained

---

## Key Technical Notes
- CancerHotspots v2 API: `https://www.cancerhotspots.org/api/hotspots/single`
- GDC open-access MAF: `https://api.gdc.cancer.gov/files` + `/data` bulk POST
- NVIDIA OpenFold2 function ID: `fb8a1686-a8c7-4302-a898-b9e150d9f0f6` (ACTIVE as of Oct 2026)
- 3Dmol.js requires real WebGL context — sandboxed iframes (Claude artifacts) only show spheres; use the Vite dev server
- EGFR full kinase (495 aa) OOM-kills OpenFold2 inference; use 169 aa kinase domain fragment centered on L858R
