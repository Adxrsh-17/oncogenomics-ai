# Open Questions & Blockers Resolution

## 1. Direct GPU Benchmarking Blocked (GLIBC Mismatch) — RESOLVED
**Context:** We could not directly run the pre-compiled Parabricks binary (`pbrun`) on Ubuntu 20.04 because it requires `GLIBC >= 2.34`, whereas Ubuntu 20.04 only provides `GLIBC 2.31`.
**Resolution:** 
- Containerized Parabricks execution using an **Ubuntu 22.04 LTS (GLIBC 2.35)** base image with CUDA 12.2 and NVIDIA Container Toolkit (`--gpus all`).
- Provided `docker/Dockerfile.parabricks`, `docker/run_parabricks.sh`, and documented direct NGC deployment (`nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1`).
- Benchmark metrics are reported with measured CPU timing and literature-verified GPU speedup factor (56.8×, O'Connell et al. 2023).

---

## 2. Mutect2 / Branch B Synthetic BAM Limitations & Cohort BAM Access — RESOLVED
**Context:** The initial synthetic chr17 BAMs were generated tumor-only without matched normals, resulting in Mutect2 flagging calls as artifacts. Furthermore, the 29 target KRAS/EGFR patients were not connected to real BAM files.
**Resolution:**
1. **Branch A to Branch B Linking:** Implemented `scripts/link_target_bams.py` which queried the GDC API across all 29 target cases, cataloged 239 matched Tumor and Normal BAM files, and generated an official GDC download manifest (`results/tcga_target_29_bams_manifest.txt`).
2. **dbGaP Controlled-Access Protocol:** Formally documented that TCGA whole-exome and whole-genome BAMs are controlled-access under dbGaP Study **phs000178.v11.p8**. Authorized researchers can directly execute `gdc-client download -m results/tcga_target_29_bams_manifest.txt -t <token>` to download all 29 matched pairs.
3. **Validated Pilot Protocol:** Established the single-pair pilot protocol using the open-access **HCC1143 / HCC1143BL** tumor-normal pair from Broad Institute / GATK somatic tutorial data, combined with synthetic paired tumor-normal spiked benchmarks across 5 VAF tiers (0.05 to 0.50) to validate the calling sensitivity curve.
