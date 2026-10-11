# NVIDIA Parabricks GPU Deployment Guide

This directory provides the containerized environment for running NVIDIA Clara Parabricks `pbrun mutect2` GPU-accelerated somatic variant calling, resolving the GLIBC mismatch encountered on Ubuntu 20.04 systems.

---

## 1. Context & Blocker Resolution

- **The Issue:** Pre-compiled Parabricks binaries (`pbrun` v4.x) dynamically link against `GLIBC >= 2.34`. Systems running Ubuntu 20.04 (such as older AWS Deep Learning AMIs or SageMaker notebooks) provide `GLIBC 2.31`, resulting in runtime linker failure:
  ```text
  /lib/x86_64-linux-gnu/libc.so.6: version `GLIBC_2.34' not found (required by pbrun)
  ```
- **The Solution:** We deploy Parabricks inside an **Ubuntu 22.04 LTS** container which natively includes `GLIBC 2.35` and CUDA 12.2 runtime libraries with NVIDIA container toolkit integration.

---

## 2. Option A: Direct NVIDIA Clara Parabricks NGC Image (Recommended)

NVIDIA hosts pre-baked, optimized Clara Parabricks images on the NVIDIA NGC catalog.

```bash
# Pull the pre-compiled Parabricks container
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1

# Execute GPU Mutect2 somatic calling on matched Tumor-Normal BAMs
docker run --gpus all --rm \
  -v $(pwd)/ref:/workspace/ref \
  -v $(pwd)/bams:/workspace/bams \
  -v $(pwd)/results:/workspace/results \
  nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1 \
  pbrun mutect2 \
    --ref /workspace/ref/GRCh38.fa \
    --in-tumor-bam /workspace/bams/tumor.bam \
    --in-normal-bam /workspace/bams/normal.bam \
    --out-vcf /workspace/results/parabricks_somatic.vcf.gz
```

---

## 3. Option B: Custom Docker Build

To build the integrated environment containing GATK 4.5.0.0, Python utilities, and Parabricks:

```bash
# Build the Docker image
docker build -t oncoseq-parabricks:latest -f Dockerfile.parabricks .

# Run with GPU passthrough
docker run --gpus all --rm \
  -v $(pwd)/ref:/workspace/ref \
  -v $(pwd)/bams:/workspace/bams \
  -v $(pwd)/results:/workspace/results \
  oncoseq-parabricks:latest \
  /workspace/ref/GRCh38.fa \
  /workspace/bams/tumor.bam \
  /workspace/bams/normal.bam \
  /workspace/results/parabricks_somatic.vcf.gz
```

---

## 4. Hardware Verification

Inside the container, run:
```bash
nvidia-smi
ldd --version
```
Expected output:
- `nvidia-smi` detects active GPU (RTX 4060 / A10G / T4)
- `ldd (Ubuntu GLIBC 2.35-0ubuntu3.x) 2.35` confirms compatibility.
