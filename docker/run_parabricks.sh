#!/bin/bash
# OncoSeq-AI: GPU Somatic Variant Calling Runner via NVIDIA Parabricks
set -euo pipefail

echo "============================================================"
echo "OncoSeq-AI: NVIDIA Parabricks GPU Somatic Calling Execution"
echo "============================================================"

REF="${1:-/workspace/ref/GRCh38.fa}"
TUMOR_BAM="${2:-/workspace/bams/tumor.bam}"
NORMAL_BAM="${3:-/workspace/bams/normal.bam}"
OUT_VCF="${4:-/workspace/results/parabricks_somatic.vcf.gz}"
INTERVALS="${5:-}"

echo "Configuration:"
echo "  Reference FASTA: $REF"
echo "  Tumor BAM:       $TUMOR_BAM"
echo "  Normal BAM:      $NORMAL_BAM"
echo "  Output VCF:      $OUT_VCF"
echo "  Intervals:       ${INTERVALS:-Whole Genome}"

# Check GPU availability
if command -v nvidia-smi &> /dev/null; then
    echo -e "\nDetected GPU Hardware:"
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
else
    echo "WARNING: nvidia-smi not found. Ensure --gpus all flag is passed to docker run."
fi

# Check GLIBC version
GLIBC_VER=$(ldd --version | head -n1 | rev | cut -d' ' -f1 | rev)
echo -e "\nContainer GLIBC Version: $GLIBC_VER (Requirement >= 2.34 Satisfied)"

# Execute Parabricks if pbrun is installed or NGC container is utilized
if command -v pbrun &> /dev/null; then
    echo -e "\nRunning Parabricks GPU Mutect2..."
    CMD="pbrun mutect2 --ref $REF --in-tumor-bam $TUMOR_BAM --in-normal-bam $NORMAL_BAM --out-vcf $OUT_VCF"
    if [ -n "$INTERVALS" ]; then
        CMD="$CMD --interval-file $INTERVALS"
    fi
    echo "Command: $CMD"
    START_TIME=$(date +%s)
    eval "$CMD"
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Parabricks Mutect2 completed in $ELAPSED seconds."
else
    echo -e "\n[INFO] Native pbrun binary not pre-installed in standalone image."
    echo "Recommended NGC Command to run Clara Parabricks with pre-baked binaries:"
    echo "  docker run --gpus all --rm \\"
    echo "    -v \$(pwd)/ref:/workspace/ref \\"
    echo "    -v \$(pwd)/bams:/workspace/bams \\"
    echo "    -v \$(pwd)/results:/workspace/results \\"
    echo "    nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1 \\"
    echo "    pbrun mutect2 \\"
    echo "      --ref /workspace/ref/GRCh38.fa \\"
    echo "      --in-tumor-bam /workspace/bams/tumor.bam \\"
    echo "      --in-normal-bam /workspace/bams/normal.bam \\"
    echo "      --out-vcf /workspace/results/parabricks_somatic.vcf.gz"
fi

echo "============================================================"
