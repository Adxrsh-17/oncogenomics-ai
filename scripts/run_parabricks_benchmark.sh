#!/bin/bash
# Parabricks 4.7.1 GPU benchmark on SageMaker A10G
# Run this as: ec2-user@ip-172-28-243-126 (notebook host, not inside container)
# Prerequisites: docker pull nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1 (done)
set -euo pipefail

WORK=/home/ec2-user/SageMaker/parabricks_bench
BAM_DIR=$WORK/bams
REF_DIR=$WORK/ref
OUT_DIR=$WORK/output
LOG_DIR=$WORK/logs
mkdir -p "$BAM_DIR" "$REF_DIR" "$OUT_DIR" "$LOG_DIR"

PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG_DIR/benchmark.log"; }

log "=== Parabricks 4.7.1 Benchmark — SageMaker A10G ==="
log "Work dir: $WORK"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | tee -a "$LOG_DIR/benchmark.log"

# -------------------------------------------------------
# 1. Download HCC1143 tumor-normal BAMs (Broad public data)
#    chr17 subset, validated breast cancer cell line
# -------------------------------------------------------
TUMOR_BAM=$BAM_DIR/HCC1143_tumor_chr17.bam
NORMAL_BAM=$BAM_DIR/HCC1143_normal_chr17.bam

if [ ! -f "$TUMOR_BAM" ]; then
    log "Downloading HCC1143 tumor BAM (chr17, ~400MB)..."
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-test-data/mutect/tutorial_9183/HCC1143_tumor_chr17_dupmark.bam" \
        -O "$TUMOR_BAM" 2>&1 | tee -a "$LOG_DIR/benchmark.log" || \
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-tutorials/workshop_2002/somatic_cnv/bams/HCC1143.tumor.21.19M-29M.bam" \
        -O "$TUMOR_BAM" 2>&1 | tee -a "$LOG_DIR/benchmark.log" || \
    { log "ERROR: Could not download tumor BAM"; exit 1; }
else
    log "Tumor BAM already present: $(du -sh $TUMOR_BAM | cut -f1)"
fi

if [ ! -f "$NORMAL_BAM" ]; then
    log "Downloading HCC1143 normal BAM (chr17, ~400MB)..."
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-test-data/mutect/tutorial_9183/HCC1143_normal_chr17_dupmark.bam" \
        -O "$NORMAL_BAM" 2>&1 | tee -a "$LOG_DIR/benchmark.log" || \
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-tutorials/workshop_2002/somatic_cnv/bams/HCC1143.normal.21.19M-29M.bam" \
        -O "$NORMAL_BAM" 2>&1 | tee -a "$LOG_DIR/benchmark.log" || \
    { log "ERROR: Could not download normal BAM"; exit 1; }
else
    log "Normal BAM already present: $(du -sh $NORMAL_BAM | cut -f1)"
fi

# Index BAMs if needed
log "Indexing BAMs..."
for BAM in "$TUMOR_BAM" "$NORMAL_BAM"; do
    if [ -f "$BAM" ] && [ ! -f "${BAM}.bai" ]; then
        samtools sort -o "${BAM%.bam}_sorted.bam" "$BAM" && \
        mv "${BAM%.bam}_sorted.bam" "$BAM" && \
        samtools index "$BAM" && \
        log "Indexed: $BAM" || \
        # fallback: run samtools from inside Parabricks container
        docker run --rm -v "$BAM_DIR:/bams" "$PB_IMAGE" \
            bash -c "samtools sort -o /bams/$(basename ${BAM%.bam})_s.bam /bams/$(basename $BAM) && \
                     mv /bams/$(basename ${BAM%.bam})_s.bam /bams/$(basename $BAM) && \
                     samtools index /bams/$(basename $BAM)"
    fi
done

# -------------------------------------------------------
# 2. Download GRCh38 chr17 reference
# -------------------------------------------------------
REF_FA=$REF_DIR/hg38_chr17.fa

if [ ! -f "$REF_FA" ]; then
    log "Downloading GRCh38 chr17 reference (~90MB)..."
    wget -q --show-progress \
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz" \
        -O "${REF_FA}.gz" 2>&1 | tee -a "$LOG_DIR/benchmark.log"
    gunzip "${REF_FA}.gz"
    log "Reference downloaded: $(du -sh $REF_FA | cut -f1)"
else
    log "Reference already present: $(du -sh $REF_FA | cut -f1)"
fi

# Index reference using Parabricks container
if [ ! -f "${REF_FA}.fai" ]; then
    log "Indexing reference (bwa index + samtools fai + dict)..."
    docker run --rm -v "$REF_DIR:/ref" "$PB_IMAGE" bash -c "
        samtools faidx /ref/hg38_chr17.fa && \
        bwa index /ref/hg38_chr17.fa && \
        samtools dict /ref/hg38_chr17.fa > /ref/hg38_chr17.dict
    " 2>&1 | tee -a "$LOG_DIR/benchmark.log"
    log "Reference indexed"
fi

# -------------------------------------------------------
# 3. CPU baseline: GATK4 Mutect2 on chr17
# -------------------------------------------------------
CPU_VCF=$OUT_DIR/cpu_mutect2_chr17.vcf.gz
GATK_JAR=/home/ec2-user/SageMaker/gatk/gatk-package-4.5.0.0-local.jar

# Download GATK4 if not present
if [ ! -f "$GATK_JAR" ]; then
    log "Downloading GATK 4.5.0.0..."
    mkdir -p /home/ec2-user/SageMaker/gatk
    wget -q --show-progress \
        "https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip" \
        -O /home/ec2-user/SageMaker/gatk/gatk.zip
    cd /home/ec2-user/SageMaker/gatk && unzip -q gatk.zip
    cp /home/ec2-user/SageMaker/gatk/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar "$GATK_JAR"
    log "GATK4 installed"
fi

if [ ! -f "$CPU_VCF" ]; then
    log "--- CPU Mutect2 START ---"
    CPU_START=$(date +%s%3N)
    java -jar "$GATK_JAR" Mutect2 \
        -R "$REF_FA" \
        -I "$TUMOR_BAM" -tumor HCC1143 \
        -I "$NORMAL_BAM" -normal HCC1143BL \
        -L chr17 \
        -O "$CPU_VCF" \
        --native-pair-hmm-threads 4 \
        2>&1 | tee "$LOG_DIR/cpu_mutect2.log"
    CPU_END=$(date +%s%3N)
    CPU_MS=$(( CPU_END - CPU_START ))
    CPU_S=$(echo "scale=1; $CPU_MS/1000" | bc)
    log "CPU Mutect2 done: ${CPU_S}s"
    echo "cpu_mutect2_sec=$CPU_S" >> "$LOG_DIR/timing.txt"
else
    log "CPU VCF already exists, skipping"
    CPU_S="(cached)"
fi

# -------------------------------------------------------
# 4. GPU: Parabricks Mutect2 on chr17
# -------------------------------------------------------
GPU_VCF=$OUT_DIR/gpu_mutect2_chr17.vcf

if [ ! -f "$GPU_VCF" ]; then
    log "--- GPU Parabricks Mutect2 START ---"
    GPU_START=$(date +%s%3N)
    docker run --rm --gpus all \
        -v "$BAM_DIR:/bams" \
        -v "$REF_DIR:/ref" \
        -v "$OUT_DIR:/output" \
        -v "$LOG_DIR:/logs" \
        "$PB_IMAGE" \
        pbrun mutect2 \
            --ref /ref/hg38_chr17.fa \
            --in-tumor-bam /bams/$(basename "$TUMOR_BAM") \
            --tumor-name HCC1143 \
            --in-normal-bam /bams/$(basename "$NORMAL_BAM") \
            --normal-name HCC1143BL \
            --interval chr17 \
            --out-vcf /output/gpu_mutect2_chr17.vcf \
            --tmp-dir /tmp \
        2>&1 | tee "$LOG_DIR/gpu_mutect2.log"
    GPU_END=$(date +%s%3N)
    GPU_MS=$(( GPU_END - GPU_START ))
    GPU_S=$(echo "scale=1; $GPU_MS/1000" | bc)
    log "GPU Parabricks done: ${GPU_S}s"
    echo "gpu_mutect2_sec=$GPU_S" >> "$LOG_DIR/timing.txt"
else
    log "GPU VCF already exists, skipping"
    GPU_S="(cached)"
fi

# -------------------------------------------------------
# 5. Compare VCF concordance
# -------------------------------------------------------
log "--- VCF concordance ---"
if [ -f "$CPU_VCF" ] && [ -f "$GPU_VCF" ]; then
    CPU_COUNT=$(zcat "$CPU_VCF" 2>/dev/null | grep -v "^#" | wc -l || grep -v "^#" "$CPU_VCF" | wc -l)
    GPU_COUNT=$(grep -v "^#" "$GPU_VCF" | wc -l)
    log "CPU variants: $CPU_COUNT"
    log "GPU variants: $GPU_COUNT"
    echo "cpu_variants=$CPU_COUNT" >> "$LOG_DIR/timing.txt"
    echo "gpu_variants=$GPU_COUNT" >> "$LOG_DIR/timing.txt"
fi

# -------------------------------------------------------
# 6. GPU: Parabricks fq2bam (BWA-MEM equivalent) — measure separately
# -------------------------------------------------------
# fq2bam requires FASTQ input, skip if we only have BAM
log "Checking for FASTQ for fq2bam test..."
FASTQ_R1=$BAM_DIR/HCC1143_R1.fastq.gz
if [ ! -f "$FASTQ_R1" ]; then
    log "No FASTQ available — generating from tumor BAM (chr17 reads)..."
    docker run --rm \
        -v "$BAM_DIR:/bams" \
        "$PB_IMAGE" bash -c "
            samtools fastq -1 /bams/HCC1143_R1.fastq.gz \
                           -2 /bams/HCC1143_R2.fastq.gz \
                           /bams/$(basename "$TUMOR_BAM")
        " 2>&1 | tee -a "$LOG_DIR/benchmark.log"
    log "FASTQ generated"
fi

if [ -f "$FASTQ_R1" ]; then
    GPU_FQ2BAM_OUT=$OUT_DIR/fq2bam_out.bam
    log "--- GPU fq2bam (BWA-MEM equivalent) START ---"
    FQ_START=$(date +%s%3N)
    docker run --rm --gpus all \
        -v "$BAM_DIR:/bams" \
        -v "$REF_DIR:/ref" \
        -v "$OUT_DIR:/output" \
        "$PB_IMAGE" \
        pbrun fq2bam \
            --ref /ref/hg38_chr17.fa \
            --in-fq /bams/HCC1143_R1.fastq.gz /bams/HCC1143_R2.fastq.gz \
            --out-bam /output/fq2bam_out.bam \
            --tmp-dir /tmp \
        2>&1 | tee "$LOG_DIR/gpu_fq2bam.log"
    FQ_END=$(date +%s%3N)
    FQ_S=$(echo "scale=1; $(( FQ_END - FQ_START ))/1000" | bc)
    log "GPU fq2bam done: ${FQ_S}s"
    echo "gpu_fq2bam_sec=$FQ_S" >> "$LOG_DIR/timing.txt"
fi

# -------------------------------------------------------
# 7. Summary
# -------------------------------------------------------
log ""
log "=== BENCHMARK SUMMARY ==="
cat "$LOG_DIR/timing.txt" | tee -a "$LOG_DIR/benchmark.log"
log ""
log "Output VCFs:"
ls -lh "$OUT_DIR"/*.vcf* 2>/dev/null | tee -a "$LOG_DIR/benchmark.log"
log ""
log "Log dir: $LOG_DIR"
log "Copy results back with:"
log "  scp -P 17022 -r ec2-user@localhost:$LOG_DIR/ ./results/sagemaker/"
