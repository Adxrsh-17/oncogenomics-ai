#!/bin/bash
# Setup GATK4 + try pip Parabricks + download real tumor-normal BAMs
# Working dir: /opt/ml/input/data/notebook/oncoseq/
set -euo pipefail

WORK=/opt/ml/input/data/notebook/oncoseq
GATK_DIR=$WORK/gatk
TOOLS=$WORK/tools
BAM_DIR=$WORK/bams
REF_DIR=$WORK/ref
mkdir -p $GATK_DIR $TOOLS $BAM_DIR $REF_DIR

log() { echo "[$(date +%H:%M:%S)] $*"; }

# -------------------------------------------------------
# 1. Download GATK4 jar if not present
# -------------------------------------------------------
GATK_JAR=$GATK_DIR/gatk.jar
if [ ! -f "$GATK_JAR" ]; then
    log "Downloading GATK4 4.5.0.0..."
    wget -q --show-progress \
        "https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip" \
        -O $GATK_DIR/gatk.zip
    cd $GATK_DIR && unzip -q gatk.zip
    cp $GATK_DIR/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar $GATK_JAR
    log "GATK4 installed: $GATK_JAR"
else
    log "GATK4 already present"
fi

# -------------------------------------------------------
# 2. Try nvidia-clara-parabricks pip install
# -------------------------------------------------------
log "Trying pip install nvidia-clara-parabricks..."
pip3 install nvidia-clara-parabricks 2>&1 | tail -10 || log "pip parabricks failed (expected if not on NGC)"
which pbrun 2>/dev/null && log "pbrun available: $(which pbrun)" || log "pbrun not in PATH after pip install"

# -------------------------------------------------------
# 3. Download HCC1143 tumor-normal BAMs (public, Broad GATK tutorial data)
#    These are chromosome-subset BAMs (chr17, ~1GB each) from a validated
#    breast cancer cell line — standard for somatic calling benchmarks
# -------------------------------------------------------
log "Downloading HCC1143 tumor BAM (chr17 subset)..."

TUMOR_BAM=$BAM_DIR/HCC1143_tumor.bam
NORMAL_BAM=$BAM_DIR/HCC1143_normal.bam

if [ ! -f "$TUMOR_BAM" ]; then
    # Broad public bucket — GATK4 somatic tutorial data
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-test-data/mutect/tutorial_9183/HCC1143_normal_chr17_dupmark.bam" \
        -O $NORMAL_BAM 2>/dev/null || \
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-tutorials/workshop_2002/somatic_cnv/bams/HCC1143_small_normal.bam" \
        -O $NORMAL_BAM 2>/dev/null || \
    log "WARNING: Could not download normal BAM from primary URL, trying backup..."

    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-test-data/mutect/tutorial_9183/HCC1143_tumor_chr17_dupmark.bam" \
        -O $TUMOR_BAM 2>/dev/null || \
    wget -q --show-progress \
        "https://storage.googleapis.com/gatk-tutorials/workshop_2002/somatic_cnv/bams/HCC1143_small.bam" \
        -O $TUMOR_BAM 2>/dev/null || \
    log "WARNING: Could not download tumor BAM"
fi

# BAI indexes
log "Indexing BAMs (if files downloaded)..."
for BAM in $TUMOR_BAM $NORMAL_BAM; do
    if [ -f "$BAM" ] && [ ! -f "${BAM}.bai" ]; then
        python3 -c "import pysam; pysam.sort('-o', '${BAM%.bam}_sorted.bam', '$BAM'); pysam.index('${BAM%.bam}_sorted.bam')" 2>/dev/null || \
        java -jar $GATK_JAR SortSam -I $BAM -O ${BAM%.bam}_sorted.bam -SORT_ORDER coordinate && \
        java -jar $GATK_JAR BuildBamIndex -I ${BAM%.bam}_sorted.bam || true
    fi
done

# -------------------------------------------------------
# 4. Download GRCh38 chr17 reference (needed for HCC1143)
#    and chr12+chr7 for KRAS/EGFR targets
# -------------------------------------------------------
log "Downloading GRCh38 chr17 reference sequence..."
REF_FA=$REF_DIR/GRCh38_chr17.fa
if [ ! -f "$REF_FA" ]; then
    wget -q --show-progress \
        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz" \
        -O ${REF_FA}.gz 2>/dev/null
    if [ -f "${REF_FA}.gz" ]; then
        gunzip ${REF_FA}.gz
        log "chr17 reference: $(du -sh $REF_FA | cut -f1)"
    else
        log "WARNING: chr17 reference download failed"
    fi
fi

# Index reference
if [ -f "$REF_FA" ] && [ ! -f "${REF_FA}.fai" ]; then
    python3 -c "import pysam; pysam.faidx('$REF_FA')"
    java -jar $GATK_JAR CreateSequenceDictionary -R $REF_FA || true
    log "Reference indexed"
fi

log "=== Setup DONE ==="
ls -lh $BAM_DIR/ 2>/dev/null
ls -lh $REF_DIR/ 2>/dev/null
