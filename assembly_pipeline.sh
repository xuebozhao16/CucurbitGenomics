#!/bin/bash
#
# cucumber_assembly_pipeline.sh
#
# A clean and safe pipeline for genome assembly of Cucurbitaceae species.
# This script summarizes commonly used steps in HiFi-based assembly,
# including hifiasm assembly, duplicate removal, scaffolding, and QC.
#
# Author: Xuebo Zhao
# Repository: CucurbitGenomics
#

###############################################
# 0. Input arguments
###############################################

READS=$1          # HiFi fastq.gz file
SAMPLE=$2         # sample prefix
REF=$3            # optional reference for RagTag scaffolding

THREADS=48

if [[ -z "$READS" || -z "$SAMPLE" ]]; then
    echo "Usage: bash cucumber_assembly_pipeline.sh <reads.fastq.gz> <sample> [ref.fa]"
    exit 1
fi

###############################################
# 1. Run Hifiasm Assembly
###############################################

echo "[Step 1] Running hifiasm..."
mkdir -p 01_hifiasm

hifiasm -t ${THREADS} \
    -o 01_hifiasm/${SAMPLE} \
    ${READS}

###############################################
# 2. Extract primary contigs (from GFA)
###############################################

echo "[Step 2] Extracting contigs from GFA..."
mkdir -p 02_contigs

awk '/^S/{print ">"$2"\n"$3}' \
    01_hifiasm/${SAMPLE}.bp.p_ctg.gfa \
    > 02_contigs/${SAMPLE}.p_ctg.fa

###############################################
# 3. Purge duplicates with purge_dups
###############################################

echo "[Step 3] Running purge_dups..."
mkdir -p 03_purge_dups

# align HiFi reads
minimap2 -x map-hifi -t ${THREADS} \
    02_contigs/${SAMPLE}.p_ctg.fa ${READS} \
    | gzip -c > 03_purge_dups/${SAMPLE}.paf.gz

cd 03_purge_dups
pbcstat ${SAMPLE}.paf.gz
calcuts PB.stat > cutoffs

# self align
split_fa ../02_contigs/${SAMPLE}.p_ctg.fa > ${SAMPLE}.split.fa

minimap2 -x asm5 -DP -t ${THREADS} \
    ${SAMPLE}.split.fa ${SAMPLE}.split.fa \
    | gzip > ${SAMPLE}.self.paf.gz

purge_dups -2 -T cutoffs \
    -c PB.base.cov \
    ${SAMPLE}.self.paf.gz \
    > dups.bed

get_seqs -e dups.bed \
    ../02_contigs/${SAMPLE}.p_ctg.fa \
    > ${SAMPLE}.purged.fa
cd ..

###############################################
# 4. Optional: scaffolding with RagTag
###############################################

if [[ ! -z "$REF" ]]; then
    echo "[Step 4] Running RagTag scaffolding..."
    mkdir -p 04_ragtag

    ragtag.py scaffold \
        -t ${THREADS} \
        -o 04_ragtag \
        ${REF} \
        03_purge_dups/${SAMPLE}.purged.fa
else
    echo "[INFO] No reference genome provided; skipping RagTag."
fi

###############################################
# 5. Summary statistics
###############################################

echo "[Step 5] Generating assembly statistics..."
mkdir -p 05_stats

stats.sh in=02_contigs/${SAMPLE}.p_ctg.fa \
    > 05_stats/${SAMPLE}.raw.stats.txt

stats.sh in=03_purge_dups/${SAMPLE}.purged.fa \
    > 05_stats/${SAMPLE}.purged.stats.txt

echo "--------------------------------------------------"
echo " Assembly pipeline finished!"
echo " Results:"
echo " - Hifiasm         : 01_hifiasm/"
echo " - Contigs         : 02_contigs/${SAMPLE}.p_ctg.fa"
echo " - Purge_dups      : 03_purge_dups/${SAMPLE}.purged.fa"
echo " - RagTag (optional): 04_ragtag/"
echo " - Stats           : 05_stats/"
echo "--------------------------------------------------"
