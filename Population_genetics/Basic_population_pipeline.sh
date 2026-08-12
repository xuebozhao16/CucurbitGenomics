#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Bottle gourd population genomics analysis pipeline
# ==============================================================================
# Author: Xuebo Zhao
#
# Purpose
# -------
# This script summarizes the major analyses used for population genomic
# characterization of bottle gourd. It is intended as a reproducible workflow
# template rather than a single push-button script.
#
# Main analyses
# -------------
#  1. Variant filtering
#  2. SNP functional annotation
#  3. Phylogenetic inference
#  4. PCA and genetic-distance analysis
#  5. ADMIXTURE analysis
#  6. Linkage disequilibrium decay
#  7. Nucleotide diversity and Tajima's D
#  8. Pairwise population differentiation (FST)
#  9. IBS genetic-distance matrix
# 10. Demographic inference with SMC++
# 11. Population migration analysis with TreeMix
#
# Notes
# -----
# - Replace example file names and sample IDs with the files used in your study.
# - Population sample lists should contain one sample ID per line.
# - Chromosome names must match those in the VCF.
# - Several format-conversion scripts are project-specific and are listed as
#   external dependencies below.
# ==============================================================================


# ==============================================================================
# 0. Software and input files
# ==============================================================================

# Required command-line programs:
#
#   bcftools
#   vcftools
#   PLINK
#   ADMIXTURE
#   PopLDdecay
#   SnpEff
#   RAxML
#   Java
#   SMC++
#   TreeMix
#   bedtools
#
# Project-specific utilities:
#
#   get_depth_stats.py
#   WGS.jar
#   C23_fasta2phy.jar
#   vcf2treemix.sh
#
# Recommended: record exact software versions before publication.

# ------------------------------
# Input files
# ------------------------------

VCF_RAW="SNP_final.vcf.gz"
VCF_FILTERED="SNP_depth_filtered.vcf.gz"

# Population sample lists:
#   group_Africa.txt
#   group_America.txt
#   group_Europe.txt
#   group_South_Asia.txt
#   group_East_Asia.txt

groups=(Africa America Europe South_Asia East_Asia)

# Create analysis directories.
mkdir -p \
    01_filter \
    02_annotation \
    03_phylogeny \
    04_pca \
    05_admixture \
    06_ld \
    07_diversity \
    08_fst \
    09_ibs \
    10_smc \
    11_treemix


# ==============================================================================
# 1. VCF filtering
# ==============================================================================

# Initial site-level filtering.
#
# --max-missing 0.5
#     Retain variants genotyped in at least 50% of samples.
#
# --mac 2
#     Retain variants with a minor allele count >= 2.
#
# --minQ 30
#     Retain variants with site QUAL >= 30.
#
# --recode-INFO-all
#     Preserve INFO fields in the output VCF.

vcftools \
    --gzvcf "${VCF_RAW}" \
    --max-missing 0.5 \
    --mac 2 \
    --minQ 30 \
    --recode \
    --recode-INFO-all \
    --out 01_filter/raw_filtered


# Compress and index the filtered VCF.
bgzip -f 01_filter/raw_filtered.recode.vcf
tabix -f -p vcf 01_filter/raw_filtered.recode.vcf.gz


# ------------------------------------------------------------------------------
# 1.1 Depth-based filtering
# ------------------------------------------------------------------------------

# Calculate depth statistics for each variant.
#
# The project-specific Python script should summarize per-site depth and produce
# a file from which sites outside the desired depth range can be excluded.

python get_depth_stats.py \
    01_filter/raw_filtered.recode.vcf.gz \
    01_filter/depth_stats.txt


# "sites_depth_filtered.txt" should contain the genomic regions/sites retained
# after depth filtering in a format compatible with `bcftools view -R`.

bcftools view \
    -R 01_filter/sites_depth_filtered.txt \
    01_filter/raw_filtered.recode.vcf.gz \
    -Oz \
    -o "${VCF_FILTERED}"

tabix -f -p vcf "${VCF_FILTERED}"


# Optional quality-control summary.
bcftools stats "${VCF_FILTERED}" > 01_filter/SNP_depth_filtered.stats.txt


# ==============================================================================
# 2. SNP functional annotation with SnpEff
# ==============================================================================

# Build the custom SnpEff genome database only once.
#
# Example:
#
#   snpEff build -gff3 -v BG
#
# The genome FASTA and GFF3 annotation must be configured in the SnpEff data
# directory before building the database.

snpEff eff \
    BG \
    "${VCF_FILTERED}" \
    > 02_annotation/SNP_annotated.vcf

bgzip -f 02_annotation/SNP_annotated.vcf
tabix -f -p vcf 02_annotation/SNP_annotated.vcf.gz


# ==============================================================================
# 3. Phylogenetic analysis with RAxML
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Convert VCF to FASTA
# ------------------------------------------------------------------------------

java -jar WGS.jar \
    --model vcf \
    --type toFasta \
    --file "${VCF_FILTERED}" \
    --out 03_phylogeny/SNP.fa


# ------------------------------------------------------------------------------
# 3.2 Convert FASTA to PHYLIP
# ------------------------------------------------------------------------------

java -jar C23_fasta2phy.jar \
    --file1 03_phylogeny/SNP.fa \
    --out 03_phylogeny/SNP.phy


# ------------------------------------------------------------------------------
# 3.3 Maximum-likelihood phylogenetic inference
# ------------------------------------------------------------------------------

# -f a       rapid bootstrap analysis followed by best-tree search
# -m         nucleotide substitution model
# -p         random seed for parsimony starting tree
# -x         random seed for rapid bootstrapping
# -# 100     number of bootstrap replicates
# -T 20      number of CPU threads

cd 03_phylogeny

raxmlHPC-PTHREADS \
    -f a \
    -m GTRGAMMA \
    -p 12345 \
    -x 12345 \
    -# 100 \
    -T 20 \
    -s SNP.phy \
    -n SNP_tree

cd ..


# ==============================================================================
# 4. PCA and genetic-distance analyses
# ==============================================================================

# Convert VCF to PLINK PED/MAP format.

vcftools \
    --gzvcf "${VCF_FILTERED}" \
    --plink \
    --out 04_pca/BG


# Principal-component analysis.

plink \
    --file 04_pca/BG \
    --pca 20 header tabs \
    --out 04_pca/PCA


# Pairwise IBS-derived distance matrix.

plink \
    --file 04_pca/BG \
    --distance square 1-ibs \
    --out 04_pca/MDS


# ==============================================================================
# 5. Population structure with ADMIXTURE
# ==============================================================================

# Filter variants before ADMIXTURE analysis.
#
# --maf 0.05
#     Retain variants with MAF >= 0.05.
#
# --hwe 0.0001
#     Remove variants showing extreme deviation from Hardy-Weinberg equilibrium.

plink \
    --file 04_pca/BG \
    --maf 0.05 \
    --hwe 0.0001 \
    --make-bed \
    --out 05_admixture/ADM


# Evaluate K = 1-10 using cross-validation.

for K in {1..10}; do
    admixture \
        --cv \
        05_admixture/ADM.bed \
        "${K}" \
        | tee "05_admixture/log_K${K}.out"
done


# Extract cross-validation errors for model comparison.

grep -h "CV error" 05_admixture/log_K*.out \
    > 05_admixture/CV_error_summary.txt


# ==============================================================================
# 6. Linkage disequilibrium decay
# ==============================================================================

# Genome-wide LD decay.

PopLDdecay \
    -InVCF "${VCF_FILTERED}" \
    -OutStat 06_ld/LD_all


# Population-specific LD decay.

for group in "${groups[@]}"; do
    PopLDdecay \
        -InVCF "${VCF_FILTERED}" \
        -OutStat "06_ld/LD_${group}" \
        -SubPop "group_${group}.txt"
done


# ==============================================================================
# 7. Nucleotide diversity and Tajima's D
# ==============================================================================

# ------------------------------------------------------------------------------
# 7.1 Nucleotide diversity (pi)
# ------------------------------------------------------------------------------

# Calculate nucleotide diversity in non-overlapping 100-kb windows.

for group in "${groups[@]}"; do
    vcftools \
        --gzvcf "${VCF_FILTERED}" \
        --keep "group_${group}.txt" \
        --window-pi 100000 \
        --out "07_diversity/Pi_${group}"
done


# ------------------------------------------------------------------------------
# 7.2 Tajima's D
# ------------------------------------------------------------------------------

for group in "${groups[@]}"; do
    vcftools \
        --gzvcf "${VCF_FILTERED}" \
        --keep "group_${group}.txt" \
        --TajimaD 100000 \
        --out "07_diversity/TajimaD_${group}"
done


# ==============================================================================
# 8. Pairwise population differentiation (Weir and Cockerham FST)
# ==============================================================================

# Calculate pairwise FST for every population pair in 100-kb windows.

for ((i=0; i<${#groups[@]}; i++)); do
    for ((j=i+1; j<${#groups[@]}; j++)); do

        g1="${groups[$i]}"
        g2="${groups[$j]}"

        vcftools \
            --gzvcf "${VCF_FILTERED}" \
            --weir-fst-pop "group_${g1}.txt" \
            --weir-fst-pop "group_${g2}.txt" \
            --fst-window-size 100000 \
            --out "08_fst/Fst_${g1}_vs_${g2}"

    done
done


# ==============================================================================
# 9. IBS genetic-distance matrix
# ==============================================================================

# Convert the filtered VCF to PLINK format.

vcftools \
    --gzvcf "${VCF_FILTERED}" \
    --plink \
    --out 09_ibs/IBS


# Calculate the square pairwise distance matrix using 1 - IBS.

plink \
    --file 09_ibs/IBS \
    --distance square 1-ibs \
    --out 09_ibs/IBS_matrix


# ==============================================================================
# 10. Demographic history with SMC++
# ==============================================================================

# SMC++ requires population-specific sample definitions.
#
# Example:
#
#   Population:sample1,sample2,sample3
#
# The example below should be replaced with the actual sample IDs and chromosome
# names used in the bottle gourd VCF.


# ------------------------------------------------------------------------------
# 10.1 Convert VCF to SMC format
# ------------------------------------------------------------------------------

mkdir -p 10_smc/smc

for chr in {00..11}; do
    smc++ vcf2smc \
        "${VCF_FILTERED}" \
        "10_smc/smc/chr${chr}.smc.gz" \
        "chr${chr}" \
        "Population:sample1,sample2"
done


# ------------------------------------------------------------------------------
# 10.2 Estimate effective population size through time
# ------------------------------------------------------------------------------

# Mutation rate used in the original workflow:
#   6.96e-9 substitutions/site/generation
#
# Confirm that this value is appropriate for the final analysis before release.

smc++ estimate \
    6.96e-9 \
    10_smc/smc/chr*.smc.gz \
    --cores 10 \
    -o 10_smc/smc_out


# Optional visualization:
#
# smc++ plot 10_smc/SMC_history.pdf 10_smc/smc_out/model.final.json


# ==============================================================================
# 11. TreeMix analysis
# ==============================================================================

# TreeMix is used to model population relationships and historical migration
# edges. An outgroup should be included in the input data if the tree is rooted.


# ------------------------------------------------------------------------------
# 11.1 Prepare TreeMix input
# ------------------------------------------------------------------------------

plink \
    --vcf SNP_with_outgroup.vcf.gz \
    --make-bed \
    --out 11_treemix/TM


# Project-specific conversion from PLINK to TreeMix input.

./vcf2treemix.sh 11_treemix/TM


# ------------------------------------------------------------------------------
# 11.2 Run TreeMix
# ------------------------------------------------------------------------------

# -m 3 specifies three migration edges.
# Multiple values of m should normally be compared during exploratory analysis.

treemix \
    -i treemix.frq.gz \
    -o 11_treemix/treemix_run \
    -m 3


# ==============================================================================
# End
# ==============================================================================

echo "Bottle gourd population genomics pipeline completed."
