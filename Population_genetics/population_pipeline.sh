#!/bin/bash
# ============================================================
# Bottle gourd population genomics analysis pipeline
# Author: Xuebo Zhao
# Description:
#   This script outlines the major analysis steps used for
#   bottle gourd population structure and diversity analyses.
#   It serves as a template pipeline for others to follow.
# ============================================================

# -----------------------------
# 0. Software requirements
# -----------------------------
# bcftools, vcftools, plink, PopLDdecay, snpEff, RAxML, WGS,
# smc++, treemix, bedtools


# ============================================================
# 1. VCF filtering
# ============================================================

# Basic filtering: quality, missing rate, MAC
vcftools --gzvcf SNP_final.vcf.gz \
  --max-missing 0.5 \
  --mac 2 \
  --minQ 30 \
  --recode --recode-INFO-all \
  --out raw_filtered

# Depth statistics (Python script required)
python get_depth_stats.py raw_filtered.recode.vcf depth_stats.txt

# Filter sites by depth range
bcftools view -R sites_depth_filtered.txt raw_filtered.recode.vcf.gz \
  -o SNP_depth_filtered.vcf


# ============================================================
# 2. SNP annotation (snpEff)
# ============================================================

# Build snpEff database (done only once)
# snpEff build -gff3 -v BG

# Annotate
snpEff eff BG SNP_depth_filtered.vcf > SNP_annotated.vcf


# ============================================================
# 3. Phylogenetic tree (RAxML)
# ============================================================

# Convert VCF → FASTA
java -jar WGS.jar --model vcf --type toFasta \
  --file SNP_depth_filtered.vcf.gz \
  --out SNP.fa

# Convert FASTA → PHYLIP
java -jar C23_fasta2phy.jar \
  --file1 SNP.fa \
  --out SNP.phy

# Run RAxML
raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 \
  -# 100 -T 20 \
  -s SNP.phy -n SNP_tree


# ============================================================
# 4. PCA / MDS (plink)
# ============================================================

vcftools --gzvcf SNP_depth_filtered.vcf.gz --plink --out BG
plink --file BG --pca 20 header tabs --out PCA
plink --file BG --distance square 1-ibs --out MDS


# ============================================================
# 5. ADMIXTURE
# ============================================================

plink --file BG --maf 0.05 --hwe 0.0001 --make-bed --out ADM
for k in {1..10}; do
  admixture --cv ADM.bed $k | tee log${k}.out
done


# ============================================================
# 6. Linkage Disequilibrium (PopLDdecay)
# ============================================================

PopLDdecay -InVCF SNP_depth_filtered.vcf.gz -OutStat LD_all

# Sub-populations
for group in Africa America Europe South_Asia East_Asia; do
  PopLDdecay -InVCF SNP_depth_filtered.vcf.gz \
    -OutStat LD_${group} \
    -SubPop group_${group}.txt
done


# ============================================================
# 7. π diversity & Tajima's D （vcftools）
# ============================================================

for group in Africa America Europe South_Asia East_Asia; do
  vcftools --gzvcf SNP_depth_filtered.vcf.gz \
    --keep group_${group}.txt \
    --window-pi 100000 \
    --out Pi_${group}
done

# Tajima's D
for group in Africa America Europe South_Asia East_Asia; do
  vcftools --gzvcf SNP_depth_filtered.vcf.gz \
    --keep group_${group}.txt \
    --TajimaD 100000 \
    --out TajimaD_${group}
done


# ============================================================
# 8. Fst (vcftools)
# ============================================================

groups=(Africa America Europe South_Asia East_Asia)
for i in {0..4}; do
  for j in $(seq $((i+1)) 4); do
    g1=${groups[$i]}
    g2=${groups[$j]}
    vcftools --gzvcf SNP_depth_filtered.vcf.gz \
      --weir-fst-pop group_${g1}.txt \
      --weir-fst-pop group_${g2}.txt \
      --fst-window-size 100000 \
      --out Fst_${g1}_vs_${g2}
  done
done


# ============================================================
# 9. IBS distance matrix
# ============================================================

vcftools --gzvcf SNP_depth_filtered.vcf.gz --plink --out IBS
plink --file IBS --distance square 1-ibs --out IBS_matrix


# ============================================================
# 10. Demographic history (smc++)
# ============================================================

# Convert VCF → SMC
for chr in {00..11}; do
  smc++ vcf2smc SNP_depth_filtered.vcf.gz \
    smc/chr${chr}.smc.gz \
    chr${chr} Population:sample1,sample2
done

# Estimate Ne history
smc++ estimate 6.96e-9 smc/chr*.smc.gz \
  --cores 10 -o smc_out


# ============================================================
# 11. Treemix
# ============================================================

# Convert → plink → treemix input
plink --vcf SNP_with_outgroup.vcf.gz --make-bed --out TM
./vcf2treemix.sh TM

# Run treemix
treemix -i treemix.frq.gz -o treemix_run -m 3


echo "Pipeline completed."
