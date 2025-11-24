## 1. Overview
This document describes the complete workflow for constructing a pangenome graph for cucumber using **Minigraph-Cactus**, followed by **variant genotyping** and preparation of **VCF files suitable for PanGenie**.

---

## 2. Prepare Assemblies by Chromosome
```bash
cd genomes36
cat ../taxa_list.txt | while read line; do
    for i in {1..7}; do
        samtools faidx /assembly/${line}.fa chr${i} | seqtk seq -l0 > ${line}.chr${i}.fa
    done
done
```

---

## 3. Create seqFile for Minigraph-Cactus
Example for chr1:
```
AM716   AM716.chr1.fa
AM117   AM117.chr1.fa
...
```

Generate seqFiles for chr2–chr7:
```bash
for i in {2..7}; do
    sed "s/chr1/chr${i}/g" chr1.evolver.txt > chr${i}.evolver.txt
done
```

---

## 4. Run Minigraph-Cactus per Chromosome
```bash
cactus-pangenome ./js.chr1 ./chr1.evolver.txt  --reference AM716 --mapCores 90 --mgCores 64   --maxLen 10000 --clip 10000 --filter 7 --permissiveContigFilter 0.1  --outDir chr1.out --outName chr1  --vcf --giraffe clip filter --gfa clip full --gbz clip filter full --odgi --logFile chr1.mc.log &
```

Batch run chr2–chr7:
```bash
for i in {2..7}; do
    nohup cactus-pangenome ./js.chr${i} ./chr${i}.evolver.txt  --reference AM716    --mapCores 20 --mgCores 64  --maxLen 10000 --clip 10000  --filter 7 --permissiveContigFilter 0.1 --outDir chr${i}.out --outName chr${i}  --vcf --giraffe clip filter --gfa clip full --gbz clip filter full  --odgi --logFile chr${i}.mc.log &
done
```

---

## 5. Merge VCF Files Across Chromosomes
```bash
cd result
for i in {1..7}; do
    cp ../chr${i}.out/chr${i}.vcf.gz .
done

vcf-concat chr1.vcf.gz chr2.vcf.gz ... chr7.vcf.gz | bgzip -c > chr_all.vcf.gz
tabix chr_all.vcf.gz
```

---

## 6. Remove Nested Variants
```bash
for i in {1..7}; do
    vcfbub -l 0 -r 100000 --input chr${i}.vcf.gz > vcfbub.chr${i}.vcf
done
```

---

## 7. Normalize VCF (split multi-allelic)
```bash
for i in {1..7}; do
    bcftools norm -m -any vcfbub.chr${i}.vcf > normalized_bi_chr${i}.vcf
done
```

---

## 8. Convert Genotypes to Phased 0|0 and 1|1
```bash
for i in {1..7}; do
    grep -v "#" normalized_bi_chr${i}.vcf |
    awk '{
        for(i=10;i<=NF;i++){
            if($i=="0")$i="0|0";
            else if($i=="1")$i="1|1";
            else if($i==".")$i=".|.";
        }}1' > chr${i}.modified_file.vcf

    grep "#" normalized_bi_chr${i}.vcf > header${i}.txt
    cat header${i}.txt chr${i}.modified_file.vcf > chr${i}.modified_file2.vcf
done
```

---

## 9. Fill Missing Genotypes and Merge Overlapping Variants
```bash
for i in {1..7}; do
    perl fillmiss_VCF_v1.pl $0 chr${i}.modified_file2.vcf > chr${i}_fillMiss.vcf
    bgzip -c chr${i}_fillMiss.vcf > chr${i}_fillMiss.vcf.gz
    bcftools index chr${i}_fillMiss.vcf.gz
done

for i in {1..7}; do
    bcftools view chr${i}_fillMiss.vcf.gz  | perl update_svlen.pl | perl merge_ovl_var_forPanGenie.pl AM716.chr${i}.fa  | bcftools view -Oz -o v1chr${i}_fillMiss.vcf.gz
done
```

---

## 10. Convert to Haploid Genotype (PanGenie)
```bash
for i in {1..7}; do
    bcftools annotate -x INFO v1chr${i}_fillMiss.vcf.gz | perl merge_diplo2haplo.pl > forpangeniechr${i}_fillMiss.vcf
done
```

---

## 11. Merge Final VCF Files
```bash
vcf-concat v1chr1_fillMiss.vcf.gz ... v1chr7_fillMiss.vcf.gz > forpangeniechrall_allhaplo_fillMiss.vcf
```

---

## 12. SV Classification
```bash
bcftools filter --include 'strlen(REF)<strlen(ALT)'   pangenie_index_filter_bi.vcf > ins
bcftools filter --include 'strlen(REF)>strlen(ALT)'   pangenie_index_filter_bi.vcf > del
bcftools filter --include 'strlen(REF)==1 && strlen(ALT)==1' pangenie_index_filter_bi.vcf > snp
```
---
