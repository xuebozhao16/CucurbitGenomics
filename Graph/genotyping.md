## 1. Overview
This document describes the genotyping workflow for cucumber, starting from PanGenie-compatible VCFs generated from a Minigraph-Cactus pangenome graph.

The pipeline includes:

1. Preparing the reference haplotype panel  
2. Running PanGenie on short-read samples  
3. Post-processing and filtering  
4. SV/SNP classification  
5. Genomic feature annotation  

---

## 2. Prepare Input Files for PanGenie

PanGenie requires:

- **Haplotype panel VCF (phased haploid)**  
- **Reference genome (FASTA)**  
- **Paired-end FASTQ reads**  
- **Tabix‑indexed VCF**

Prepare final panel:

```bash
bgzip -c forpangeniechrall_allhaplo_fillMiss.vcf > panel.vcf.gz
tabix panel.vcf.gz
```

---

## 3. Run PanGenie for Genotyping

Example for one sample:

```bash
PanGenie -i panel.vcf.gz  -r reference.fa -f sample_R1.fastq.gz sample_R2.fastq.gz -t 32     -o sample
```

Batch mode:

```bash
cat sample_list.txt | while read id; do
    PanGenie -i panel.vcf.gz -r reference.fa  -f fq/${id}_1.fq.gz fq/${id}_2.fq.gz  -t 32 -o ${id} &
done
wait
```

---

## 4. Merge Genotyped Samples

```bash
ls *.genotypes.vcf > vcf.list
bcftools merge -m all -Oz -o merged.genotypes.vcf.gz -l vcf.list
tabix merged.genotypes.vcf.gz
```

---

## 5. Filter Variants

Quality filter:

```bash
bcftools view -i 'QUAL>20 && INFO/AC>0' merged.genotypes.vcf.gz -Oz -o filtered.genotypes.vcf.gz
```

Missing rate filter:

```bash
bcftools view -i 'F_MISSING<0.2' filtered.genotypes.vcf.gz -Oz -o filtered2.genotypes.vcf.gz
```

---

## 6. Classify Variant Types

Insertions:

```bash
bcftools filter --include 'strlen(REF)<strlen(ALT)' filtered2.genotypes.vcf.gz > ins.vcf
```

Deletions:

```bash
bcftools filter --include 'strlen(REF)>strlen(ALT)' filtered2.genotypes.vcf.gz > del.vcf
```

SNPs:

```bash
bcftools filter --include 'strlen(REF)==1 && strlen(ALT)==1' filtered2.genotypes.vcf.gz > snp.vcf
```

---

## 7. Convert VCF to BED

```bash
vcftools --vcf ins.vcf --positions --out ins_pos
vcftools --vcf del.vcf --positions --out del_pos
vcftools --vcf snp.vcf --positions --out snp_pos
```

---

## 8. Genomic Feature Annotation

Example annotation:

```bash
bedtools intersect -a ins_pos -b gene_regions_CDS.txt -wa -u > ins_CDS.txt
bedtools intersect -a ins_pos -b upstream_regions.txt -wa -u > ins_upstream.txt
bedtools intersect -a ins_pos -b downstream_regions.txt -wa -u > ins_downstream.txt
```

Loop for all variant types:

```bash
for type in ins del snp; do
    bedtools intersect -a ${type}_pos -b gene_regions_CDS.txt -wa -u > ${type}_CDS.txt
    bedtools intersect -a ${type}_pos -b upstream_regions.txt -wa -u > ${type}_upstream.txt
    bedtools intersect -a ${type}_pos -b downstream_regions.txt -wa -u > ${type}_downstream.txt
done
```

---

## 9. Summary Counts

```bash
wc -l ins_CDS.txt ins_upstream.txt ins_downstream.txt
wc -l del_CDS.txt del_upstream.txt del_downstream.txt
wc -l snp_CDS.txt snp_upstream.txt snp_downstream.txt
```

---
