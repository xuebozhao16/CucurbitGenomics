# GERP Conservation Analysis in Cucumber

## Overview

This workflow describes the calculation of genomic evolutionary constraint in cucumber using GERP++. The cucumber accession **AM716** is used as the reference genome. Whole-genome alignments are generated between AM716 and representative cucurbit genomes, combined into a multiple-sequence alignment, and used to calculate nucleotide-level GERP rejected-substitution (RS) scores and conserved elements.

The workflow was adapted from the conservation analysis used for the maize NAM genomes and modified for the cucumber dataset.

## Software

The workflow requires:

- LAST (`lastdb`, `last-train`, `lastal`, `maf-convert`, `last-split`, `maf-swap`, `last-postmask`)
- Kent utilities (`axtChain`, `chainMergeSort`, `faSize`, `chainPreNet`, `chainNet`, `faToTwoBit`, `netToAxt`, `axtToMaf`, `mafSplit`, `mafFilter`, `faSplit`)
- MULTIZ
- PHAST (`msa_view`)
- GERP++
- `seqtk`
- `bedtools`
- Perl
- R

---

## 1. Prepare reference and query genomes

AM716 was used as the cucumber reference genome. The chromosome-level reference sequence and repeat-masked reference sequence were prepared as follows.

```bash
# AM716 reference genome
head -n 14 AM716.assembly.final.allchr0.fa | \
    sed 's/chr/Chr0/g' \
    > AM716.genome.fa

# Repeat-masked AM716 reference genome
seqtk seq -l0 AM716.reTE.fa | \
    head -n 14 | \
    sed 's/chr/Chr0/g' \
    > AM716.rm.genome.fa
```

Representative cucumber genomes were also included.

```bash
head -n 14 AM001.assembly.final.allchr0.fa | sed 's/chr/Chr0/g' > AM001.genome.fa
head -n 14 AM717.assembly.final.allchr0.fa | sed 's/chr/Chr0/g' > AM717.genome.fa
head -n 14 AM739.assembly.final.allchr0.fa | sed 's/chr/Chr0/g' > AM739.genome.fa
head -n 14 AM070.assembly.final.allchr0.fa | sed 's/chr/Chr0/g' > AM070.genome.fa
head -n 14 AM746.assembly.final.allchr0.fa | sed 's/chr/Chr0/g' > AM746.genome.fa
```

Additional cucurbit genomes used for comparative alignment included watermelon, melon, squash, bitter gourd, bottle gourd, sponge gourd, wax gourd, chayote, snake gourd, and *Cucumis hystrix*. Chromosome/scaffold identifiers were standardized before alignment so that sequence names were unique and consistent.

A text file containing the query genome prefixes was used throughout the pipeline:

```text
taxa_12.txt
```

Each entry corresponds to `<taxon>.genome.fa`.

---

## 2. Build the LAST reference database

Build a LAST database from the AM716 reference genome.

```bash
mkdir -p 00_lastdb

lastdb \
    -P 90 \
    -uMAM4 \
    -R01 \
    00_lastdb/MAM4 \
    AM716.genome.fa
```

---

## 3. Generate pairwise whole-genome alignments

### 3.1 Train substitution and gap parameters

```bash
mkdir -p 01_mat

while read -r taxon; do
    last-train \
        --revsym \
        --matsym \
        --gapsym \
        -E0.05 \
        -C2 \
        00_lastdb/MAM4 \
        "${taxon}.genome.fa" \
        > "01_mat/ref_${taxon}.mat" &
done < taxa_12.txt

wait
```

### 3.2 Generate LAST alignments

```bash
mkdir -p 02_maf

while read -r taxon; do
    lastal \
        -P90 \
        -m50 \
        -E0.05 \
        -C2 \
        -p "01_mat/ref_${taxon}.mat" \
        00_lastdb/MAM4 \
        "${taxon}.genome.fa" \
        > "02_maf/ref_${taxon}.maf" &
done < taxa_12.txt

wait
```

---

## 4. Convert pairwise alignments to one-to-one MAF files

### 4.1 Convert MAF to AXT

```bash
mkdir -p 03_axt

while read -r taxon; do
    maf-convert axt \
        "02_maf/ref_${taxon}.maf" \
        > "03_axt/ref_${taxon}.axt" &
done < taxa_12.txt

wait
```

### 4.2 Build chains

```bash
mkdir -p 04_chain

while read -r taxon; do
    axtChain \
        "03_axt/ref_${taxon}.axt" \
        AM716.genome.fa \
        "${taxon}.genome.fa" \
        "04_chain/ref_${taxon}.chain" \
        -linearGap=loose \
        -faQ \
        -faT &
done < taxa_12.txt

wait
```

Merge chains:

```bash
while read -r taxon; do
    chainMergeSort \
        "04_chain/ref_${taxon}.chain" \
        > "04_chain/ref_${taxon}.merged_chain" &
done < taxa_12.txt

wait
```

### 4.3 Generate chromosome-size files and chain nets

```bash
mkdir -p 05_net

faSize AM716.genome.fa -detailed > 05_net/ref_size

while read -r taxon; do
    faSize "${taxon}.genome.fa" -detailed \
        > "05_net/query_${taxon}_size"
done < taxa_12.txt
```

```bash
while read -r taxon; do
    chainPreNet \
        "04_chain/ref_${taxon}.merged_chain" \
        05_net/ref_size \
        "05_net/query_${taxon}_size" \
        "05_net/ref_${taxon}.chain_prenet"

    chainNet \
        "05_net/ref_${taxon}.chain_prenet" \
        05_net/ref_size \
        "05_net/query_${taxon}_size" \
        "05_net/ref_${taxon}.target_net" \
        "05_net/ref_${taxon}.query_net"
done < taxa_12.txt
```

### 4.4 Convert genomes to 2bit and generate net AXT files

```bash
faToTwoBit AM716.genome.fa 05_net/ref.2bit

while read -r taxon; do
    faToTwoBit \
        "${taxon}.genome.fa" \
        "05_net/${taxon}.query.2bit"

    netToAxt \
        "05_net/ref_${taxon}.target_net" \
        "05_net/ref_${taxon}.chain_prenet" \
        05_net/ref.2bit \
        "05_net/${taxon}.query.2bit" \
        "05_net/ref_${taxon}.net.axt"

    axtToMaf \
        "05_net/ref_${taxon}.net.axt" \
        05_net/ref_size \
        "05_net/query_${taxon}_size" \
        "05_net/ref_${taxon}.net.maf"
done < taxa_12.txt
```

### 4.5 Retain one-to-one alignments

```bash
mkdir -p 06_onemaf

while read -r taxon; do
    head -n 29 "02_maf/ref_${taxon}.maf" \
        > "06_onemaf/${taxon}.net_maf_w_header"

    cat "05_net/ref_${taxon}.net.maf" \
        >> "06_onemaf/${taxon}.net_maf_w_header"

    last-split -m1 "06_onemaf/${taxon}.net_maf_w_header" | \
        maf-swap | \
        awk -v q="$taxon" -v r="ref" \
            '/^s/ {$2 = (++s % 2 ? q "." : r ".") $2} 1' | \
        last-split -m1 | \
        maf-swap | \
        last-postmask \
        > "06_onemaf/ref_${taxon}.one_maf" &
done < taxa_12.txt

wait
```

---

## 5. Construct the multiple-sequence alignment

Combine the one-to-one pairwise MAF files using MULTIZ.

```bash
mkdir -p 07_combinemaf
cd 07_combinemaf

maf_array=(../06_onemaf/*.one_maf)
combined_maf=combined.maf

sed -n '/##maf version=1 scoring=blastz/,$p' \
    "${maf_array[0]}" > "${maf_array[0]}_tmp"

sed -n '/##maf version=1 scoring=blastz/,$p' \
    "${maf_array[1]}" > "${maf_array[1]}_tmp"

multiz \
    "${maf_array[0]}_tmp" \
    "${maf_array[1]}_tmp" \
    1 \
    > "$combined_maf"

for maf in "${maf_array[@]:2}"; do
    sed -n '/##maf version=1 scoring=blastz/,$p' \
        "$maf" > "${maf}_tmp"

    multiz \
        "$combined_maf" \
        "${maf}_tmp" \
        1 \
        > "${combined_maf}_tmp"

    mv "${combined_maf}_tmp" "$combined_maf"
done
```

For the cucumber GERP analysis, alignment blocks were required to contain the reference sequence, span at least 20 bp, and contain at least eight aligned sequences.

```bash
mafFilter \
    -minCol=20 \
    -minRow=8 \
    -needComp="ref" \
    combined.maf \
    > combined.maf.filtered
```

---

## 6. Convert the MAF alignment to chromosome-level FASTA alignments

Split the combined MAF by AM716 chromosome.

```bash
mkdir -p split_maf
mafSplit \
    -byTarget dummy.bed \
    split_maf/ \
    combined.maf.filtered \
    -useFullSequenceName
```

Split the unmasked and repeat-masked AM716 reference genomes by chromosome.

```bash
mkdir -p ref_split ref_split_rm

faSplit byname \
    ../AM716.genome.fa \
    ref_split/

faSplit byname \
    ../AM716.rm.genome.fa \
    ref_split_rm/
```

Generate chromosome-level multiple-sequence FASTA alignments and transfer the repeat masking from AM716 to the alignment.

```bash
mkdir -p msa_fasta

for maf_file in split_maf/*Chr*.maf; do
    chr=$(basename "$maf_file" .maf)

    msa_view \
        "$maf_file" \
        -f \
        -G 1 \
        --refseq "ref_split/${chr}.fa" \
        --missing-as-indels \
        > "msa_fasta/${chr}.fa"

    sed -E 's/> />/g' \
        "msa_fasta/${chr}.fa" \
        > "msa_fasta/${chr}.tmp.fa"

    mv "msa_fasta/${chr}.tmp.fa" \
       "msa_fasta/${chr}.fa"

    perl matchMasking.pl \
        --ref "ref_split_rm/${chr}.fa" \
        --fasta "msa_fasta/${chr}.fa" \
        --out "msa_fasta/${chr}_rm.fa"
done
```

The seven AM716 chromosomes are analyzed independently:

```text
Chr01
Chr02
Chr03
Chr04
Chr05
Chr06
Chr07
```

---

## 7. Neutral phylogenetic tree

The phylogenetic tree used by GERP++ should contain exactly the taxa retained in the multiple-sequence alignment. In the cucumber analysis, the reference sequence is represented by `ref`.

The tree used in the analysis was:

```text
(((((ref:0.00006536,AM746:0.00000000):0.01877181,
hystrix:0.02127075):0.01347120,
melon:0.04597876):0.09391855,
(bottlegourd:0.03704591,watermelon:0.02851631):0.04448881):0.03321731,
(chayote:0.14282523,squash:0.12217638):0.02214371,
bittergourd:0.23075412);
```

Save the final tree as:

```text
neutral_tree.txt
```

Bootstrap/support labels should be removed before the tree is supplied to GERP++.

---

## 8. Calculate nucleotide-level GERP scores

Convert each chromosome alignment to single-line FASTA and run `gerpcol`.

```bash
mkdir -p 08_gerp

for i in {1..7}; do
    seqtk seq -l0 \
        "07_combinemaf/msa_fasta/Chr0${i}.fa" | \
        sed 's/> />/g' \
        > "08_gerp/Chr0${i}.gerp.fa"
done
```

```bash
for i in {1..7}; do
    gerpcol \
        -f "08_gerp/Chr0${i}.gerp.fa" \
        -t neutral_tree.txt \
        -v \
        -e ref \
        -j \
        -a &
done

wait
```

`-e ref` excludes the AM716 reference lineage when calculating the neutral expectation.

GERP++ produces a `.rates` file for each chromosome. The second column contains the **RS (Rejected Substitutions) score**, which quantifies evolutionary constraint. Positive RS scores indicate fewer substitutions than expected under neutrality and therefore stronger evolutionary conservation.

Convert nucleotide-level scores to BED format:

```bash
for i in {1..7}; do
    awk -v chr="Chr0${i}" \
        'BEGIN {OFS="\t"} {print chr, NR-1, NR, $1, $2}' \
        "08_gerp/Chr0${i}.gerp.fa.rates" \
        > "08_gerp/Chr0${i}.gerp.rates.bed"
done
```

Extract sites with positive RS scores:

```bash
for i in {1..7}; do
    awk '$5 > 0' \
        "08_gerp/Chr0${i}.gerp.rates.bed" \
        > "08_gerp/Chr0${i}.gerp.positive.bed"
done
```

Combine chromosome-level GERP scores:

```bash
cat 08_gerp/Chr0*.gerp.rates.bed \
    > cucumber_AM716_GERP_scores.bed
```

---

## 9. Identify conserved elements

Run `gerpelem` using the chromosome-level `.rates` files.

```bash
for i in {1..7}; do
    gerpelem \
        -f "08_gerp/Chr0${i}.gerp.fa.rates" &
done

wait
```

Convert the conserved-element output to BED format.

```bash
for i in {1..7}; do
    awk -v chr="Chr0${i}" \
        'BEGIN {OFS="\t"} {
            start=$2-1;
            if (start < 0) start=0;
            print chr, start, $3, $4, $5, $6, $7, $8
        }' \
        "08_gerp/Chr0${i}.gerp.fa.rates.elems" | \
        sort -k1,1 -k2,2n \
        > "08_gerp/Chr0${i}.gerp.elements.bed"
done

cat 08_gerp/Chr0*.gerp.elements.bed \
    > cucumber_AM716_GERP_elements.bed
```

---

## 10. Intersect conserved elements with cucumber gene annotations

Using the AM716 gene annotation:

```bash
grep -v "chr0" AM716.final.gff3 | \
    awk '$3=="CDS"' \
    > AM716.CDS.gff3

grep -v "chr0" AM716.final.gff3 | \
    awk '$3=="gene"' \
    > AM716.gene.gff3
```

Standardize chromosome names in the GERP element file if required.

```bash
sed 's/Chr/chr/g' \
    cucumber_AM716_GERP_elements.bed \
    > cucumber_AM716_GERP_elements.rename.bed
```

Identify conserved elements overlapping coding regions:

```bash
bedtools intersect \
    -a cucumber_AM716_GERP_elements.rename.bed \
    -b AM716.CDS.gff3 \
    > cucumber_AM716_GERP_elements_CDS.bed
```

Identify conserved elements overlapping genes:

```bash
bedtools intersect \
    -a cucumber_AM716_GERP_elements.rename.bed \
    -b AM716.gene.gff3 \
    > cucumber_AM716_GERP_elements_genic.bed
```

---

## Output files

The principal outputs of the cucumber GERP analysis are:

```text
cucumber_AM716_GERP_scores.bed
cucumber_AM716_GERP_elements.bed
cucumber_AM716_GERP_elements_CDS.bed
cucumber_AM716_GERP_elements_genic.bed
```

- `cucumber_AM716_GERP_scores.bed`: nucleotide-level GERP scores across the seven cucumber chromosomes.
- `cucumber_AM716_GERP_elements.bed`: GERP-defined conserved elements.
- `cucumber_AM716_GERP_elements_CDS.bed`: conserved elements overlapping annotated CDS regions.
- `cucumber_AM716_GERP_elements_genic.bed`: conserved elements overlapping annotated genes.

---

## Notes

1. All genome and chromosome identifiers should be standardized before the pairwise alignments are generated.
2. The taxon names in the neutral tree must exactly match the sequence names in the multiple-sequence alignment.
3. The AM716 reference sequence is named `ref` in the multiple alignment and is excluded from the GERP neutral-rate calculation using `-e ref`.
4. The workflow operates on the seven chromosome-level sequences of AM716 and excludes unplaced sequence (`chr0`) from the final chromosome-level summaries.
5. Local paths shown in the original analysis have been replaced by relative paths where possible to make the workflow suitable for a public GitHub repository.
