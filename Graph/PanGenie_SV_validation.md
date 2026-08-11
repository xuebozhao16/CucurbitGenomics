# Validation of PanGenie-based structural variant genotyping in *Cucurbita pepo*

This document contains the complete reproducible workflow used to evaluate the accuracy of PanGenie-based structural variant (SV) genotyping. The validation uses a leave-one-assembly-out strategy: for each accession, a pangenome graph is constructed without that accession, short reads are generated from the held-out assembly, and PanGenie genotypes the held-out accession against the reduced graph. Assembly-derived SV genotypes are treated as the reference truth.

## Validation design

Nine *C. pepo* assemblies were used:

`C31`, `C38`, `C39`, `JOL`, `OVF`, `SRP`, `TRF`, `VSP`, and `114A`.

For each validation replicate:

1. One assembly is excluded from graph construction.
2. The remaining eight assemblies are used to construct a Minigraph-Cactus pangenome graph.
3. Paired-end short reads are generated from the held-out assembly.
4. PanGenie is used to genotype the held-out accession.
5. Insertions and deletions ≥20 bp are retained.
6. PanGenie calls are compared with assembly-derived SV calls.
7. Accuracy is additionally evaluated for TE-associated SVs and allele-frequency classes.

## Software

The workflow requires Minigraph-Cactus (`cactus-pangenome`), PanGenie, `bcftools`, `samtools`, `seqtk`, `vcfbub`, `vcf-concat`, `bedtools`, Python 3, and the helper scripts used to prepare the PanGenie input VCF.

## Definition of structural variants

After decomposition of multiallelic records, sequence-resolved variants were classified by REF/ALT length. Insertions satisfy `length(ALT) > length(REF)` and deletions satisfy `length(REF) > length(ALT)`. Only variants with an absolute length difference of at least 20 bp were retained.

## Definition of TE regions

TE regions are taken directly from the `repeat_region` features in the C39 EDTA annotation. GFF3 coordinates are converted to BED coordinates as `[start-1, end)` and overlapping TE intervals are merged. **No flanking sequence is added to the TE intervals.**

```bash
awk 'BEGIN{OFS="\t"}
     $0 !~ /^#/ && $3=="repeat_region" {
         print $1, $4-1, $5
     }' C39.assembly.final.allchr0.fa.mod.EDTA.TEanno.gff3 \
> TE.repeat_region.bed

sort -k1,1V -k2,2n TE.repeat_region.bed \
    | bedtools merge -i - \
    > TE.repeat_region.merged.bed
```

## Accuracy metrics

At the site level, variants are matched using exact `CHROM + POS + REF + ALT` identity, with the assembly-based call set treated as truth:

- **TP:** variant present in both assembly and PanGenie call sets.
- **FP:** variant present only in the PanGenie call set.
- **FN:** variant present only in the assembly call set.
- **TN:** not defined for a variant-only VCF.

For genotype-level evaluation at common exact sites, diploid genotypes are additionally converted to SV presence/absence:

- `0/0` → SV absent
- `0/1`, `1/0`, or `1/1` → SV present
- missing genotype → excluded from valid genotype pairs

Precision, recall, F1 score, accuracy, false-positive rate, false-negative rate, exact genotype concordance, binary genotype concordance, and missing rate are reported.

---

# graph leave one out
## minigraph-cactus 
```bash
source /data/zhaojiantao/tools/cactus/cactus_env/bin/activate 
#!/bin/bash
GENOMEDIR=/data/xuebo/pepo/Assembly/07_final2/final_add_chr0
genomes=("C31" "C38" "C39" "JOL" "OVF" "SRP" "TRF" "VSP" "114A")
for leave in "${genomes[@]}"
do
    mkdir -p MC_no_${leave}/genomes
    for sample in "${genomes[@]}"
    do
        [[ "$sample" == "$leave" ]] && continue
        for i in {1..20}
        do
            samtools faidx ${GENOMEDIR}/${sample}.assembly.final.allchr0.fa chr${i} \
                | seqtk seq -l0 \
                > MC_no_${leave}/genomes/${sample}.chr${i}.fa
        done
    done
done
######
#!/bin/bash
genomes=("C31" "C38" "C39" "JOL" "OVF" "SRP" "TRF" "VSP" "114A")
for leave in "${genomes[@]}"
do
    dir="MC_no_${leave}"
    outfile="${dir}/chr1.evolver.txt"
    > "$outfile"
    for sample in "${genomes[@]}"
    do
        [[ "$sample" == "$leave" ]] && continue

        echo -e "${sample}\tgenomes/${sample}.chr1.fa" >> "$outfile"
    done
done
######
for leave in {C31,C38,C39,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    for i in {2..20}
    do
        sed "s/chr1/chr${i}/g" chr1.evolver.txt > chr${i}.evolver.txt
    done
    cd ..
done
####runing
source /data/zhaojiantao/tools/cactus/cactus_env/bin/activate ## 134 134 
## First process /data/xuebo/pepo/R1/graph/MC_no_C39， using --reference C39 
export PATH=/data/zhaojiantao/tools/cactus/build-pangenome-tools/samtools-1.11:$PATH
for i in {1..10}
do
    nohup cactus-pangenome ./js.chr${i} ./chr${i}.evolver.txt  --mapCores 7 --permissiveContigFilter 0.1 --maxLen 10000 --clip 10000 --outDir chr${i}.out --outName chr${i} --reference TRF  --vcf --giraffe --gfa --gbz --odgi > nohupout_chr${i} 2>& 1 &
done
for i in {11..20}
do
    nohup cactus-pangenome ./js.chr${i} ./chr${i}.evolver.txt  --mapCores 7 --permissiveContigFilter 0.1 --maxLen 10000 --clip 10000 --outDir chr${i}.out --outName chr${i} --reference TRF  --vcf --giraffe --gfa --gbz --odgi > nohupout_chr${i} 2>& 1 &
done
#########run the other eight leave-one-out graphs using C39 as reference
## First process /data/xuebo/pepo/R1/graph/MC_no_C39， using --reference C39 
export PATH=/data/zhaojiantao/tools/cactus/build-pangenome-tools/samtools-1.11:$PATH
hash -r
BASE_DIR=$(pwd)
for leave in C31 C38 JOL OVF SRP TRF VSP 114A
do
    workdir="${BASE_DIR}/MC_no_${leave}"
    for i in {1..20}
    do
        # wait when 10 jobs are already running
        while [ "$(jobs -rp | wc -l)" -ge 10 ]
        do
            sleep 30
        done
        nohup bash -c "
            cd '${workdir}' || exit 1
            cactus-pangenome \
                ./js.chr${i} \
                ./chr${i}.evolver.txt \
                --mapCores 7 \
                --permissiveContigFilter 0.1 \
                --maxLen 10000 \
                --clip 10000 \
                --outDir chr${i}.out \
                --outName chr${i} \
                --reference C39 \
                --vcf \
                --giraffe \
                --gfa \
                --gbz \
                --odgi
        " > "${workdir}/nohupout_chr${i}" 2>&1 &
        echo "Started MC_no_${leave}/chr${i}, PID=$!"
    done
done
wait
echo "All cactus-pangenome jobs finished."  ## Transfer results to the analysis server if required

#####Collect results into one directory
for leave in {C31,C38,C39,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    mkdir result
    cd result
    for i in {1..20}
    do
        cp ../chr${i}.out/chr${i}.vcf.gz* ./
        #cp ../chr${i}.out/chr${i}.stats.tgz ./
        /data/xuebo/software/vcfbub-0.1.0/vcfbub -l 0 -r 100000 --input chr${i}.vcf.gz > vcfbub.chr${i}.vcf 
        bcftools norm -m -any vcfbub.chr${i}.vcf > normalized_bi_chr${i}.vcf
        grep -v "#" normalized_bi_chr${i}.vcf | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == "0") $i = "0|0"; else if ($i == "1") $i = "1|1"; else if ($i == ".") $i = ".|." } print }' > chr${i}.modified_file.vcf
        grep "#" normalized_bi_chr${i}.vcf > header${i}.txt
        cat header${i}.txt chr${i}.modified_file.vcf > chr${i}.modified_file2.vcf
    done
    cd ..
    cd ..
done

###### Prepare PanGenie input VCFs for leave-one-out graphs except C39
conda deactivate
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    rm -rf forpangenie && mkdir forpangenie
    cd forpangenie
    for i in {1..20}
    do
        #cp ../result/chr${i}.modified_file2.vcf ./
        perl /data/xuebo/cucumber/graph/tools/fillmiss_VCF_v1.pl $0 ../result/chr${i}.modified_file2.vcf > chr${i}_fillMiss.vcf
        bgzip -c chr${i}_fillMiss.vcf > chr${i}_fillMiss.vcf.gz
        bcftools index --threads 50 --tbi chr${i}_fillMiss.vcf.gz
        bcftools view --threads 50 chr${i}_fillMiss.vcf.gz | perl /data/xuebo/cucumber/graph/tools/update_svlen.pl | perl /data/xuebo/cucumber/graph/tools/merge_ovl_var_forPanGenie.pl /data/xuebo/pepo/R1/graph/MC_no_C31/genomes/C39.chr${i}.fa | bcftools view --threads 60 -Oz -o v1chr${i}_fillMiss.vcf.gz
        bcftools index --threads 60 --tbi v1chr${i}_fillMiss.vcf.gz
    done
    cd ..
    cd ..
done
conda activate
##contig=<ID=chr1,length=26126445>
##contig=<ID=chr2,length=18012633>
##contig=<ID=chr3,length=21830965>
##contig=<ID=chr4,length=14687242>
##contig=<ID=chr5,length=16890741>
##contig=<ID=chr6,length=21953437>
##contig=<ID=chr7,length=16309332>
##contig=<ID=chr8,length=18684281>
##contig=<ID=chr9,length=17695715>
##contig=<ID=chr10,length=19487942>
##contig=<ID=chr11,length=14076622>
##contig=<ID=chr12,length=20760685>
##contig=<ID=chr13,length=18348356>
##contig=<ID=chr14,length=15169833>
##contig=<ID=chr15,length=12448478>
##contig=<ID=chr16,length=22437600>
##contig=<ID=chr17,length=17714425>
##contig=<ID=chr18,length=12288125>
##contig=<ID=chr19,length=16558619>
##contig=<ID=chr20,length=19596639>
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    #mkdir forpangenie
    cd forpangenie
    vcf-concat v1chr1_fillMiss.vcf.gz  v1chr2_fillMiss.vcf.gz  v1chr3_fillMiss.vcf.gz  v1chr4_fillMiss.vcf.gz  v1chr5_fillMiss.vcf.gz  v1chr6_fillMiss.vcf.gz  v1chr7_fillMiss.vcf.gz  v1chr8_fillMiss.vcf.gz  v1chr9_fillMiss.vcf.gz  v1chr10_fillMiss.vcf.gz  v1chr11_fillMiss.vcf.gz  v1chr12_fillMiss.vcf.gz  v1chr13_fillMiss.vcf.gz  v1chr14_fillMiss.vcf.gz  v1chr15_fillMiss.vcf.gz  v1chr16_fillMiss.vcf.gz  v1chr17_fillMiss.vcf.gz  v1chr18_fillMiss.vcf.gz  v1chr19_fillMiss.vcf.gz  v1chr20_fillMiss.vcf.gz > xxx
    cp /data/xuebo/pepo/R1/graph/contig_C39.txt  ./
    grep -v "##contig=" xxx | sed '11r contig_C39.txt' > forpangeniechrall_allhaplo_fillMiss.vcf
    cd ..
    cd ..
done
##grep -v "#" forpangeniechrall_allhaplo_fillMiss.vcf | awk '{for(i=10;i<=NF;i++) count[$i]++} END {for(num in count) print num, count[num]}'

###### Prepare the PanGenie input VCF for the C39 leave-one-out graph
conda deactivate
    cd MC_no_C39
    rm -rf forpangenie && mkdir forpangenie
    cd forpangenie
    for i in {1..20}
    do
        #cp ../result/chr${i}.modified_file2.vcf ./
        perl /data/xuebo/cucumber/graph/tools/fillmiss_VCF_v1.pl $0 ../result/chr${i}.modified_file2.vcf > chr${i}_fillMiss.vcf
        bgzip -c chr${i}_fillMiss.vcf > chr${i}_fillMiss.vcf.gz
        bcftools index --threads 50 --tbi chr${i}_fillMiss.vcf.gz
        bcftools view --threads 50 chr${i}_fillMiss.vcf.gz | perl /data/xuebo/cucumber/graph/tools/update_svlen.pl | perl /data/xuebo/cucumber/graph/tools/merge_ovl_var_forPanGenie.pl /data/xuebo/pepo/R1/graph/MC_no_C39/genomes/TRF.chr${i}.fa | bcftools view --threads 60 -Oz -o v1chr${i}_fillMiss.vcf.gz
        bcftools index --threads 60 --tbi v1chr${i}_fillMiss.vcf.gz
    done
#####
conda activate
##contig=<ID=chr1,length=22615013>
##contig=<ID=chr2,length=18016878>
##contig=<ID=chr3,length=21387330>
##contig=<ID=chr4,length=14745873>
##contig=<ID=chr5,length=14257695>
##contig=<ID=chr6,length=18168248>
##contig=<ID=chr7,length=13809297>
##contig=<ID=chr8,length=16644938>
##contig=<ID=chr9,length=18544679>
##contig=<ID=chr10,length=18593965>
##contig=<ID=chr11,length=18575241>
##contig=<ID=chr12,length=19152503>
##contig=<ID=chr13,length=20545316>
##contig=<ID=chr14,length=11371129>
##contig=<ID=chr15,length=13071278>
##contig=<ID=chr16,length=15521614>
##contig=<ID=chr17,length=18188132>
##contig=<ID=chr18,length=13409077>
##contig=<ID=chr19,length=14631124>
##contig=<ID=chr20,length=15979361>
    cd MC_no_C39
    #mkdir forpangenie
    cd forpangenie
    vcf-concat v1chr1_fillMiss.vcf.gz  v1chr2_fillMiss.vcf.gz  v1chr3_fillMiss.vcf.gz  v1chr4_fillMiss.vcf.gz  v1chr5_fillMiss.vcf.gz  v1chr6_fillMiss.vcf.gz  v1chr7_fillMiss.vcf.gz  v1chr8_fillMiss.vcf.gz  v1chr9_fillMiss.vcf.gz  v1chr10_fillMiss.vcf.gz  v1chr11_fillMiss.vcf.gz  v1chr12_fillMiss.vcf.gz  v1chr13_fillMiss.vcf.gz  v1chr14_fillMiss.vcf.gz  v1chr15_fillMiss.vcf.gz  v1chr16_fillMiss.vcf.gz  v1chr17_fillMiss.vcf.gz  v1chr18_fillMiss.vcf.gz  v1chr19_fillMiss.vcf.gz  v1chr20_fillMiss.vcf.gz > xxx
    cp /data/xuebo/pepo/R1/graph/contig_TRF.txt  ./
    grep -v "##contig=" xxx | sed '11r contig_TRF.txt' > forpangeniechrall_allhaplo_fillMiss.vcf
```

## pangenie
```bash
## Build PanGenie indexes using C39 as the reference genome
conda activate pangenie
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    rm -rf pangenie && mkdir pangenie
    cd pangenie
    /data/xuebo/software/pangenie/build/src/PanGenie-index -v ../forpangenie/forpangeniechrall_allhaplo_fillMiss.vcf -r /data/xuebo/pepo/graph/minigraph_cactus/genomes/C39.chr_all.fa -t 110 -o pangenie14.index
    cd ..
    cd ..
done
## Build the PanGenie index using TRF as the reference genome
cat TRF.chr1.fa  TRF.chr2.fa  TRF.chr3.fa  TRF.chr4.fa  TRF.chr5.fa  TRF.chr6.fa  TRF.chr7.fa  TRF.chr8.fa  TRF.chr9.fa  TRF.chr10.fa  TRF.chr11.fa  TRF.chr12.fa  TRF.chr13.fa  TRF.chr14.fa  TRF.chr15.fa  TRF.chr16.fa  TRF.chr17.fa  TRF.chr18.fa  TRF.chr19.fa  TRF.chr20.fa > TRF.chr_all.fa
rm -rf pangenie && mkdir pangenie
cd pangenie
nohup /data/xuebo/software/pangenie/build/src/PanGenie-index -v ../forpangenie/forpangeniechrall_allhaplo_fillMiss.vcf  -r /data/xuebo/pepo/R1/graph/MC_no_C39/genomes/TRF.chr_all.fa -t 110 -o pangenie14.index &

######## Generate paired-end short reads from the assemblies
for leave in {C31,C38,C39,JOL,OVF,SRP,TRF,VSP,114A}
do
    WGS --model fasta --type toPEfastq --file /data/xuebo/pepo/Assembly/07_final2/final_add_chr0/${leave}.assembly.final.allchr0.fa --out ${leave}.short_reads
done
########pangenie
conda activate pangenie
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    cd pangenie
    /data/xuebo/software/pangenie/build/src/PanGenie -f pangenie14.index -i /data/xuebo/pepo/R1/graph/short_reads/${leave}.short_reads_all.fastq -s no_${leave} -j 30 -t 30 -o no_${leave}
    cd ..
    cd ..
done
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    cd pangenie
    grep "#" no_${leave}_genotyping.vcf > header.txt
    grep -v "#" no_${leave}_genotyping.vcf | awk 'BEGIN {OFS="\t"} {for(i=1; i<=NF; i++) sub(/:.*/, "", $i)}1' > temp
    cat header.txt temp  > PanGenie_simplified.vcf
    sed '11r /data/xuebo/pepo/R1/graph/contig_C39.txt' PanGenie_simplified.vcf > PanGenie_simplified.vcf2
    bcftools norm -m -any  PanGenie_simplified.vcf2 >  PanGenie_simplified_bi.vcf
    bcftools view PanGenie_simplified_bi.vcf | bcftools filter --include 'strlen(REF)<strlen(ALT)' | bcftools view -H > ins
    bcftools view PanGenie_simplified_bi.vcf | bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H > del 
    grep -v "#" ins | awk '{if (length($5)-length($4) > 19) print}' | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == ".") $i = "./." } print }'> temp
    grep "#" PanGenie_simplified_bi.vcf > header2.txt
    cat header2.txt temp > ins_no_${leave}.vcf 
    grep -v "#" del | awk '{if (length($4)-length($5) > 19) print}' | awk 'BEGIN {OFS="\t"} { for (i=10; i<=NF; i++) { if ($i == ".") $i = "./." } print }' > temp
    cat header2.txt temp > del_no_${leave}.vcf 
    cd ..
    cd ..
done
##### Sort VCF records
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    cd MC_no_${leave}
    cd pangenie
    (
    grep "^##" ins_no_${leave}.vcf 
    grep "^#CHROM" ins_no_${leave}.vcf | sed 's/no_//g' 
    grep -v "^#" ins_no_${leave}.vcf  | sort -k1,1V -k2,2n 
    ) > /data/xuebo/pepo/R1/graph/F1_score_leave1/ins_no_${leave}.vcf 
    (
    grep "^##" del_no_${leave}.vcf 
    grep "^#CHROM" del_no_${leave}.vcf  | sed 's/no_//g' 
    grep -v "^#" del_no_${leave}.vcf  | sort -k1,1V -k2,2n
    ) > /data/xuebo/pepo/R1/graph/F1_score_leave1/del_no_${leave}.vcf
    cd ..
    cd ..
done
```

## F1 score
```bash
######## Compare assembly-derived genotypes with leave-one-out PanGenie genotypes
cp /data/xuebo/pepo/R1/graph/MC_all/F1_score/Ass9_* ./
gunzip Ass9_ins.vcf.gz 
gunzip Ass9_del.vcf.gz
#################################
#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import OrderedDict
from typing import Dict, List, Optional, Set, Tuple, TextIO, Any


VariantKey = Tuple[str, int, str, str]


def open_text(filename: str) -> TextIO:
    """
    Support both plain-text and gzip-compressed VCF files.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")

    return open(filename, "r")


def chromosome_sort_key(chromosome: str) -> Tuple[int, Any]:
    """
    Use natural ordering for chr1, chr2, ..., chr20.
    Place other chromosome names after numbered chromosomes.
    """
    name = chromosome

    if name.lower().startswith("chr"):
        name = name[3:]

    try:
        return 0, int(name)
    except ValueError:
        return 1, chromosome


def variant_sort_key(key: VariantKey) -> Tuple:
    """
    Sort by chromosome, position, REF, and ALT.
    """
    chrom, pos, ref, alt = key

    return (
        chromosome_sort_key(chrom),
        pos,
        ref,
        alt
    )


def extract_gt(format_field: str, sample_field: str) -> str:
    """
    Extract GT from the FORMAT field.

    Example:
        FORMAT = GT:GQ:GL
        sample = 1/1:99:...

    Returns:
        1/1

    Return "." if GT is unavailable.
    """
    if format_field in {"", "."}:
        return "."

    format_keys = format_field.split(":")
    sample_values = sample_field.split(":")

    if "GT" not in format_keys:
        return "."

    gt_index = format_keys.index("GT")

    if gt_index >= len(sample_values):
        return "."

    gt = sample_values[gt_index]

    if gt == "":
        return "."

    return gt


def normalize_gt(gt: str) -> Optional[str]:
    """
    Normalize diploid genotypes.

    Rules:
        0|0 -> 0/0
        0|1 -> 0/1
        1|0 -> 0/1
        1/0 -> 0/1
        1|1 -> 1/1

    Return None for any genotype containing ".".
    """
    if gt is None:
        return None

    gt = gt.split(":", 1)[0]
    gt = gt.replace("|", "/")

    if "." in gt:
        return None

    alleles = gt.split("/")

    if len(alleles) != 2:
        return None

    # Only alleles 0 and 1 are expected in this dataset
    if any(allele not in {"0", "1"} for allele in alleles):
        return None

    alleles = sorted(alleles)

    return "/".join(alleles)


def gt_to_binary(gt: Optional[str]) -> Optional[int]:
    """
    Convert genotype to SV presence/absence.

        0/0 -> 0
        0/1 -> 1
        1/1 -> 1
        missing -> None
    """
    if gt is None:
        return None

    if gt == "0/0":
        return 0

    return 1


def format_metric(value: Optional[float]) -> str:
    """
     
    """
    if value is None:
        return "NA"

    return f"{value:.6f}"


def read_vcf(
    filename: str,
    target_sample: str
) -> Tuple[List[str], OrderedDict, List[str]]:
    """
    Read a VCF file.

     
        CHROM + POS + REF + ALT
    as the unique variant key.

    Returns:
        meta_headers
        records
        sample_names
    """
    meta_headers: List[str] = []
    records: OrderedDict = OrderedDict()
    sample_names: List[str] = []

    sample_column_index: Optional[int] = None
    found_chrom_header = False

    with open_text(filename) as infile:

        for line_number, line in enumerate(infile, start=1):

            if line.startswith("##"):
                meta_headers.append(line.rstrip("\n"))
                continue

            if line.startswith("#CHROM"):
                found_chrom_header = True

                header_fields = line.rstrip("\n").split("\t")

                if len(header_fields) < 10:
                    raise ValueError(
                        f"{filename} contains no sample columns."
                    )

                sample_names = header_fields[9:]

                if target_sample not in sample_names:
                    raise ValueError(
                        f"Sample '{target_sample}' is not present in {filename}.\n"
                        f"Available samples: {', '.join(sample_names)}"
                    )

                sample_column_index = (
                    9 + sample_names.index(target_sample)
                )

                continue

            if line.startswith("#"):
                continue

            if not found_chrom_header:
                raise ValueError(
                    f"{filename} does not contain a #CHROM header."
                )

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                print(
                    f"Warning: skip malformed line "
                    f"{line_number} in {filename}; "
                    f"only {len(fields)} columns.",
                    file=sys.stderr
                )
                continue

            if sample_column_index is None:
                raise ValueError(
                    f"Unable to determine the sample column in {filename}."
                )

            if sample_column_index >= len(fields):
                print(
                    f"Warning: sample column missing at line "
                    f"{line_number} in {filename}.",
                    file=sys.stderr
                )
                continue

            chrom = fields[0]

            try:
                pos = int(fields[1])
            except ValueError:
                print(
                    f"Warning: invalid POS at line "
                    f"{line_number} in {filename}.",
                    file=sys.stderr
                )
                continue

            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            format_field = fields[8]
            sample_field = fields[sample_column_index]

            gt = extract_gt(format_field, sample_field)

            key: VariantKey = (
                chrom,
                pos,
                ref,
                alt
            )

            record = {
                "chrom": chrom,
                "pos": pos,
                "id": variant_id,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filt,
                "info": info,
                "format": format_field,
                "sample_field": sample_field,
                "gt": gt
            }

            if key in records:
                print(
                    f"Warning: duplicate exact variant in "
                    f"{filename}: "
                    f"{chrom}:{pos}:{ref}:{alt}; "
                    f"keeping the first record.",
                    file=sys.stderr
                )
                continue

            records[key] = record

    return meta_headers, records, sample_names


def calculate_site_metrics(
    assembly_keys: Set[VariantKey],
    shortread_keys: Set[VariantKey]
) -> Dict[str, Any]:
    """
    Site-level comparison.

    Exact matching criterion:
        CHROM + POS + REF + ALTmust be identical

    Assembly calls are treated as truth:
        TP = present in both VCFs
        FP = short-read/PanGenie only
        FN = assembly only
        TN = cannot be defined from variant-only VCFs
    """
    common = assembly_keys & shortread_keys
    assembly_only = assembly_keys - shortread_keys
    shortread_only = shortread_keys - assembly_keys

    tp = len(common)
    fp = len(shortread_only)
    fn = len(assembly_only)

    precision = (
        tp / (tp + fp)
        if (tp + fp) > 0
        else None
    )

    recall = (
        tp / (tp + fn)
        if (tp + fn) > 0
        else None
    )

    if precision is None or recall is None:
        f1 = None
    elif precision + recall > 0:
        f1 = (
            2 * precision * recall
            / (precision + recall)
        )
    else:
        f1 = 0.0

    return {
        "AssemblyTotal": len(assembly_keys),
        "ShortreadTotal": len(shortread_keys),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": None,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "CommonKeys": common,
        "AssemblyOnlyKeys": assembly_only,
        "ShortreadOnlyKeys": shortread_only
    }


def calculate_genotype_metrics(
    common_keys: List[VariantKey],
    assembly_records: OrderedDict,
    shortread_records: OrderedDict
) -> Dict[str, Any]:
    """
     CHROM POS REF ALTmust be identical 
    compare assembly and short-read genotypes.

    Binary classification:
        0/0 -> SV absent
        0/1、1/0、1/1 -> SV present

    Assembly truth 
    """
    tp = 0
    fp = 0
    fn = 0
    tn = 0

    both_missing = 0
    assembly_missing_only = 0
    shortread_missing_only = 0

    valid_pairs = 0

    exact_gt_correct = 0
    binary_gt_correct = 0

    exact_00 = 0
    exact_01 = 0
    exact_11 = 0

    for key in common_keys:

        assembly_gt = normalize_gt(
            assembly_records[key]["gt"]
        )

        shortread_gt = normalize_gt(
            shortread_records[key]["gt"]
        )

        # both genotypes missing
        if (
            assembly_gt is None
            and shortread_gt is None
        ):
            both_missing += 1
            continue

        # assembly genotype missing only
        if assembly_gt is None:
            assembly_missing_only += 1
            continue

        # short-read genotype missing only
        if shortread_gt is None:
            shortread_missing_only += 1
            continue

        valid_pairs += 1

        # exact genotype match
        if assembly_gt == shortread_gt:
            exact_gt_correct += 1

            if assembly_gt == "0/0":
                exact_00 += 1
            elif assembly_gt == "0/1":
                exact_01 += 1
            elif assembly_gt == "1/1":
                exact_11 += 1

        truth = gt_to_binary(assembly_gt)
        prediction = gt_to_binary(shortread_gt)

        # binary genotype match
        if truth == prediction:
            binary_gt_correct += 1

        if truth == 1 and prediction == 1:
            tp += 1

        elif truth == 0 and prediction == 1:
            fp += 1

        elif truth == 1 and prediction == 0:
            fn += 1

        elif truth == 0 and prediction == 0:
            tn += 1

    precision = (
        tp / (tp + fp)
        if (tp + fp) > 0
        else None
    )

    recall = (
        tp / (tp + fn)
        if (tp + fn) > 0
        else None
    )

    if precision is None or recall is None:
        f1 = None
    elif precision + recall > 0:
        f1 = (
            2 * precision * recall
            / (precision + recall)
        )
    else:
        f1 = 0.0

    accuracy = (
        (tp + tn) / valid_pairs
        if valid_pairs > 0
        else None
    )

    exact_gt_concordance = (
        exact_gt_correct / valid_pairs
        if valid_pairs > 0
        else None
    )

    binary_gt_concordance = (
        binary_gt_correct / valid_pairs
        if valid_pairs > 0
        else None
    )

    total_common = len(common_keys)

    missing_pairs = (
        both_missing
        + assembly_missing_only
        + shortread_missing_only
    )

    missing_rate = (
        missing_pairs / total_common
        if total_common > 0
        else None
    )

    false_negative_rate = (
        fn / (tp + fn)
        if (tp + fn) > 0
        else None
    )

    false_positive_rate = (
        fp / (fp + tn)
        if (fp + tn) > 0
        else None
    )

    return {
        "CommonSites": total_common,
        "ValidPairs": valid_pairs,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "Accuracy": accuracy,
        "FNR": false_negative_rate,
        "FPR": false_positive_rate,
        "ExactGTCorrect": exact_gt_correct,
        "ExactGTConcordance": exact_gt_concordance,
        "BinaryGTCorrect": binary_gt_correct,
        "BinaryGTConcordance": binary_gt_concordance,
        "Exact00": exact_00,
        "Exact01": exact_01,
        "Exact11": exact_11,
        "BothMissing": both_missing,
        "AssemblyMissingOnly": assembly_missing_only,
        "ShortreadMissingOnly": shortread_missing_only,
        "MissingPairs": missing_pairs,
        "MissingRate": missing_rate
    }


def write_common_vcf(
    output_file: str,
    assembly_headers: List[str],
    assembly_records: OrderedDict,
    shortread_records: OrderedDict,
    common_keys: List[VariantKey]
) -> None:
    """
    Write a two-sample VCF for exact common sites.

    Samples:
        C31_assembly
        C31_shortread
    """
    with open(output_file, "w") as out:

        has_gt_header = False

        for header in assembly_headers:
            if header.startswith("##FORMAT=<ID=GT,"):
                has_gt_header = True

            out.write(header + "\n")

        if not has_gt_header:
            out.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,'
                'Description="Genotype">\n'
            )

        out.write(
            "##comparisonMethod=Exact match using "
            "CHROM,POS,REF,ALT\n"
        )

        out.write(
            "##assemblySample=C31_assembly\n"
        )

        out.write(
            "##predictionSample=C31_shortread\n"
        )

        out.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            "\tFORMAT\tC31_assembly\tC31_shortread\n"
        )

        for key in common_keys:

            assembly_record = assembly_records[key]
            shortread_record = shortread_records[key]

            # Use GT only to avoid differences in FORMAT definitions between VCFs
            assembly_gt = assembly_record["gt"]
            shortread_gt = shortread_record["gt"]

            out.write(
                "\t".join([
                    assembly_record["chrom"],
                    str(assembly_record["pos"]),
                    assembly_record["id"],
                    assembly_record["ref"],
                    assembly_record["alt"],
                    assembly_record["qual"],
                    assembly_record["filter"],
                    assembly_record["info"],
                    "GT",
                    assembly_gt,
                    shortread_gt
                ]) + "\n"
            )


def write_variant_table(
    output_file: str,
    keys: List[VariantKey],
    records: OrderedDict
) -> None:
    """
    Write a table of dataset-specific variants.
    """
    with open(output_file, "w") as out:

        out.write(
            "CHROM\tPOS\tID\tREF\tALT\tGT\n"
        )

        for key in keys:

            record = records[key]

            out.write(
                "\t".join([
                    record["chrom"],
                    str(record["pos"]),
                    record["id"],
                    record["ref"],
                    record["alt"],
                    record["gt"]
                ]) + "\n"
            )


def write_genotype_detail_table(
    output_file: str,
    common_keys: List[VariantKey],
    assembly_records: OrderedDict,
    shortread_records: OrderedDict
) -> None:
    """
    Write genotype comparison details for every common exact site.
    """
    with open(output_file, "w") as out:

        out.write(
            "CHROM\tPOS\tREF\tALT"
            "\tAssemblyGT\tShortreadGT"
            "\tExactGTMatch\tBinaryTruth"
            "\tBinaryPrediction\tClassification\n"
        )

        for key in common_keys:

            assembly_raw_gt = assembly_records[key]["gt"]
            shortread_raw_gt = shortread_records[key]["gt"]

            assembly_gt = normalize_gt(assembly_raw_gt)
            shortread_gt = normalize_gt(shortread_raw_gt)

            exact_match = "NA"
            truth_text = "NA"
            prediction_text = "NA"
            classification = "NA"

            if (
                assembly_gt is None
                and shortread_gt is None
            ):
                exact_match = "BothMissing"
                classification = "BothMissing"

            elif assembly_gt is None:
                exact_match = "No"
                classification = "AssemblyMissing"

            elif shortread_gt is None:
                exact_match = "No"
                classification = "ShortreadMissing"

            else:
                exact_match = (
                    "Yes"
                    if assembly_gt == shortread_gt
                    else "No"
                )

                truth = gt_to_binary(assembly_gt)
                prediction = gt_to_binary(shortread_gt)

                truth_text = str(truth)
                prediction_text = str(prediction)

                if truth == 1 and prediction == 1:
                    classification = "TP"
                elif truth == 0 and prediction == 1:
                    classification = "FP"
                elif truth == 1 and prediction == 0:
                    classification = "FN"
                elif truth == 0 and prediction == 0:
                    classification = "TN"

            chrom, pos, ref, alt = key

            out.write(
                "\t".join([
                    chrom,
                    str(pos),
                    ref,
                    alt,
                    assembly_raw_gt,
                    shortread_raw_gt,
                    exact_match,
                    truth_text,
                    prediction_text,
                    classification
                ]) + "\n"
            )


def write_stats(
    output_file: str,
    site_metrics: Dict[str, Any],
    genotype_metrics: Dict[str, Any]
) -> None:
    """
    Write complete validation statistics.
    """
    with open(output_file, "w") as out:

        out.write("[Exact_site_matching]\n")

        out.write(
            f"Assembly_total\t"
            f"{site_metrics['AssemblyTotal']}\n"
        )

        out.write(
            f"Shortread_total\t"
            f"{site_metrics['ShortreadTotal']}\n"
        )

        out.write(
            f"Site_TP_common_exact\t"
            f"{site_metrics['TP']}\n"
        )

        out.write(
            f"Site_FP_shortread_only\t"
            f"{site_metrics['FP']}\n"
        )

        out.write(
            f"Site_FN_assembly_only\t"
            f"{site_metrics['FN']}\n"
        )

        # Genome-wide TN cannot be defined from a variant-only VCF
        out.write("Site_TN\tNA\n")

        out.write(
            f"Site_Precision\t"
            f"{format_metric(site_metrics['Precision'])}\n"
        )

        out.write(
            f"Site_Recall\t"
            f"{format_metric(site_metrics['Recall'])}\n"
        )

        out.write(
            f"Site_F1\t"
            f"{format_metric(site_metrics['F1'])}\n"
        )

        out.write("\n")

        out.write(
            "[Genotype_comparison_at_exact_common_sites]\n"
        )

        out.write(
            f"Common_exact_sites\t"
            f"{genotype_metrics['CommonSites']}\n"
        )

        out.write(
            f"GT_ValidPairs\t"
            f"{genotype_metrics['ValidPairs']}\n"
        )

        out.write(
            f"GT_TP\t"
            f"{genotype_metrics['TP']}\n"
        )

        out.write(
            f"GT_FP\t"
            f"{genotype_metrics['FP']}\n"
        )

        out.write(
            f"GT_FN\t"
            f"{genotype_metrics['FN']}\n"
        )

        out.write(
            f"GT_TN\t"
            f"{genotype_metrics['TN']}\n"
        )

        out.write(
            f"GT_Precision\t"
            f"{format_metric(genotype_metrics['Precision'])}\n"
        )

        out.write(
            f"GT_Recall\t"
            f"{format_metric(genotype_metrics['Recall'])}\n"
        )

        out.write(
            f"GT_F1\t"
            f"{format_metric(genotype_metrics['F1'])}\n"
        )

        out.write(
            f"GT_Accuracy\t"
            f"{format_metric(genotype_metrics['Accuracy'])}\n"
        )

        out.write(
            f"GT_FalseNegativeRate\t"
            f"{format_metric(genotype_metrics['FNR'])}\n"
        )

        out.write(
            f"GT_FalsePositiveRate\t"
            f"{format_metric(genotype_metrics['FPR'])}\n"
        )

        out.write(
            f"Exact_GT_correct\t"
            f"{genotype_metrics['ExactGTCorrect']}\n"
        )

        out.write(
            f"Exact_GT_concordance\t"
            f"{format_metric(genotype_metrics['ExactGTConcordance'])}\n"
        )

        out.write(
            f"Binary_GT_correct\t"
            f"{genotype_metrics['BinaryGTCorrect']}\n"
        )

        out.write(
            f"Binary_GT_concordance\t"
            f"{format_metric(genotype_metrics['BinaryGTConcordance'])}\n"
        )

        out.write(
            f"Exact_0_0_matches\t"
            f"{genotype_metrics['Exact00']}\n"
        )

        out.write(
            f"Exact_0_1_matches\t"
            f"{genotype_metrics['Exact01']}\n"
        )

        out.write(
            f"Exact_1_1_matches\t"
            f"{genotype_metrics['Exact11']}\n"
        )

        out.write(
            f"Both_missing\t"
            f"{genotype_metrics['BothMissing']}\n"
        )

        out.write(
            f"Assembly_missing_only\t"
            f"{genotype_metrics['AssemblyMissingOnly']}\n"
        )

        out.write(
            f"Shortread_missing_only\t"
            f"{genotype_metrics['ShortreadMissingOnly']}\n"
        )

        out.write(
            f"Missing_pairs\t"
            f"{genotype_metrics['MissingPairs']}\n"
        )

        out.write(
            f"Missing_rate\t"
            f"{format_metric(genotype_metrics['MissingRate'])}\n"
        )


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Compare an assembly-based C31 VCF with a "
            "leave-one-out PanGenie short-read C31 VCF."
        )
    )

    parser.add_argument(
        "--assembly-vcf",
        required=True,
        help="Assembly-based VCF"
    )

    parser.add_argument(
        "--shortread-vcf",
        required=True,
        help="Leave-one-out PanGenie short-read VCF"
    )

    parser.add_argument(
        "--assembly-sample",
        default="C31",
        help=(
            "Sample name in assembly VCF "
            "(default: C31)"
        )
    )

    parser.add_argument(
        "--shortread-sample",
        default="C31",
        help=(
            "Sample name in short-read VCF "
            "(default: C31)"
        )
    )

    parser.add_argument(
        "--prefix",
        default="C31_compare",
        help=(
            "Output prefix "
            "(default: C31_compare)"
        )
    )

    args = parser.parse_args()

    print("Reading assembly VCF...", file=sys.stderr)

    (
        assembly_headers,
        assembly_records,
        _
    ) = read_vcf(
        args.assembly_vcf,
        args.assembly_sample
    )

    print("Reading short-read VCF...", file=sys.stderr)

    (
        _,
        shortread_records,
        _
    ) = read_vcf(
        args.shortread_vcf,
        args.shortread_sample
    )

    assembly_keys = set(assembly_records.keys())
    shortread_keys = set(shortread_records.keys())

    site_metrics = calculate_site_metrics(
        assembly_keys,
        shortread_keys
    )

    # Sort common sites using natural chromosome order
    common_keys = sorted(
        site_metrics["CommonKeys"],
        key=variant_sort_key
    )

    assembly_only_keys = sorted(
        site_metrics["AssemblyOnlyKeys"],
        key=variant_sort_key
    )

    shortread_only_keys = sorted(
        site_metrics["ShortreadOnlyKeys"],
        key=variant_sort_key
    )

    genotype_metrics = calculate_genotype_metrics(
        common_keys,
        assembly_records,
        shortread_records
    )

    stats_file = f"{args.prefix}.stats.txt"
    common_vcf = f"{args.prefix}.common.vcf"

    assembly_only_file = (
        f"{args.prefix}.assembly_only.tsv"
    )

    shortread_only_file = (
        f"{args.prefix}.shortread_only.tsv"
    )

    genotype_detail_file = (
        f"{args.prefix}.genotype_details.tsv"
    )

    write_stats(
        stats_file,
        site_metrics,
        genotype_metrics
    )

    write_common_vcf(
        common_vcf,
        assembly_headers,
        assembly_records,
        shortread_records,
        common_keys
    )

    write_variant_table(
        assembly_only_file,
        assembly_only_keys,
        assembly_records
    )

    write_variant_table(
        shortread_only_file,
        shortread_only_keys,
        shortread_records
    )

    write_genotype_detail_table(
        genotype_detail_file,
        common_keys,
        assembly_records,
        shortread_records
    )

    print("\nExact site matching:")
    print(
        f"  Assembly total       : "
        f"{site_metrics['AssemblyTotal']}"
    )
    print(
        f"  Short-read total     : "
        f"{site_metrics['ShortreadTotal']}"
    )
    print(
        f"  TP (common exact)    : "
        f"{site_metrics['TP']}"
    )
    print(
        f"  FP (short-read only) : "
        f"{site_metrics['FP']}"
    )
    print(
        f"  FN (assembly only)   : "
        f"{site_metrics['FN']}"
    )
    print("  TN                    : NA")
    print(
        f"  Precision             : "
        f"{format_metric(site_metrics['Precision'])}"
    )
    print(
        f"  Recall                : "
        f"{format_metric(site_metrics['Recall'])}"
    )
    print(
        f"  F1                    : "
        f"{format_metric(site_metrics['F1'])}"
    )

    print("\nGenotype comparison at common exact sites:")
    print(
        f"  Common sites          : "
        f"{genotype_metrics['CommonSites']}"
    )
    print(
        f"  Valid GT pairs        : "
        f"{genotype_metrics['ValidPairs']}"
    )
    print(
        f"  TP                    : "
        f"{genotype_metrics['TP']}"
    )
    print(
        f"  FP                    : "
        f"{genotype_metrics['FP']}"
    )
    print(
        f"  FN                    : "
        f"{genotype_metrics['FN']}"
    )
    print(
        f"  TN                    : "
        f"{genotype_metrics['TN']}"
    )
    print(
        f"  Precision             : "
        f"{format_metric(genotype_metrics['Precision'])}"
    )
    print(
        f"  Recall                : "
        f"{format_metric(genotype_metrics['Recall'])}"
    )
    print(
        f"  F1                    : "
        f"{format_metric(genotype_metrics['F1'])}"
    )
    print(
        f"  Accuracy              : "
        f"{format_metric(genotype_metrics['Accuracy'])}"
    )
    print(
        f"  Exact GT concordance  : "
        f"{format_metric(genotype_metrics['ExactGTConcordance'])}"
    )
    print(
        f"  Binary concordance    : "
        f"{format_metric(genotype_metrics['BinaryGTConcordance'])}"
    )

    print("\nOutput files:")
    print(f"  {stats_file}")
    print(f"  {common_vcf}")
    print(f"  {assembly_only_file}")
    print(f"  {shortread_only_file}")
    print(f"  {genotype_detail_file}")


if __name__ == "__main__":
    main()

#####
python compare_two_C31_vcf.py \
    --assembly-vcf Ass9_ins.vcf \
    --shortread-vcf ins_no_C31.vcf \
    --assembly-sample C31 \
    --shortread-sample C31 \
    --prefix C31_exact_compare

##########################
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf Ass9_ins.vcf --shortread-vcf ins_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_INS
done
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf Ass9_del.vcf --shortread-vcf del_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_DEL
done
mkdir stat
cp *_INS.stats.txt ./stat
cp *_DEL.stats.txt ./stat
```

## F1 score (TE region)
```bash
#### Extract TE-associated SVs and evaluate validation metrics
awk 'BEGIN{OFS="\t"}
     $0 !~ /^#/ && $3=="repeat_region" {
         print $1, $4-1, $5
     }' C39.assembly.final.allchr0.fa.mod.EDTA.TEanno.gff3 \
> TE.repeat_region.bed
sort -k1,1V -k2,2n TE.repeat_region.bed | bedtools merge -i - > TE.repeat_region.merged.bed
cp ../F1_score_leave1/Ass9_ins.vcf ./
cp ../F1_score_leave1/Ass9_del.vcf ./
## Swap REF and ALT where needed for interval overlap
awk 'BEGIN{FS=OFS="\t"}
/^#/{
    print
    next
}
{
    tmp=$4
    $4=$5
    $5=tmp
    print
}' Ass9_ins.vcf > Ass9_ins.swap.vcf
awk 'BEGIN{FS=OFS="\t"}
/^#/{
    print
    next
}
{
    tmp=$4
    $4=$5
    $5=tmp
    print
}' Ass9_del.vcf > Ass9_del.swap.vcf
####
bgzip -c Ass9_ins.vcf > Ass9_ins.vcf.gz
tabix Ass9_ins.vcf.gz
bgzip -c Ass9_ins.swap.vcf > Ass9_ins.swap.vcf.gz
tabix Ass9_ins.swap.vcf.gz
bgzip -c Ass9_del.vcf > Ass9_del.vcf.gz
tabix Ass9_del.vcf.gz
bgzip -c Ass9_del.swap.vcf > Ass9_del.swap.vcf.gz
tabix Ass9_del.swap.vcf.gz
bcftools view -R TE.repeat_region.merged.bed Ass9_ins.vcf.gz -Ov -o Ass9_ins.TE.vcf
bcftools view -R TE.repeat_region.merged.bed Ass9_ins.swap.vcf.gz -Ov -o Ass9_ins.swap.TE.vcf
bcftools view -R TE.repeat_region.merged.bed Ass9_del.vcf.gz -Ov -o Ass9_del.TE.vcf
bcftools view -R TE.repeat_region.merged.bed Ass9_del.swap.vcf.gz -Ov -o Ass9_del.swap.TE.vcf
####
awk 'BEGIN{FS=OFS="\t"}
/^#/{
    print
    next
}
{
    tmp=$4
    $4=$5
    $5=tmp
    print
}' Ass9_del.swap.TE.vcf > Ass9_del.swap2.TE.vcf
vcf-concat Ass9_del.TE.vcf Ass9_del.swap2.TE.vcf > xxx
(
    bcftools view -h xxx
    bcftools view -H xxx | sort -u
) > merged_del.TE.unique.vcf
####
awk 'BEGIN{FS=OFS="\t"}
/^#/{
    print
    next
}
{
    tmp=$4
    $4=$5
    $5=tmp
    print
}' Ass9_ins.swap.TE.vcf > Ass9_ins.swap2.TE.vcf
vcf-concat Ass9_ins.TE.vcf Ass9_ins.swap2.TE.vcf > xxxx
(
    bcftools view -h xxxx
    bcftools view -H xxxx | sort -u
) > merged_ins.TE.unique.vcf

####################
cp ../F1_score_leave1/compare_two_C31_vcf.py ./
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf merged_ins.TE.unique.vcf --shortread-vcf ../F1_score_leave1/ins_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_INS
done
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf merged_del.TE.unique.vcf --shortread-vcf ../F1_score_leave1/del_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_DEL
done
mkdir stat_TE
cp *_INS.stats.txt ./stat_TE
cp *_DEL.stats.txt ./stat_TE 
```

## F1 score rare site
```bash
######## Select rare SVs using allele frequencies from the population dataset
Input files: /data/xuebo/pepo/pop_genetics/00_207_data/final_SV2_change_name/Insersion_sorted.vcf && /data/xuebo/pepo/pop_genetics/00_207_data/final_SV2_change_name/Deletion_sorted.vcf
## Define rare insertion SVs using the population VCF
bcftools +fill-tags Insersion_sorted.vcf -Ov  -o Insersion_sorted.fillAF.vcf -- -t AC,AN,AF
## Extract matching variants from the assembly-derived VCF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' Insersion_sorted.rare.vcf > rare.keys.tsv
awk 'BEGIN {
         FS=OFS="\t"
     }
     NR==FNR {
         key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4
         rare[key]=1
         next
     }
     /^#/ {
         print
         next
     }
     {
         key=$1 SUBSEP $2 SUBSEP $4 SUBSEP $5

         if (key in rare) {
             print
         }
     }
' \
rare.keys.tsv ../F1_score_leave1/Ass9_ins.vcf > Ass9_ins.rare.vcf
## Define rare deletion SVs using the population VCF
bcftools +fill-tags Deletion_sorted.vcf -Ov  -o Deletion_sorted.fillAF.vcf -- -t AC,AN,AF
## Extract matching variants from the assembly-derived VCF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' Deletion_sorted.rare.vcf > rare.keys.del.tsv
awk 'BEGIN {
         FS=OFS="\t"
     }
     NR==FNR {
         key=$1 SUBSEP $2 SUBSEP $3 SUBSEP $4
         rare[key]=1
         next
     }
     /^#/ {
         print
         next
     }
     {
         key=$1 SUBSEP $2 SUBSEP $4 SUBSEP $5

         if (key in rare) {
             print
         }
     }
' \
rare.keys.del.tsv ../F1_score_leave1/Ass9_del.vcf > Ass9_del.rare.vcf

#########################
cp ../F1_score_leave1/compare_two_C31_vcf.py ./
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf Ass9_ins.rare.vcf --shortread-vcf ../F1_score_leave1/ins_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_INS
done
for leave in {C31,C38,JOL,OVF,SRP,TRF,VSP,114A}
do
    python compare_two_C31_vcf.py --assembly-vcf Ass9_del.rare.vcf --shortread-vcf ../F1_score_leave1/del_no_${leave}.vcf --assembly-sample ${leave} --shortread-sample ${leave}   --prefix ${leave}_DEL
done
mkdir stat_rare
cp *_INS.stats.txt ./stat_rare
cp *_DEL.stats.txt ./stat_rare 
```












