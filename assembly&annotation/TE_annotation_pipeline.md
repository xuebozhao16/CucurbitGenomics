# TE Annotation Pipeline

Following the NAM pipeline:
https://github.com/HuffordLab/NAM-genomes/tree/master/te-annotation

This workflow describes the transposable element (TE) annotation process for cucumber genomes:

1. TE identification and classification using EDTA
2. TE masking using RepeatMasker
3. Construction of a pan-genome TE library
4. Summary of TE annotation results
5. TE presence/absence variation (PAV) profiling

## 1. Install EDTA and Dependencies
```bash
git clone https://github.com/oushujun/EDTA.git
cd EDTA
conda env create -f EDTA.yml
conda activate EDTA
conda install -c bioconda repeatmasker bedtools seqtk
```

## 2. Run EDTA on Cucumber Assemblies
```bash
#!/bin/bash
max_jobs=4
count=0

while read genome; do
    mkdir -p EDTA_${genome}
    (
        cd EDTA_${genome}
        perl ../EDTA.pl --genome ../genomes/${genome}.fa --cds ../reference/ChineseLong_CDS.fa --anno 1 --sensitive 1 --threads 40
    ) &
    ((count++))

    if [ "$count" -eq "$max_jobs" ]; then
        wait
        count=0
    fi
done < accession_list.txt

wait
```

## 3. Extract TE Regions
```bash
for g in $(cat accession_list.txt); do
    awk 'NR>6 {print $1"\t"$4-1"\t"$5}' EDTA_${g}/${g}.mod.EDTA.TEanno.gff3 > ${g}.TE.bed
done
```

Mask TE regions:
```bash
bedtools maskfasta -mc N -fi genomes/${g}.fa -bed ${g}.TE.bed -fo ${g}.masked.fa
```

## 4. Build Pan-Genome TE Library
```bash
for g in $(cat accession_list.txt); do
  perl output_by_list.pl 1 <(grep -v "#" ${g}.TElib.fa) 1 <(awk '{print $10}' ${g}.mod.out | sort | uniq) -FA > ${g}.TElib.novel.fa
done
```

Merge & rename:
```bash
perl rename_TE.pl *.novel.fa > panTE.raw.fa
```

Remove nested and redundant TEs:
```bash
perl cleanup_nested.pl -in panTE.raw.fa -cov 0.95 -minlen 80 -miniden 80 -t 90 > panEDTA.TElib.fa
```

## 5. RepeatMasker Annotation Using Pan-TE Library
```bash
RepeatMasker -pa 36 -q -div 40 -lib panEDTA.TElib.fa -gff genome.fa
```

## 6. TE Summary Statistics
```bash
for g in $(cat accession_list.txt); do
  perl buildSummary.pl -maxDiv 60 -stats ${g}.stats ${g}.out > ${g}.TEanno.sum
done

perl transpose3.pl *.TEanno.sum > TE_summary_matrix.txt
```

## 7. TE PAV (Presence/Absence Variation)
```bash
paste TE_name.txt number2_*.txt > PAV_matrix.txt

awk '{
  c=0;
  for(i=2;i<=NF;i++) if($i==1) c++;
  print $1"\t"c;
}' PAV_matrix.txt > TE_frequency.txt
```
