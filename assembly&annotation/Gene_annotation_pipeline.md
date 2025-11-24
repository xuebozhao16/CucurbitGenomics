# Gene Annotation Pipeline

## Overview
This pipeline annotates genes across 39 *Cucumis sativus* genomes using:

- Repeat masking (RepeatModeler + RepeatMasker)  
- MAKER annotation (AUGUSTUS, SNAP, Spaln, RNA-seq evidence)  
- Gene model transfer (Liftoff)  
- Deep learning annotation (Helixer)  
- Final integration (EVidenceModeler, EVM)  
- BUSCO completeness assessment  

---

## 1. Repeat Identification & Masking (TE_annotation_pipeline.md)

## 2. MAKER Annotation Pipeline

### Transcript Evidence (RNA‑seq from root, stem, leaf, flower, fruit)
```
Trinity --seqType fq --left R1.fq --right R2.fq --CPU 40 --max_memory 200G
```

### Protein Homology Evidence  
Sources include:
- Arabidopsis TAIR10  
- Cucumis melo  
- C. hystrix  
- C. hardwickii  
- Cucumber 9930, Gy14  
- UniProt Swiss-Prot plants  

Aligned using Spaln:
```
spaln -Q7 -O0 -t4 -M 15000 genome.fa proteins.fa > spaln_output.gff3
```

### MAKER Setup  
```
maker -CTL
```

Edit `maker_opts.ctl`:
```
genome=genome.masked.fa
est=trinity.fasta
protein=plant_proteins.fa
snaphmm=snap.hmm
augustus_species=cucumber
keep_preds=1
```

Run MAKER:
```
mpirun -n 40 maker
gff3_merge -d genome_master_datastore_index.log
```

---

## 3. Liftoff Gene Model Transfer
```
liftoff  -g maker_reference.gff3 -o liftoff.gff3 -u unmapped.txt genome.fa reference.fa    -p 40
```

---

## 4. Deep Learning Annotation with Helixer
```
helixer --download-models helixer_models/

helixer    --fasta genome.fa    --model helixer_models/Helixer_Dicot    --out-prefix helixer    --threads 40
```

---

## 5. Ab Initio Predictors: SNAP & AUGUSTUS

### SNAP
```
maker2zff -n maker.gff3
fathom genome.ann genome.dna -categorize 1000
fathom genome.ann genome.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl cucumber_snap hmm > snap.hmm
```

### AUGUSTUS
```
new_species.pl --species=cucumber
etraining --species=cucumber training.gb
augustus --species=cucumber genome.fa > augustus.gff3
```

---

## 6. EVidenceModeler Integration (EVM)

### Weight File Example
```
ABINITIO_PREDICTION    augustus    1
ABINITIO_PREDICTION    snap        1
TRANSCRIPT_ALIGNMENT   pasa        10
PROTEIN_ALIGNMENT      spaln       5
OTHER_PREDICTION       helixer     2
OTHER_PREDICTION       liftoff     8
MAKER_PREDICTION       maker       6
```

### Run EVM  
```
partition_EVM_inputs.pl --genome genome.fa    --gene_predictions *.gff3    --protein_alignments spaln_protein.gff3    --transcript_alignments pasa_assemblies.gff3    --segmentSize 100000 --overlapSize 10000    --partition_listing partitions.list
```

Generate commands:
```
write_EVM_commands.pl --weights weights.txt    --output_file_name evm.out    --partitions partitions.list > commands.list
```

Run:
```
execute_EVM_commands.pl commands.list
```

Recombine:
```
recombine_EVM_partial_outputs.pl --partitions partitions.list  --output_file_name evm.gff3
```

Convert to GFF3:
```
convert_EVM_outputs_to_GFF3.pl  --genome genome.fa   --evm_output_file evm.gff3  --output evm.final.gff3
```

## 7. BUSCO Completeness

```
busco -i final.proteins.fa -l embryophyta_odb10 -m proteins -o busco_out
```
