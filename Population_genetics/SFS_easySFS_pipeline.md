# Estimating the Site-Frequency Spectrum (SFS) with easySFS

## Overview

The **site-frequency spectrum (SFS)** summarizes the distribution of allele frequencies in a population or across multiple populations. It is widely used in demographic inference because different demographic histories leave characteristic signatures in the frequency distribution of genetic variants.

An SFS is commonly used as input for demographic inference programs such as:

- `fastsimcoal2`
- `dadi`

One convenient way to generate an SFS from a filtered VCF file is to use [`easySFS`](https://github.com/isaacovercast/easySFS). `easySFS` is a Python utility built on top of the `dadi` libraries and can generate SFS files in formats suitable for both `fastsimcoal2` and `dadi`.

This workflow describes how to:

1. prepare the working directory;
2. define the input VCF;
3. create the population assignment file;
4. preview projection values;
5. choose projection sizes;
6. generate the SFS;
7. understand the output files;
8. perform basic quality checks before demographic inference.

---

# 1. Software requirements

The workflow requires:

```text
bcftools
Python
easySFS
dadi
```

Before running the analysis, confirm that each program is available:

```bash
bcftools --version
python --version
```

For `easySFS`, the exact command used to check installation depends on how it was installed. A common test is:

```bash
easySFS.py -h
```

or:

```bash
python easySFS.py -h
```

---

# 2. Input data requirements

The primary input is a filtered VCF file containing the individuals to be included in the demographic analysis.

For example:

```text
cichlid_subset_filtered.vcf.gz
```

The VCF should ideally already have undergone quality filtering.

Typical filtering considerations include:

```text
site quality
genotype quality
depth
missing data
minor allele count
biallelic SNP filtering
sample-level missingness
```

For demographic inference, it is especially important that the filtering strategy be consistent across populations.

---

# 3. Create a working directory

It is good practice to keep SFS calculations in a separate analysis directory.

```bash
mkdir -p SFS
cd SFS
```

Copy or link the filtered VCF into the working directory if desired.

```bash
cp cichlid_subset_filtered.vcf.gz ./
```

Define the VCF as a shell variable:

```bash
VCF="cichlid_subset_filtered.vcf.gz"
```

Using a variable makes later commands easier to read and reduces the chance of typographical errors.

Confirm the file exists:

```bash
ls -lh "$VCF"
```

If the VCF is compressed with `bgzip`, make sure it is indexed:

```bash
bcftools index -f "$VCF"
```

This usually generates:

```text
cichlid_subset_filtered.vcf.gz.csi
```

or:

```text
cichlid_subset_filtered.vcf.gz.tbi
```

depending on the indexing method.

---

# 4. Inspect samples in the VCF

Before constructing the population file, list all samples contained in the VCF.

```bash
bcftools query -l "$VCF"
```

The output will contain one sample ID per line.

Example:

```text
Mak.01
Mak.02
Mak.03
Mak.04
Other.01
Other.02
```

Check the number of samples:

```bash
bcftools query -l "$VCF" | wc -l
```

---

# 5. Create the population assignment file

`easySFS` requires a population file assigning each sample to a population.

The format is:

```text
sample_ID    population_ID
```

For example:

```text
Mak.01    Mak
Mak.02    Mak
Mak.03    Mak
Mak.04    Mak
```

In the example dataset, only individuals from Makobe are retained.

The following command extracts samples containing `Mak` and assigns them to a population based on the second field of the sample name:

```bash
bcftools query -l "$VCF" | \
    grep "Mak" | \
    awk '{
        split($0,a,".");
        print $1,a[2]
    }' \
    > pop_file
```

However, the exact parsing command depends on the sample naming convention.

Always inspect the resulting population file:

```bash
cat pop_file
```

Count the number of samples:

```bash
wc -l pop_file
```

---

# 6. General population-file format

For a two-population analysis, the population file might look like:

```text
sample1    Pop1
sample2    Pop1
sample3    Pop1
sample4    Pop1
sample5    Pop2
sample6    Pop2
sample7    Pop2
sample8    Pop2
```

The sample IDs must exactly match the sample IDs in the VCF.

A useful check is:

```bash
cut -f1 pop_file | sort > pop_samples.txt
bcftools query -l "$VCF" | sort > vcf_samples.txt
```

You can then inspect whether population-file samples are present in the VCF.

---

# 7. Why projection is needed

The SFS cannot be computed from sites where the required number of chromosomes is not observed because of missing genotype data.

`easySFS` solves this by allowing **projection**.

Projection means reducing the number of chromosomes sampled from a population so that more genomic sites can contribute to the SFS.

There is a tradeoff:

```text
larger projection
    → more chromosomes retained
    → potentially fewer sites retained

smaller projection
    → fewer chromosomes retained
    → potentially more sites retained
```

The goal is usually to choose a projection that balances:

```text
sample size
and
number of retained sites
```

---

# 8. Preview possible projection values

Use the `--preview` option to inspect how many sites would be retained under different projections.

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview
```

If `easySFS.py` is not in the system `PATH`, run:

```bash
python easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview
```

The preview output shows the number of segregating sites retained for different projection values.

This should be examined before generating the final SFS.

---

# 9. Interpreting projection size

For diploid organisms:

```text
number of sampled chromosomes = 2 × number of diploid individuals
```

For example:

```text
4 diploid individuals = 8 chromosomes
```

Therefore:

```text
--proj 8
```

corresponds to projecting to 8 chromosomes, not necessarily 8 individuals.

For two populations:

```bash
--proj 8,8
```

means:

```text
Population 1: 8 chromosomes
Population 2: 8 chromosomes
```

The number and order of projection values must match the populations in the population file.

---

# 10. Choose the final projection

If missing data are very low, little or no downsampling may be necessary.

For example, if two populations can each support a projection of 8 chromosomes:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8,8
```

This generates the SFS using the selected projection values.

---

# 11. Recommended output directory

It is useful to explicitly specify an output directory.

For example:

```bash
mkdir -p easySFS_output
```

Then run:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8,8 \
    -o easySFS_output
```

This keeps all SFS-related output together.

---

# 12. Folded versus unfolded SFS

An SFS can be either:

```text
folded
or
unfolded
```

## Folded SFS

A folded SFS uses the **minor allele frequency** rather than distinguishing ancestral from derived alleles.

It does not require ancestral-state information.

This is often the safer option when no reliable outgroup or ancestral allele assignment is available.

## Unfolded SFS

An unfolded SFS distinguishes:

```text
ancestral allele
derived allele
```

This requires reliable ancestral-state inference, typically using one or more outgroups.

If ancestral states are incorrectly assigned, the SFS can be biased.

Therefore, the decision to use a folded or unfolded SFS should be made before demographic model fitting.

---

# 13. Missing data considerations

`easySFS` works well for VCF datasets with reasonably confident genotype calls and relatively low missingness.

When missing data are high, projection can recover additional sites, but very low-coverage datasets may be better analyzed using genotype-likelihood-based approaches.

For low-coverage or high-uncertainty data, a program such as:

```text
ANGSD
```

may be preferable because it can estimate the SFS directly from genotype likelihoods rather than relying only on called genotypes.

---

# 14. Single-population SFS

For a single population, the population file could be:

```text
sample1    Pop1
sample2    Pop1
sample3    Pop1
sample4    Pop1
```

Preview projection:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview
```

Generate the final SFS:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8 \
    -o easySFS_output
```

---

# 15. Two-population joint SFS

For two populations:

```text
sample1    Pop1
sample2    Pop1
sample3    Pop1
sample4    Pop1
sample5    Pop2
sample6    Pop2
sample7    Pop2
sample8    Pop2
```

Generate the joint SFS:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8,8 \
    -o easySFS_output
```

The resulting two-dimensional SFS can be used for demographic models involving:

```text
population divergence
population-size change
migration
gene flow
bottlenecks
expansion
```

---

# 16. Three-population joint SFS

The same framework can be extended to three populations.

Example population file:

```text
sample1    Pop1
sample2    Pop1
sample3    Pop2
sample4    Pop2
sample5    Pop3
sample6    Pop3
```

Example:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 6,6,6 \
    -o easySFS_output
```

The number of projection values must equal the number of populations.

---

# 17. Inspect easySFS output

After the run completes:

```bash
find easySFS_output -maxdepth 3 -type f
```

`easySFS` can generate output files suitable for demographic programs such as:

```text
fastsimcoal2
dadi
```

The exact filenames depend on the easySFS version and analysis settings.

Before proceeding, inspect the generated files:

```bash
ls -lh easySFS_output
```

---

# 18. fastsimcoal2-compatible SFS

For `fastsimcoal2`, the SFS is typically represented in files corresponding to the number and type of populations included.

Common SFS types include:

```text
1D SFS
2D joint SFS
multidimensional SFS
```

These files are used together with:

```text
.tpl demographic model file
.est parameter file
```

for demographic parameter estimation.

A conceptual fastsimcoal2 workflow is:

```text
Filtered VCF
    ↓
easySFS
    ↓
Observed SFS
    ↓
fastsimcoal2 model
    ↓
Likelihood optimization
    ↓
Demographic parameter estimates
```

---

# 19. dadi-compatible SFS

`easySFS` can also generate SFS output suitable for `dadi`.

The general dadi workflow is:

```text
Filtered VCF
    ↓
easySFS
    ↓
dadi SFS
    ↓
Demographic model
    ↓
Parameter optimization
    ↓
Model comparison
```

---

# 20. Important issue: SNP-only versus all-sites SFS

Before demographic inference, decide whether the SFS represents:

```text
only variable sites
```

or:

```text
both variable and invariant sites
```

This choice matters because demographic parameter estimation can depend on whether monomorphic sites are represented.

If only SNPs are used, model configuration must be consistent with an SNP-only SFS.

If the total number of callable sites is required by the demographic model, it should be estimated from the genomic regions that passed the same filtering criteria used for the VCF.

---

# 21. Callable sequence length

For many demographic analyses, the total analyzed sequence length is important.

The relevant value is not necessarily the total physical genome size.

Instead, it should represent the number of sites that could have been observed after filtering.

Conceptually:

```text
callable sequence length
=
sites passing coverage filters
+
sites passing missing-data filters
+
sites passing quality filters
+
sites in the genomic regions included in the analysis
```

This quantity may later be used to scale demographic estimates.

---

# 22. Mutation rate and generation time

Programs such as `fastsimcoal2` and `dadi` estimate demographic parameters in units that often need to be converted using:

```text
mutation rate
generation time
```

Before interpreting the final demographic model, record the values used for:

```text
μ = mutation rate per site per generation
g = generation time
```

These values do not affect the easySFS calculation itself, but they are important for interpreting demographic parameters such as:

```text
effective population size
divergence time
migration timing
```

---

# 23. Recommended data checks before running easySFS

## Check the number of samples

```bash
bcftools query -l "$VCF" | wc -l
```

## Check sample IDs

```bash
bcftools query -l "$VCF"
```

## Check the number of variants

```bash
bcftools view -H "$VCF" | wc -l
```

## Check missingness

A missingness summary can be generated with:

```bash
vcftools \
    --gzvcf "$VCF" \
    --missing-site \
    --out missingness
```

## Check allele count distribution

```bash
bcftools +fill-tags "$VCF" \
    -Oz \
    -o VCF_with_AF.vcf.gz \
    -- -t AC,AN,AF
```

---

# 24. Recommended preprocessing for demographic inference

Before generating the SFS, it is often useful to restrict the input to high-confidence biallelic SNPs.

Example:

```bash
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    "$VCF" \
    -Oz \
    -o demographic_snps.vcf.gz
```

Index:

```bash
bcftools index -f demographic_snps.vcf.gz
```

Then:

```bash
VCF="demographic_snps.vcf.gz"
```

This makes the demographic input more explicit and reproducible.

---

# 25. Optional: remove selected or problematic genomic regions

If the goal is neutral demographic inference, users may choose to exclude genomic regions likely to violate neutral-model assumptions.

Examples may include:

```text
coding regions
strongly selected regions
large structural variants
centromeres
low-recombination regions
poorly mappable regions
repeats
```

The exact exclusion strategy depends on the study design.

If regions are excluded, the callable sequence length should be calculated from the same retained regions.

---

# 26. Reproducible one-population example

```bash
#!/usr/bin/env bash
set -euo pipefail

mkdir -p SFS
cd SFS

VCF="population.filtered.vcf.gz"

bcftools index -f "$VCF"

bcftools query -l "$VCF" \
    | awk '{print $1, "Pop1"}' \
    OFS="\t" \
    > pop_file

easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview

easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 20 \
    -o easySFS_output
```

---

# 27. Reproducible two-population example

Assume the VCF contains two populations:

```text
PopA
PopB
```

Population file:

```text
sample01    PopA
sample02    PopA
sample03    PopA
sample04    PopA
sample05    PopB
sample06    PopB
sample07    PopB
sample08    PopB
```

Preview:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview
```

Final SFS:

```bash
easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8,8 \
    -o easySFS_output
```

---

# 28. Choosing projection values systematically

Projection values should not simply be chosen based on the total number of available chromosomes.

Instead, compare:

```text
projection size
versus
number of retained sites
```

A practical strategy is:

```text
1. Run --preview.
2. Examine the number of retained SNPs for each projection.
3. Identify where increasing projection size causes a large loss of sites.
4. Choose a projection before that sharp decline.
5. Use biologically comparable projections across populations when possible.
```

For example:

```text
Projection    Retained sites
6             950,000
8             930,000
10            910,000
12            700,000
14            400,000
```

A projection of 10 may provide a reasonable compromise in this hypothetical example.

---

# 29. Important consistency checks

Before demographic inference, verify:

```text
Population order in the SFS
Projection order
Population labels
Sample IDs
Folded/unfolded status
Number of retained sites
Number of analyzed populations
Chromosome/scaffold inclusion
Callable sequence length
Mutation rate
Generation time
```

Population order is particularly important in multidimensional SFS analyses because demographic model parameters are usually linked to population positions.

---

# 30. Suggested output structure

A clean analysis directory could be:

```text
SFS/
├── input/
│   └── demographic_snps.vcf.gz
├── pop/
│   └── pop_file
├── preview/
│   └── projection_preview.txt
├── output/
│   └── easySFS_output/
└── README.md
```

Large VCF files do not need to be stored in a public GitHub repository. Instead, the repository can contain:

```text
workflow
population-file template
projection parameters
small example files
software versions
```

---

# 31. Suggested manuscript reporting

For reproducibility, report at least:

```text
software and version
input VCF filtering criteria
number of individuals per population
population assignment
number of populations
projection values
folded or unfolded SFS
number of retained sites
whether invariant sites were included
callable sequence length
```

An example Methods-style description is:

> The site-frequency spectrum was generated from the filtered SNP dataset using easySFS. Individuals were assigned to populations using a two-column population file containing sample and population identifiers. Projection values were selected after evaluating the number of retained segregating sites across alternative projections with the easySFS preview function. The final projected SFS was generated in formats compatible with downstream demographic inference.

---

# 32. Complete workflow summary

```text
Filtered VCF
     ↓
Inspect samples
     ↓
Create population assignment file
     ↓
Restrict to demographic SNP set
     ↓
easySFS --preview
     ↓
Choose projection values
     ↓
Generate projected SFS
     ↓
Inspect output
     ↓
Record callable sequence length
     ↓
fastsimcoal2 / dadi
     ↓
Demographic model fitting
```

---

# 33. Final example command set

```bash
mkdir -p SFS
cd SFS

VCF="demographic_snps.vcf.gz"

bcftools index -f "$VCF"

bcftools query -l "$VCF" | \
    grep "Mak" | \
    awk '{
        split($0,a,".");
        print $1,a[2]
    }' \
    OFS="\t" \
    > pop_file

easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --preview

easySFS.py \
    -i "$VCF" \
    -p pop_file \
    -a \
    -f \
    --proj 8,8 \
    -o easySFS_output
```

The resulting SFS files can then be used as input for downstream demographic inference with `fastsimcoal2` or `dadi`.
