# Demographic History Inference: fastsimcoal2, PSMC, MSMC/MSMC2, and SMC++

## Overview

This document provides an integrated workflow for demographic-history inference using **fastsimcoal2**, **PSMC**, **MSMC/MSMC2**, and **SMC++**.

These methods answer complementary demographic questions. `fastsimcoal2` fits explicit demographic models to the site-frequency spectrum (SFS), whereas PSMC, MSMC/MSMC2, and SMC++ use patterns of heterozygosity and coalescence along the genome to reconstruct changes in effective population size through time and, for multi-population analyses, population separation or cross-coalescence patterns.

The workflow covers:

1. preparing observed SFS files;
2. defining demographic models in `.tpl` files;
3. defining parameter ranges in `.est` files;
4. running repeated fastsimcoal2 optimizations;
5. selecting the best run for each model;
6. comparing alternative demographic models;
7. evaluating model fit;
8. estimating likelihood distributions;
9. performing block bootstrap analyses;
10. obtaining confidence intervals for demographic parameters.

The examples below use a two-population system, but the same framework can be extended to more populations.

---

# 1. Software requirements

The core workflow requires:

```text
fastsimcoal2
easySFS
bcftools
R
```

Additional helper utilities may be used for:

```text
best-run selection
AIC calculation
SFS visualization
demographic model plotting
bootstrap summarization
```

Before starting, record the software version used in the final analysis.

Example:

```bash
fastsimcoal2 -h
bcftools --version
R --version
```

---

# 2. Required input files

For each demographic model, fastsimcoal2 requires three files sharing the same prefix.

For example, with:

```bash
PREFIX="early_geneflow"
```

the files are:

```text
early_geneflow.tpl
early_geneflow.est
early_geneflow_jointMAFpop1_0.obs
```

or, for an unfolded SFS:

```text
early_geneflow_jointDAFpop1_0.obs
```

The three required input types are:

```text
Observed SFS
Template file (.tpl)
Estimation file (.est)
```

---

# 3. Observed SFS

## 3.1 Folded and unfolded SFS

The observed SFS can be either:

```text
folded
or
unfolded
```

A folded SFS is represented as a minor-allele-frequency spectrum.

Typical naming:

```text
PREFIX_MAFpop0.obs
PREFIX_jointMAFpop1_0.obs
```

An unfolded SFS is represented as a derived-allele-frequency spectrum.

Typical naming:

```text
PREFIX_DAFpop0.obs
PREFIX_jointDAFpop1_0.obs
```

For more than two populations, a multidimensional SFS can be used.

Typical suffixes include:

```text
_DSFS.obs
_MSFS.obs
```

and fastsimcoal2 is run with:

```bash
-multiSFS
```

when appropriate.

---

## 3.2 Generate the SFS

The observed SFS can be generated using tools such as:

```text
easySFS
ANGSD
Arlequin
```

For the workflow here, `easySFS` is assumed.

Example:

```bash
easySFS.py \
    -i demographic_snps.vcf.gz \
    -p pop_file \
    -a \
    -f \
    --proj 8,8 \
    -o SFS_output
```

After creating the SFS, copy or rename the relevant observed SFS file so that its prefix matches the demographic model.

Example:

```bash
cp SFS_output/dadi/pop1-pop2.sfs \
   early_geneflow_jointMAFpop1_0.obs
```

The exact source filename depends on the easySFS version and output format.

---

# 4. The template file

## 4.1 Purpose of the `.tpl` file

The template file defines the demographic model.

It includes information such as:

```text
number of populations
current effective population sizes
sample sizes
growth rates
migration matrices
historical events
population splits
population mergers
changes in migration
changes in population size
```

Parameters that should be estimated are represented by symbolic names.

Examples:

```text
NPOP1
NPOP2
TDIV
M12
M21
TGF
```

Parameter names should be unique and should not be substrings of other parameter names or reserved expressions.

---

## 4.2 Demographic models are specified backward in time

fastsimcoal2 describes demographic history **backward in time**.

Therefore:

```text
most recent event
    ↓
older event
    ↓
oldest event
```

must be written in that order in the historical-event section of the `.tpl` file.

This is one of the most important points when constructing a demographic model.

---

# 5. Example two-population demographic models

Several alternative two-population models can be compared.

Examples include:

```text
No gene flow
Constant gene flow
Early gene flow
Recent gene flow
Different early and recent gene flow
```

These represent different biological hypotheses.

---

## 5.1 Model 1: no gene flow

Conceptual history:

```text
Present
│
│ Pop1              Pop2
│  │                  │
│  │                  │
│  └────────┬─────────┘
│           │
│        divergence
│
Past
```

After divergence, no migration occurs between populations.

Typical model parameters:

```text
N1
N2
NANC
TDIV
```

---

## 5.2 Model 2: constant gene flow

Conceptual history:

```text
Present
│
│ Pop1  ⇄  Pop2
│   continuous migration
│
│  └────────┬─────────┘
│           │
│        divergence
│
Past
```

Typical parameters:

```text
N1
N2
NANC
TDIV
M12
M21
```

---

## 5.3 Model 3: early gene flow

This model represents divergence with gene flow immediately after the split, followed by reproductive isolation.

Forward-time interpretation:

```text
ancestral population
        ↓
     divergence
        ↓
 early gene flow
        ↓
 complete isolation
        ↓
     present
```

Because fastsimcoal2 is written backward in time, the isolation period is encountered first, followed by the older migration phase.

Typical parameters:

```text
N1
N2
NANC
TDIV
TISO
M12
M21
```

---

## 5.4 Model 4: recent gene flow

This model represents secondary contact.

Forward-time interpretation:

```text
ancestral population
        ↓
     divergence
        ↓
     isolation
        ↓
 secondary contact
        ↓
     present
```

Typical parameters:

```text
N1
N2
NANC
TDIV
TCONTACT
M12
M21
```

---

## 5.5 Model 5: different early and recent gene flow

This model allows two different migration phases.

Conceptually:

```text
Present
│
│ recent migration
│
│ isolation or reduced migration
│
│ ancient migration
│
│ divergence
│
Past
```

Typical parameters:

```text
N1
N2
NANC
TDIV
T1
T2
M12_RECENT
M21_RECENT
M12_OLD
M21_OLD
```

---

# 6. Example template structure

The exact syntax depends on the demographic model, but a simplified two-population template may look like:

```text
//Number of population samples
2

//Population effective sizes
NPOP1
NPOP2

//Sample sizes
8
8

//Growth rates
0
0

//Number of migration matrices
1

//Migration matrix
0 M12
M21 0

//Historical events
1 historical_event
TDIV 0 1 1 1 0 0

//Number of independent loci
1 0

//Per chromosome: number of linkage blocks
1

//Per block: data type, number of loci, recombination rate, mutation rate
FREQ 1 0 0
```

This example is only a structural illustration. The final `.tpl` file must match the exact demographic model and SFS configuration.

---

# 7. The estimation file

## 7.1 Purpose of the `.est` file

The `.est` file defines:

```text
which parameters are estimated
whether parameters are integer or continuous
prior distribution
minimum value
maximum value
whether values are written to output
```

Every symbolic parameter used in the `.tpl` file must be defined in the `.est` file.

---

## 7.2 Simple parameters

A typical simple-parameter section may look like:

```text
[PARAMETERS]

//isInt? name distribution min max output

1 NPOP1 unif 100 1000000 output
1 NPOP2 unif 100 1000000 output
1 NANC  unif 100 1000000 output

1 TDIV  unif 100 100000 output

0 M12   logunif 1e-8 1e-2 output
0 M21   logunif 1e-8 1e-2 output
```

General interpretation:

```text
1 = integer parameter
0 = continuous parameter
```

Common distributions include:

```text
unif
logunif
```

---

## 7.3 Complex parameters

Complex parameters can be derived from simple parameters.

Examples include:

```text
ratios
differences
minimum values
maximum values
time constraints
```

This is useful when one event must be older than another.

A conceptual example:

```text
TOLD = TYOUNG + DELTA
```

This prevents impossible time ordering.

---

# 8. Parameter ranges

Parameter ranges should be biologically plausible but sufficiently broad.

Examples:

```text
Population size:
10^2 – 10^7

Divergence time:
10^2 – 10^6 generations

Migration rate:
10^-8 – 10^-2
```

These ranges are only examples.

The final parameter boundaries should be chosen based on:

```text
biology of the species
previous literature
generation time
mutation rate
population structure
pilot fastsimcoal2 runs
```

If a parameter consistently reaches the minimum or maximum boundary, the range should be reconsidered.

---

# 9. Mutation rate and parameter scaling

Demographic parameters estimated by fastsimcoal2 are scaled according to the mutation model.

If a reliable mutation rate is available, it should be included in the model.

If a reliable mutation rate is not available, one demographic parameter may need to be fixed to provide a scale for the model.

For example:

```text
fix divergence time
estimate population sizes and migration
```

or:

```text
fix mutation rate
estimate divergence time
```

Without an external scale, some combinations of demographic parameters and mutation rate are not independently identifiable.

---

# 10. Create a model directory

Use one directory per demographic model.

```bash
mkdir -p fastsimcoal
cd fastsimcoal

PREFIX="early_geneflow"

mkdir -p "$PREFIX"
cd "$PREFIX"
```

Place the three required files in the directory:

```text
early_geneflow.tpl
early_geneflow.est
early_geneflow_jointMAFpop1_0.obs
```

---

# 11. Run a first fastsimcoal2 optimization

Example:

```bash
fastsimcoal2 \
    -t "${PREFIX}.tpl" \
    -e "${PREFIX}.est" \
    -m \
    -0 \
    -C 10 \
    -n 10000 \
    -L 40 \
    -s 0 \
    -M
```

---

# 12. Explanation of important options

## `-t`

Template file:

```bash
-t "${PREFIX}.tpl"
```

---

## `-e`

Parameter-estimation file:

```bash
-e "${PREFIX}.est"
```

---

## `-m`

Use a folded/minor-allele-frequency SFS.

```bash
-m
```

For an unfolded/derived SFS, use the appropriate unfolded-SFS option instead.

---

## `-0`

Ignore monomorphic sites.

```bash
-0
```

This option is appropriate when the observed SFS does not include invariant sites.

---

## `-C 10`

Pool SFS entries containing fewer than 10 SNPs.

```bash
-C 10
```

This can reduce overfitting when many SFS bins have very small counts.

---

## `-n`

Number of coalescent simulations used to approximate the expected SFS during each optimization cycle.

Tutorial-scale example:

```bash
-n 10000
```

For a final analysis, substantially larger values are preferable.

---

## `-L`

Number of ECM optimization cycles.

Tutorial-scale example:

```bash
-L 40
```

For final analyses, more cycles are often used.

---

## `-M`

Perform parameter estimation.

```bash
-M
```

---

## `-q`

Quiet mode.

```bash
-q
```

Useful when launching many replicate runs.

---

# 13. Recommended final-analysis settings

For initial debugging:

```text
n = 10,000
L = 20–40
few replicate runs
```

For final parameter estimation, a more intensive analysis may use:

```text
n = 200,000–1,000,000
L = 50–100
100 or more independent optimization replicates
```

The final values depend on computational resources and model complexity.

---

# 14. Inspect fastsimcoal2 outputs

After optimization, fastsimcoal2 creates a directory named after the prefix.

Important output files include:

```text
PREFIX.bestlhoods
PREFIX_jointMAFpop1_0.txt
PREFIX.simparam
PREFIX_maxL.par
```

---

## 14.1 `.bestlhoods`

Contains:

```text
best parameter estimates
maximum observed likelihood
maximum estimated likelihood
```

The two likelihood columns of particular interest are:

```text
MaxObsLhood
MaxEstLhood
```

A better-fitting model has a smaller difference between these two likelihoods.

---

## 14.2 Expected SFS

Example:

```text
PREFIX_jointMAFpop1_0.txt
```

This is the expected SFS generated under the best-fitting parameter estimates.

It should be compared visually with the observed SFS.

---

## 14.3 `.simparam`

This file contains the fully resolved simulation configuration.

It is useful for checking whether:

```text
population sizes
event times
migration matrices
historical events
```

were interpreted correctly.

---

## 14.4 `_maxL.par`

This contains the model with estimated parameter values substituted into the template.

It can later be used to simulate datasets under the best-fitting model.

---

# 15. Run many independent optimization replicates

One fastsimcoal2 optimization is not sufficient because the optimizer may converge on a local optimum.

Run many independent replicates.

Example with 100 runs:

```bash
PREFIX="early_geneflow"

for i in $(seq 1 100); do

    mkdir -p "run${i}"

    cp \
        "${PREFIX}.tpl" \
        "${PREFIX}.est" \
        "${PREFIX}_jointMAFpop1_0.obs" \
        "run${i}/"

    (
        cd "run${i}"

        fastsimcoal2 \
            -t "${PREFIX}.tpl" \
            -e "${PREFIX}.est" \
            -m \
            -0 \
            -C 10 \
            -n 200000 \
            -L 50 \
            -s 0 \
            -M \
            -q
    )

done
```

For computational clusters, these independent runs can be submitted as separate jobs or as a job array.

---

# 16. Select the best run

The best run is typically the one with the highest `MaxEstLhood`.

Inspect all runs:

```bash
for file in run*/"${PREFIX}"/"${PREFIX}.bestlhoods"; do
    tail -n 1 "$file"
done
```

A simple ranking command can be adapted to the exact column structure of the `.bestlhoods` file.

For example:

```bash
for i in $(seq 1 100); do

    file="run${i}/${PREFIX}/${PREFIX}.bestlhoods"

    awk -v run="$i" '
        NR==2 {
            print run, $0
        }
    ' "$file"

done > all_runs.bestlhoods.txt
```

Then sort by the `MaxEstLhood` column after confirming its position in the file.

---

# 17. Copy the best run

After identifying the best run:

```bash
BEST_RUN=42

mkdir -p bestrun

cp -r \
    "run${BEST_RUN}/${PREFIX}/"* \
    bestrun/
```

The `bestrun` directory can then be used for:

```text
model comparison
SFS visualization
likelihood re-evaluation
bootstrap analyses
parameter reporting
```

---

# 18. Alternative demographic models

Repeat the full optimization procedure for each candidate model.

For example:

```text
no_geneflow
ongoing_geneflow
early_geneflow
recent_geneflow
diff_geneflow
```

Each model should have:

```text
MODEL.tpl
MODEL.est
MODEL_jointMAFpop1_0.obs
```

The observed SFS is the same biological dataset but renamed so that its prefix matches the model.

---

# 19. Model comparison with AIC

Raw likelihood values should not be compared without accounting for differences in the number of estimated parameters.

A common criterion is:

```text
AIC = 2k - 2ln(L)
```

where:

```text
k = number of estimated parameters
L = model likelihood
```

For each model:

```text
extract best MaxEstLhood
count estimated parameters
calculate AIC
```

---

## 19.1 Example R code for AIC calculation

```r
models <- data.frame(
    Model = c(
        "no_geneflow",
        "ongoing_geneflow",
        "early_geneflow",
        "recent_geneflow",
        "diff_geneflow"
    ),

    logLik = c(
        NA,
        NA,
        NA,
        NA,
        NA
    ),

    k = c(
        NA,
        NA,
        NA,
        NA,
        NA
    )
)

models$AIC <- (
    2 * models$k
    - 2 * models$logLik
)

models$DeltaAIC <- (
    models$AIC
    - min(
        models$AIC,
        na.rm = TRUE
    )
)

models <- models[
    order(models$AIC),
]

write.table(
    models,
    "model_comparison_AIC.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

---

# 20. Akaike weights

Akaike weights can also be calculated.

```r
models$AIC_weight <- (
    exp(
        -0.5 * models$DeltaAIC
    )
    /
    sum(
        exp(
            -0.5 * models$DeltaAIC
        )
    )
)
```

The weights provide a relative measure of support among candidate models.

---

# 21. Visualize observed and expected SFS

A demographic model should not be accepted based only on likelihood or AIC.

The expected SFS should also be visually compared with the observed SFS.

Conceptually:

```text
Observed SFS
     vs
Expected SFS under best-fit model
```

Useful comparisons include:

```text
heatmap of observed SFS
heatmap of expected SFS
residual heatmap
difference matrix
```

Large systematic residuals may indicate that the model does not capture important aspects of the population history.

---

# 22. Plot the demographic model

After obtaining the best-fitting parameter estimates, plot the inferred model.

The final figure can show:

```text
effective population sizes
population divergence
migration direction
migration magnitude
migration intervals
event timing
```

This is useful both for interpretation and for detecting model-specification errors.

---

# 23. Likelihood distribution analysis

Because fastsimcoal2 approximates likelihoods using simulations, repeated evaluation of the same best-fitting model can produce slightly different likelihood values.

Therefore, likelihood distributions can be generated for each best-fitting model.

This provides another way to evaluate whether two models are meaningfully different.

---

# 24. Prepare the maximum-likelihood model

Inside the best-run directory:

```bash
PREFIX="early_geneflow"

cp \
    "${PREFIX}_jointMAFpop1_0.obs" \
    "${PREFIX}_maxL_jointMAFpop1_0.obs"
```

The observed SFS must match the prefix of the maximum-likelihood `.par` file.

---

# 25. Recalculate likelihood repeatedly

Example with 100 replicates:

```bash
> "${PREFIX}.lhoods"

for i in $(seq 1 100); do

    fastsimcoal2 \
        -i "${PREFIX}_maxL.par" \
        -n 1000000 \
        -m \
        -q \
        -0

    sed -n '2,3p' \
        "${PREFIX}_maxL/${PREFIX}_maxL.lhoods" \
        >> "${PREFIX}.lhoods"

    rm -rf "${PREFIX}_maxL"

done
```

This produces a likelihood distribution for the fitted model.

Repeat the same procedure for every competing model.

---

# 26. Compare likelihood distributions in R

Example:

```r
no_gene_flow <- scan(
    "no_geneflow.lhoods"
)

ongoing_gene_flow <- scan(
    "ongoing_geneflow.lhoods"
)

early_gene_flow <- scan(
    "early_geneflow.lhoods"
)

recent_gene_flow <- scan(
    "recent_geneflow.lhoods"
)

different_gene_flow <- scan(
    "diff_geneflow.lhoods"
)

boxplot(
    different_gene_flow,
    recent_gene_flow,
    early_gene_flow,
    ongoing_gene_flow,
    no_gene_flow,
    range = 0,
    ylab = "Likelihood",
    names = c(
        "Early + recent",
        "Recent",
        "Early",
        "Constant",
        "None"
    )
)
```

If likelihood distributions strongly overlap, the models may provide very similar fits despite differences in best-run likelihood.

---

# 27. Collect AIC values for all models

A simple shell workflow can combine per-model results.

```bash
> all_models_AIC.txt

for file in */bestrun/*.AIC; do

    model=$(basename "$file" .AIC)

    echo -e \
        "${model}\t$(tail -n 1 "$file")" \
        >> all_models_AIC.txt

done
```

The exact parsing step depends on the format of the AIC output.

---

# 28. Parameter uncertainty

The maximum-likelihood parameter estimates from one observed SFS do not provide uncertainty intervals by themselves.

A common approach is:

```text
block bootstrap of genomic data
       ↓
new bootstrap SFS
       ↓
fit best demographic model
       ↓
best parameter estimates
       ↓
repeat
       ↓
parameter distributions
       ↓
confidence intervals
```

---

# 29. Why block bootstrap?

SNPs are often linked.

Resampling individual SNPs independently would therefore underestimate uncertainty.

Instead, genomic blocks are resampled.

Conceptually:

```text
Genome
│
├── block 1
├── block 2
├── block 3
├── ...
└── block N
```

Bootstrap datasets are generated by sampling blocks with replacement.

---

# 30. Extract VCF header and variant records

Assume:

```bash
VCF="demographic.vcf.gz"
PREFIX="best_model"
```

Extract the VCF header:

```bash
bcftools view -h "$VCF" \
    > header.vcf
```

Extract variant records:

```bash
bcftools view -H "$VCF" \
    > all_sites.txt
```

---

# 31. Split variants into blocks

A simple line-based example:

```bash
split \
    -l 4338 \
    all_sites.txt \
    bootstrap_block_
```

This creates files such as:

```text
bootstrap_block_aa
bootstrap_block_ab
bootstrap_block_ac
...
```

A biologically preferable approach is to define physical genomic windows rather than arbitrary numbers of SNPs whenever chromosome coordinates are available.

For example:

```text
100-kb blocks
500-kb blocks
1-Mb blocks
```

The block size should be chosen with linkage in mind.

---

# 32. Generate bootstrap VCFs

Example with 50 bootstrap replicates and 100 sampled blocks per replicate:

```bash
PREFIX="best_model"

for i in $(seq 1 50); do

    mkdir -p "bs${i}"

    cat header.vcf \
        > "bs${i}/${PREFIX}.bs.${i}.vcf"

    for r in $(seq 1 100); do

        block=$(find . \
            -maxdepth 1 \
            -name 'bootstrap_block_*' \
            | shuf -n 1)

        cat "$block" \
            >> "bs${i}/${PREFIX}.bs.${i}.vcf"

    done

    bgzip -f \
        "bs${i}/${PREFIX}.bs.${i}.vcf"

    bcftools index -f \
        "bs${i}/${PREFIX}.bs.${i}.vcf.gz"

done
```

---

# 33. Generate an SFS for each bootstrap replicate

For every bootstrap VCF:

```bash
for i in $(seq 1 50); do

    easySFS.py \
        -i "bs${i}/${PREFIX}.bs.${i}.vcf.gz" \
        -p pop_file \
        -a \
        -f \
        --proj 8,8 \
        -o "bs${i}/SFS"

done
```

Then copy or rename the generated SFS to match the bootstrap model prefix.

Conceptually:

```text
best_model.bs.1_jointMAFpop1_0.obs
best_model.bs.2_jointMAFpop1_0.obs
...
```

---

# 34. Prepare bootstrap model files

The same best-fitting demographic model is fitted to every bootstrap SFS.

For each bootstrap replicate, copy:

```text
best model .tpl
best model .est
bootstrap observed SFS
```

and rename them consistently.

Example:

```bash
for i in $(seq 1 50); do

    cp \
        "${PREFIX}.tpl" \
        "bs${i}/${PREFIX}.bs.${i}.tpl"

    cp \
        "${PREFIX}.est" \
        "bs${i}/${PREFIX}.bs.${i}.est"

done
```

If model prefixes are used internally in scripts, verify that renaming does not alter the content required by fastsimcoal2.

---

# 35. Fit each bootstrap SFS repeatedly

Each bootstrap replicate should itself be optimized multiple times because fastsimcoal2 can converge on local optima.

Example:

```bash
for bs in $(seq 1 50); do

    cd "bs${bs}"

    BOOT_PREFIX="${PREFIX}.bs.${bs}"

    for run in $(seq 1 100); do

        mkdir -p "run${run}"

        cp \
            "${BOOT_PREFIX}.tpl" \
            "${BOOT_PREFIX}.est" \
            "${BOOT_PREFIX}_jointMAFpop1_0.obs" \
            "run${run}/"

        (
            cd "run${run}"

            fastsimcoal2 \
                -t "${BOOT_PREFIX}.tpl" \
                -e "${BOOT_PREFIX}.est" \
                -m \
                -0 \
                -C 10 \
                -n 200000 \
                -L 50 \
                -s 0 \
                -M \
                -q
        )

    done

    cd ..

done
```

For large analyses, this step is usually parallelized on a computing cluster.

---

# 36. Select the best run for every bootstrap replicate

For each bootstrap replicate:

```text
compare all optimization runs
select highest MaxEstLhood
extract parameter estimates
```

Store one parameter vector per bootstrap replicate.

For example:

```text
bootstrap_parameters.tsv
```

with columns:

```text
Bootstrap
N1
N2
NANC
TDIV
M12
M21
...
```

---

# 37. Calculate confidence intervals

Once bootstrap parameter estimates have been collected, calculate confidence intervals.

Example R code:

```r
params <- read.table(
    "bootstrap_parameters.tsv",
    header = TRUE,
    sep = "\t"
)

parameter_names <- setdiff(
    colnames(params),
    "Bootstrap"
)

CI <- data.frame()

for (parameter in parameter_names) {

    x <- params[[parameter]]

    result <- data.frame(
        Parameter = parameter,
        Median = median(
            x,
            na.rm = TRUE
        ),
        Lower_2.5 = quantile(
            x,
            0.025,
            na.rm = TRUE
        ),
        Upper_97.5 = quantile(
            x,
            0.975,
            na.rm = TRUE
        )
    )

    CI <- rbind(
        CI,
        result
    )
}

write.table(
    CI,
    "bootstrap_parameter_CI.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

---

# 38. Convert generations to years

If fastsimcoal2 estimates time in generations, convert to years using generation time:

```text
Time in years
=
Time in generations × generation time
```

Example in R:

```r
generation_time <- 1

params$TDIV_years <- (
    params$TDIV
    * generation_time
)
```

Use the generation time appropriate for the study organism.

---

# 39. Interpret migration parameters carefully

Migration parameters can be represented in different ways depending on model specification.

Possible quantities include:

```text
per-generation migration probability
number of migrants per generation
scaled migration
direction-specific migration
```

The parameter names in the `.tpl` and `.est` files should clearly indicate direction.

For example:

```text
M12 = migration from population 1 to population 2
M21 = migration from population 2 to population 1
```

Be consistent in both model files and figures.

---

# 40. Check parameter boundaries

After selecting the best run, inspect whether estimated parameters occur near prior boundaries.

For each parameter:

```text
estimate
minimum bound
maximum bound
```

If:

```text
estimate ≈ minimum
```

or:

```text
estimate ≈ maximum
```

the optimization range may be too narrow or the parameter may be poorly identifiable.

Such parameters should be investigated before final interpretation.

---

# 41. Check replicate convergence

Do not rely only on the single best run.

Inspect the likelihood distribution among replicate optimizations.

Useful diagnostics include:

```text
top 10 likelihoods
parameter values among top runs
difference between best and second-best run
frequency of convergence to similar parameter values
```

If many independent runs converge to similar likelihoods and parameter estimates, confidence in the optimum is greater.

---

# 42. Recommended model-fitting workflow

For each demographic model:

```text
1. Build .tpl file.
2. Build .est file.
3. Use the same observed SFS.
4. Perform test runs.
5. Check .simparam.
6. Increase simulation number.
7. Run ≥100 independent optimizations.
8. Select best run.
9. Inspect observed vs expected SFS.
10. Record MaxEstLhood.
11. Calculate AIC.
```

After all models:

```text
12. Compare AIC.
13. Compare likelihood distributions.
14. Select the best-supported model.
15. Perform block bootstrap.
16. Estimate parameter confidence intervals.
```

---

# 43. Recommended directory structure

A clean project layout is:

```text
Demographic_history/
├── SFS/
│   ├── pop_file
│   └── observed_SFS/
├── models/
│   ├── no_geneflow/
│   │   ├── no_geneflow.tpl
│   │   ├── no_geneflow.est
│   │   └── no_geneflow_jointMAFpop1_0.obs
│   ├── ongoing_geneflow/
│   ├── early_geneflow/
│   ├── recent_geneflow/
│   └── diff_geneflow/
├── model_comparison/
├── likelihood_distributions/
├── bootstrap/
└── README.md
```

Large VCF or SFS files do not necessarily need to be stored in a public GitHub repository.

The repository can instead contain:

```text
workflow
model files
population-file template
parameter definitions
small example SFS
helper scripts
summary results
```

---

# 44. Example batch runner for one model

```bash
#!/usr/bin/env bash
set -euo pipefail

PREFIX="early_geneflow"
N_RUNS=100
N_SIM=200000
N_ECM=50

for i in $(seq 1 "$N_RUNS"); do

    RUN_DIR="run${i}"

    mkdir -p "$RUN_DIR"

    cp \
        "${PREFIX}.tpl" \
        "${PREFIX}.est" \
        "${PREFIX}_jointMAFpop1_0.obs" \
        "$RUN_DIR/"

    (
        cd "$RUN_DIR"

        fastsimcoal2 \
            -t "${PREFIX}.tpl" \
            -e "${PREFIX}.est" \
            -m \
            -0 \
            -C 10 \
            -n "$N_SIM" \
            -L "$N_ECM" \
            -s 0 \
            -M \
            -q
    )

done
```

---

# 45. Example extraction of best-run likelihoods

```bash
PREFIX="early_geneflow"

echo -e "Run\tMaxEstLhood" \
    > run_likelihoods.tsv

for i in run*; do

    file="${i}/${PREFIX}/${PREFIX}.bestlhoods"

    if [[ -f "$file" ]]; then

        likelihood=$(awk '
            NR==2 {
                print $NF
            }
        ' "$file")

        echo -e \
            "${i}\t${likelihood}" \
            >> run_likelihoods.tsv

    fi

done
```

Because `.bestlhoods` column order can vary with model specification, always inspect the header before deciding which column contains `MaxEstLhood`.

---

# 46. Example likelihood ranking in R

```r
x <- read.table(
    "run_likelihoods.tsv",
    header = TRUE,
    sep = "\t"
)

x <- x[
    order(
        x$MaxEstLhood,
        decreasing = TRUE
    ),
]

head(
    x,
    10
)
```

The first row is the best optimization replicate if larger likelihood values indicate a better fit.

---

# 47. Example parameter stability plot

After extracting parameter estimates from the top runs:

```r
params <- read.table(
    "top_run_parameters.tsv",
    header = TRUE,
    sep = "\t"
)

boxplot(
    params$N1,
    main = "N1 across top runs",
    ylab = "N1"
)
```

The same approach can be used for:

```text
N2
NANC
TDIV
migration rates
contact times
```

Strong instability among near-optimal runs suggests weak parameter identifiability.

---

# 48. Recommended final outputs

For each model, retain:

```text
best-run parameter estimates
MaxObsLhood
MaxEstLhood
AIC
Delta AIC
Akaike weight
expected SFS
observed-vs-expected SFS plot
model diagram
```

For the selected best model, additionally retain:

```text
likelihood distribution
bootstrap parameter estimates
95% confidence intervals
parameter-converted values in generations/years
```

---

# 49. Suggested summary table

A useful final model-comparison table is:

```text
Model
Number_of_parameters
MaxEstLhood
AIC
Delta_AIC
AIC_weight
```

A useful parameter table for the best model is:

```text
Parameter
Estimate
Bootstrap_median
Lower_95CI
Upper_95CI
Units
```

---

# 50. Important interpretation notes

When interpreting demographic history, remember:

1. The best model is the best among the models tested, not necessarily the true demographic history.
2. More complex models can fit better simply because they contain more parameters; model comparison must account for this.
3. Composite likelihood assumes independence among SFS entries more strongly than may be true for linked genomic data.
4. Parameter estimates can be correlated and weakly identifiable.
5. Population structure not represented in the model can bias estimates.
6. Selection and linked selection can distort the SFS.
7. Errors in ancestral-state inference can bias unfolded SFS analyses.
8. Mutation rate and generation time directly affect scaling of population sizes and event times.
9. Bootstrap uncertainty should be reported for important parameters.
10. Observed-versus-expected SFS fit should be inspected before biological interpretation.

---

# 51. Complete demographic-history workflow

```text
Filtered VCF
     ↓
Population assignment
     ↓
Observed SFS
     ↓
Candidate demographic models
     ↓
.tpl + .est files
     ↓
fastsimcoal2 optimization
     ↓
100+ independent runs/model
     ↓
Best run/model
     ↓
Observed vs expected SFS
     ↓
AIC comparison
     ↓
Likelihood-distribution comparison
     ↓
Best-supported demographic model
     ↓
Block bootstrap
     ↓
Parameter confidence intervals
     ↓
Biological interpretation
```

---

# 52. Suggested manuscript reporting

For reproducibility, report:

```text
fastsimcoal2 version
SFS type: folded or unfolded
projection values
number of populations
number of SNPs/sites
whether invariant sites were included
mutation rate
generation time
candidate demographic models
parameter ranges
number of simulated SFS per ECM cycle
number of ECM cycles
number of independent optimization runs
criterion used to select the best run
model-comparison criterion
number of bootstrap replicates
bootstrap block definition
```

A concise Methods-style description could be:

> Demographic history was inferred from the observed site-frequency spectrum using fastsimcoal2. Alternative demographic scenarios differed in the timing and presence of gene flow following population divergence. For each model, demographic parameters were optimized across multiple independent runs using coalescent simulations of the expected SFS. The run with the highest composite likelihood was retained for each model, and alternative models were compared using AIC and likelihood distributions. Parameter uncertainty for the best-supported model was evaluated using block-bootstrap SFS replicates.

---

# 53. Final checklist

Before considering the demographic analysis complete, confirm:

```text
[ ] SFS population order is correct
[ ] folded/unfolded status is correct
[ ] .tpl and .est prefixes match
[ ] all parameters in .tpl are defined in .est
[ ] historical events are ordered backward in time
[ ] parameter ranges are biologically reasonable
[ ] .simparam file matches the intended model
[ ] at least ~100 optimization runs were performed
[ ] replicate convergence was checked
[ ] best run was selected using MaxEstLhood
[ ] model comparison accounts for parameter number
[ ] observed and expected SFS were visually compared
[ ] likelihood distributions were evaluated
[ ] bootstrap confidence intervals were calculated
[ ] mutation rate and generation time were documented
[ ] final parameters were converted to biologically interpretable units
```
---

# Part II. Sequentially Markovian Coalescent Methods

# 54. Choosing among PSMC, MSMC/MSMC2, SMC++, and fastsimcoal2

The four approaches are complementary rather than interchangeable.

```text
Method          Main input                         Typical purpose
PSMC            One diploid whole genome           Historical Ne of one population
MSMC/MSMC2      Multiple phased haplotypes          Ne and cross-coalescence history
SMC++           Multiple unphased genomes + VCF     Ne and population split inference
fastsimcoal2    Observed SFS                         Explicit demographic model testing
```

A useful integrated strategy is:

```text
PSMC / MSMC / SMC++
        ↓
Explore broad demographic patterns
        ↓
Identify plausible changes in Ne and population separation
        ↓
Construct explicit demographic hypotheses
        ↓
fastsimcoal2
        ↓
Compare alternative demographic models
        ↓
Bootstrap parameter uncertainty
```

The methods should not be expected to produce identical trajectories because they use different summaries of genomic information and have different assumptions and temporal resolution.

---

# Part III. PSMC

# 55. PSMC overview

PSMC (**Pairwise Sequentially Markovian Coalescent**) reconstructs historical effective population size from the distribution of heterozygous sites along a single diploid genome.

The basic workflow is:

```text
Reference genome
      +
Mapped reads
      ↓
BAM
      ↓
Consensus / diploid sequence
      ↓
PSMCFA
      ↓
PSMC
      ↓
Bootstrap
      ↓
Scale by mutation rate and generation time
      ↓
Historical effective population size
```

PSMC is particularly useful when only one or a few high-coverage diploid genomes are available.

Because it uses only two haplotypes at a time, its resolution for very recent demographic history is limited.

---

# 56. PSMC software requirements

Typical software:

```text
samtools
bcftools
bcftools consensus or consensus-generation utilities
PSMC
seqtk
Perl
R
```

Record versions used in the final analysis.

---

# 57. PSMC input requirements

For each individual, start with a high-quality alignment:

```text
sample.bam
reference.fa
```

The BAM should be:

```text
coordinate sorted
indexed
properly filtered
generated against the same reference genome
```

Example:

```bash
samtools index sample.bam
```

Index the reference if necessary:

```bash
samtools faidx reference.fa
```

---

# 58. Generate a diploid consensus sequence

A classical PSMC workflow generates a consensus sequence in which confidently heterozygous positions are retained and low-quality positions are masked.

The exact consensus-generation command depends on the `samtools`/`bcftools` version and the filtering strategy used in the study.

A modern starting point is to call variants:

```bash
bcftools mpileup \
    -Ou \
    -f reference.fa \
    sample.bam | \
bcftools call \
    -mv \
    -Oz \
    -o sample.vcf.gz
```

Index:

```bash
bcftools index sample.vcf.gz
```

For PSMC, however, it is important that the final diploid sequence also represents confidently callable homozygous positions and masks poorly covered or low-quality regions. Therefore, a callable-site mask should be generated and applied rather than simply converting variant sites alone.

---

# 59. Define callable regions

Coverage filters should remove positions with extremely low or extremely high depth.

Calculate depth:

```bash
samtools depth \
    -a \
    sample.bam \
    > sample.depth.txt
```

Estimate the mean or median depth and choose study-specific thresholds.

Conceptually:

```text
minimum depth = approximately 1/3 of average depth
maximum depth = approximately 2× average depth
```

These values are examples, not universal thresholds.

Create a BED file of callable positions or intervals after applying:

```text
depth filters
mapping-quality filters
base-quality filters
repeat/mappability filters
```

The same callable-region logic should be applied consistently among individuals.

---

# 60. Convert the diploid consensus to PSMCFA

Once a suitable diploid consensus FASTQ is available:

```bash
fq2psmcfa \
    -q20 \
    sample.diploid.fq.gz \
    > sample.psmcfa
```

Here:

```text
-q20
```

retains bases with consensus quality of at least 20.

The exact threshold should be documented.

---

# 61. Run PSMC

A commonly used PSMC command structure is:

```bash
psmc \
    -N25 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o sample.psmc \
    sample.psmcfa
```

Important options include:

```text
-N    maximum number of iterations
-t    initial θ/ρ ratio-related parameter
-r    initial recombination-to-mutation ratio
-p    atomic time interval pattern
```

The interval pattern:

```text
4+25*2+4+6
```

is widely used as a starting configuration, but it should not be treated as universally optimal.

---

# 62. Inspect PSMC convergence

Inspect the output:

```bash
less sample.psmc
```

Check whether parameter estimates stabilize across iterations.

If necessary, compare alternative interval patterns or starting parameters.

The final inference should not rely on a clearly unconverged run.

---

# 63. PSMC bootstrap analysis

PSMC includes utilities for bootstrap resampling.

Split the PSMCFA sequence into blocks:

```bash
splitfa \
    sample.psmcfa \
    > sample.split.psmcfa
```

Run bootstrap replicates:

```bash
mkdir -p bootstrap

for i in $(seq 1 100); do

    psmc \
        -N25 \
        -t15 \
        -r5 \
        -b \
        -p "4+25*2+4+6" \
        -o "bootstrap/sample.bootstrap.${i}.psmc" \
        sample.split.psmcfa

done
```

The `-b` option performs bootstrap resampling.

A typical final analysis uses approximately 100 bootstrap replicates.

---

# 64. Combine PSMC bootstrap results

Combine the original run and bootstrap runs:

```bash
cat \
    sample.psmc \
    bootstrap/sample.bootstrap.*.psmc \
    > sample.combined.psmc
```

This combined file can be used for plotting.

---

# 65. Scale PSMC results

PSMC output must be scaled using:

```text
mutation rate per site per generation (μ)
generation time (g)
```

A typical plotting command is:

```bash
psmc_plot.pl \
    -u MUTATION_RATE \
    -g GENERATION_TIME \
    sample \
    sample.combined.psmc
```

Replace:

```text
MUTATION_RATE
GENERATION_TIME
```

with values appropriate for the study species.

Do not copy mutation rates from unrelated species without justification.

---

# 66. PSMC interpretation

The final trajectory is usually interpreted as:

```text
x-axis = time before present
y-axis = effective population size (Ne)
```

Potential signals include:

```text
population expansion
population contraction
bottleneck
long-term decline
historical increase in Ne
```

PSMC trajectories should be interpreted as smoothed coalescent histories rather than exact census population-size histories.

---

# 67. Run PSMC for multiple individuals

When several high-coverage individuals are available, run PSMC independently for each one.

Example:

```bash
for sample in sample1 sample2 sample3 sample4; do

    psmc \
        -N25 \
        -t15 \
        -r5 \
        -p "4+25*2+4+6" \
        -o "${sample}.psmc" \
        "${sample}.psmcfa"

done
```

Overlaying individuals from the same population can reveal whether inferred histories are consistent.

---

# Part IV. MSMC and MSMC2

# 68. MSMC/MSMC2 overview

MSMC (**Multiple Sequentially Markovian Coalescent**) extends pairwise SMC approaches to multiple haplotypes.

MSMC2 is generally preferred for many modern analyses.

It can be used to infer:

```text
historical effective population size
coalescence-rate changes
population separation
relative cross-coalescence rate
```

A conceptual workflow is:

```text
Phased VCF
    +
Callable masks
    ↓
MSMC input files
    ↓
MSMC2 within-population analysis
    ↓
MSMC2 cross-population analysis
    ↓
combineCrossCoal.py
    ↓
Relative cross-coalescence rate
    ↓
Scale by mutation rate and generation time
```

---

# 69. MSMC2 software requirements

Typical requirements:

```text
bcftools
MSMC-tools
MSMC2
bgzip
tabix
Python
R
```

For population-separation analyses, phased genotypes are strongly preferred because MSMC uses haplotype information.

---

# 70. Prepare phased VCF files

A typical input is:

```text
phased.vcf.gz
```

Restrict to high-quality biallelic SNPs if necessary:

```bash
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    phased.vcf.gz \
    -Oz \
    -o phased.biallelic.vcf.gz
```

Index:

```bash
bcftools index \
    phased.biallelic.vcf.gz
```

---

# 71. Prepare callable masks for MSMC2

Each individual should have a callable-region mask.

Conceptually:

```text
sample1.mask.bed.gz
sample2.mask.bed.gz
sample3.mask.bed.gz
...
```

Masks should account for:

```text
coverage
mapping quality
base quality
repetitive sequence
poorly mappable sequence
assembly gaps
```

Using a mask prevents non-callable regions from being interpreted as long homozygous segments.

---

# 72. Generate MSMC input with generate_multihetsep.py

MSMC-tools provides:

```text
generate_multihetsep.py
```

A general command is:

```bash
generate_multihetsep.py \
    --mask sample1.mask.bed.gz \
    --mask sample2.mask.bed.gz \
    sample1.vcf.gz \
    sample2.vcf.gz \
    > chr01.multihetsep.txt
```

For chromosome-wise analysis:

```bash
for chr in chr01 chr02 chr03 chr04; do

    generate_multihetsep.py \
        --mask "${chr}.sample1.mask.bed.gz" \
        --mask "${chr}.sample2.mask.bed.gz" \
        "${chr}.sample1.vcf.gz" \
        "${chr}.sample2.vcf.gz" \
        > "${chr}.multihetsep.txt"

done
```

The exact input organization can be adapted to whether genotypes are stored in one multisample VCF or separate per-sample files.

---

# 73. Understand haplotype indices

MSMC2 uses haplotype indices.

For two diploid individuals:

```text
individual 1 → haplotypes 0,1
individual 2 → haplotypes 2,3
```

For four diploid individuals:

```text
individual 1 → 0,1
individual 2 → 2,3
individual 3 → 4,5
individual 4 → 6,7
```

The `-I` option determines which haplotypes are analyzed.

Carefully verify the sample and haplotype order before running cross-population analyses.

---

# 74. MSMC2 within-population Ne analysis

For one diploid individual represented by haplotypes `0,1`:

```bash
msmc2 \
    -I 0,1 \
    -o population1 \
    chr*.multihetsep.txt
```

For multiple haplotypes from one population:

```bash
msmc2 \
    -I 0,1,2,3 \
    -o population1 \
    chr*.multihetsep.txt
```

The exact number of haplotypes depends on the study design and computational constraints.

---

# 75. MSMC2 analysis for population 2

Example:

```bash
msmc2 \
    -I 4,5,6,7 \
    -o population2 \
    chr*.multihetsep.txt
```

This estimates within-population coalescence rates for the second population.

---

# 76. MSMC2 cross-population analysis

To estimate cross-population coalescence, specify cross-population haplotype pairs.

For example:

```bash
msmc2 \
    -I 0-4,0-5,1-4,1-5 \
    -o pop1_pop2_cross \
    chr*.multihetsep.txt
```

The exact pair specification must match the haplotype order in the input.

For larger sample sets, construct the cross-population index list systematically rather than manually.

---

# 77. Combine cross-coalescence results

MSMC-tools provides utilities such as:

```text
combineCrossCoal.py
```

A typical conceptual command is:

```bash
combineCrossCoal.py \
    population1.final.txt \
    population2.final.txt \
    pop1_pop2_cross.final.txt \
    > pop1_pop2.combined.txt
```

This combines within- and between-population coalescence estimates for relative cross-coalescence analysis.

Check the exact script interface supplied with the installed MSMC-tools version before running.

---

# 78. Relative cross-coalescence rate

The relative cross-coalescence rate (RCCR) is often interpreted on a scale from approximately:

```text
1 → populations behaving as one panmictic population
0 → complete population separation
```

The transition is usually gradual.

A value such as:

```text
RCCR = 0.5
```

is sometimes used as a descriptive midpoint for population separation, but it should not automatically be interpreted as an exact divergence time.

---

# 79. Scale MSMC2 time and Ne

MSMC2 output is scaled by the mutation rate.

For each time estimate:

```text
time in generations
=
scaled time / mutation rate
```

Then:

```text
time in years
=
time in generations × generation time
```

Effective population size is similarly scaled using the mutation rate according to the output definition.

Use the exact scaling equations corresponding to the MSMC/MSMC2 output being plotted.

---

# 80. MSMC2 bootstrap strategy

Bootstrap analysis can be performed by resampling genomic segments.

A conceptual workflow is:

```text
multihetsep input
      ↓
split genome into blocks
      ↓
sample blocks with replacement
      ↓
construct bootstrap input
      ↓
run MSMC2
      ↓
repeat 50–100 times
      ↓
confidence envelope
```

MSMC-tools includes utilities for generating bootstrap datasets.

Bootstrap block size should be sufficiently large to preserve local linkage structure.

---

# 81. MSMC2 interpretation

Within-population MSMC2 provides a historical coalescence trajectory that can be converted to effective population size.

Cross-population MSMC2 provides information about the gradual loss of shared ancestry.

Use it to explore:

```text
population separation
prolonged gene flow
secondary contact signatures
asynchronous divergence
```

However, RCCR is not itself a direct migration-rate estimate.

---

# Part V. SMC++

# 82. SMC++ overview

SMC++ combines ideas from the SFS and sequentially Markovian coalescent framework.

Important advantages include:

```text
uses multiple individuals
does not require fully phased genomes
can infer population-size history
can estimate population split histories
generally has better recent-time resolution than PSMC
```

The basic workflow is:

```text
Filtered multisample VCF
      +
Callable/masked genomic regions
      ↓
smc++ vcf2smc
      ↓
.smc.gz files
      ↓
smc++ estimate
      ↓
population-size history
      ↓
smc++ split
      ↓
population split inference
      ↓
smc++ plot
```

---

# 83. SMC++ software requirements

Required software:

```text
SMC++
bcftools
bgzip
tabix
```

Check installation:

```bash
smc++ --help
```

---

# 84. Prepare the SMC++ VCF

Start with a high-quality multisample VCF:

```text
demographic.vcf.gz
```

Restrict to biallelic SNPs if appropriate:

```bash
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    demographic.vcf.gz \
    -Oz \
    -o demographic.biallelic.vcf.gz
```

Index:

```bash
bcftools index \
    demographic.biallelic.vcf.gz
```

---

# 85. Population sample lists

Prepare one sample list per population.

Example:

```text
pop1.samples.txt
pop2.samples.txt
```

Each contains one sample ID per line.

Example:

```text
sample01
sample02
sample03
sample04
```

Check that every sample occurs in the VCF:

```bash
bcftools query -l \
    demographic.biallelic.vcf.gz
```

---

# 86. Distinguished individuals

`smc++ vcf2smc` uses one or more distinguished individuals to contribute long-range coalescent information.

A population specification has the form:

```text
Population:sample1,sample2,...
```

For example:

```text
Pop1:sample01,sample02
```

The exact selection of distinguished samples should be documented.

Running alternative distinguished individuals can be useful as a sensitivity check.

---

# 87. Convert VCF to SMC format

Create an output directory:

```bash
mkdir -p smc
```

For one chromosome:

```bash
smc++ vcf2smc \
    demographic.biallelic.vcf.gz \
    smc/chr01.Pop1.smc.gz \
    chr01 \
    Pop1:sample01,sample02
```

Repeat chromosome by chromosome.

---

# 88. Convert all chromosomes for one population

Example:

```bash
for chr in chr01 chr02 chr03 chr04 chr05; do

    smc++ vcf2smc \
        demographic.biallelic.vcf.gz \
        "smc/${chr}.Pop1.smc.gz" \
        "$chr" \
        Pop1:sample01,sample02

done
```

Use chromosome names exactly as they appear in the VCF.

Check with:

```bash
bcftools query \
    -f '%CHROM\n' \
    demographic.biallelic.vcf.gz | \
sort -u
```

---

# 89. Include all population samples

A more complete population specification can be generated from the sample list.

For example:

```bash
POP1_SAMPLES=$(paste -sd, pop1.samples.txt)
```

Then:

```bash
smc++ vcf2smc \
    demographic.biallelic.vcf.gz \
    smc/chr01.Pop1.smc.gz \
    chr01 \
    "Pop1:${POP1_SAMPLES}"
```

Depending on the intended SMC++ setup, distinguished individuals and the full sample set should be configured according to the installed SMC++ version and study design.

---

# 90. Mask inaccessible genomic regions

Regions that cannot be confidently genotyped should not contribute as ordinary homozygous sequence.

Prepare a mask such as:

```text
mask.bed.gz
```

Possible excluded regions include:

```text
assembly gaps
low coverage
excessive coverage
low mappability
repetitive sequence
poor mapping-quality regions
```

When using SMC++, ensure that missing and inaccessible regions are represented consistently in the VCF/SMC conversion workflow.

---

# 91. Estimate population-size history with SMC++

Assume a mutation rate:

```bash
MU="6.96e-9"
```

Run:

```bash
smc++ estimate \
    "$MU" \
    smc/chr*.Pop1.smc.gz \
    --cores 10 \
    -o Pop1_estimate
```

The mutation rate should be justified for the study organism.

---

# 92. Estimate another population

```bash
smc++ estimate \
    "$MU" \
    smc/chr*.Pop2.smc.gz \
    --cores 10 \
    -o Pop2_estimate
```

Use the same mutation rate and comparable filtering when populations are compared.

---

# 93. Set time boundaries

SMC++ allows the inference interval to be restricted.

A conceptual command is:

```bash
smc++ estimate \
    "$MU" \
    --timepoints TMIN TMAX \
    smc/chr*.Pop1.smc.gz \
    -o Pop1_estimate
```

Choose time boundaries appropriate for:

```text
sample size
mutation rate
genome length
expected demographic history
```

Do not overinterpret periods outside the well-supported temporal range.

---

# 94. Plot SMC++ results

A basic plot can be generated with:

```bash
smc++ plot \
    Pop1.pdf \
    Pop1_estimate/model.final.json
```

For biological time in years, specify generation time using the plotting options supported by the installed SMC++ version.

The final plot should clearly report:

```text
mutation rate
generation time
population
time units
Ne units
```

---

# 95. Compare multiple populations

If multiple model files are available, they can be plotted together.

Conceptually:

```bash
smc++ plot \
    populations.pdf \
    Pop1_estimate/model.final.json \
    Pop2_estimate/model.final.json
```

This makes it possible to compare population-specific Ne trajectories.

---

# 96. SMC++ split analysis

SMC++ can estimate population split histories using two populations.

First generate SMC files containing the required population information.

Then estimate marginal population histories before fitting the split model.

Conceptually:

```bash
smc++ split \
    -o split_output \
    Pop1_estimate/model.final.json \
    Pop2_estimate/model.final.json \
    smc/*.smc.gz
```

The exact argument order and options should be verified against the installed SMC++ version because the CLI can differ among releases.

---

# 97. Interpreting SMC++ split estimates

A split estimate describes population separation under the fitted model.

Interpret it together with:

```text
population structure
gene flow
MSMC2 RCCR
fastsimcoal2 migration models
```

If populations experienced prolonged migration after divergence, a single split-time estimate may summarize a much more complex process.

---

# 98. SMC++ cross-validation

SMC++ supports cross-validation/model-complexity control in suitable workflows.

This can help avoid overfitting the Ne trajectory.

A practical strategy is to compare results across:

```text
regularization settings
time-point settings
distinguished individuals
chromosome subsets
```

Stable major features are more convincing than small fluctuations present in only one configuration.

---

# 99. SMC++ bootstrap or chromosome resampling

Uncertainty can be explored by resampling chromosomes or large genomic blocks.

Conceptually:

```text
original chromosomes
      ↓
resample chromosomes/blocks
      ↓
generate bootstrap SMC input
      ↓
smc++ estimate
      ↓
repeat
      ↓
Ne confidence envelope
```

For organisms with relatively few chromosomes, block-based resampling may provide more bootstrap units than chromosome-only resampling.

---

# Part VI. Integrated Demographic Analysis

# 100. Recommended combined workflow

A comprehensive demographic-history analysis can use all four methods:

```text
                    Whole-genome variation
                           │
          ┌────────────────┼────────────────┐
          │                │                │
        PSMC          MSMC/MSMC2          SMC++
          │                │                │
          │                │                │
      Historical Ne   Historical Ne     Historical Ne
                           +
                     Cross-coalescence
                           │
          └────────────────┼────────────────┘
                           │
                 Exploratory demographic
                       patterns
                           │
                           ↓
                    Candidate models
                           │
                           ↓
                    fastsimcoal2
                           │
             ┌─────────────┴─────────────┐
             │                           │
        Model comparison           Parameter fitting
             │                           │
             └─────────────┬─────────────┘
                           │
                      Best model
                           │
                           ↓
                    Block bootstrap
                           │
                           ↓
                  Confidence intervals
```

---

# 101. What each method contributes

## PSMC

Best suited for:

```text
long-term Ne trajectory
single high-coverage diploid genomes
population-specific historical patterns
```

Main limitation:

```text
weak recent-time resolution
only two haplotypes per analysis
```

---

## MSMC/MSMC2

Best suited for:

```text
multiple haplotypes
population separation
relative cross-coalescence
Ne trajectories
```

Important consideration:

```text
phasing quality
callable masks
haplotype indexing
```

---

## SMC++

Best suited for:

```text
larger sample sizes
unphased genomic data
recent and ancient Ne inference
population split estimation
```

Important consideration:

```text
mutation-rate scaling
masking
distinguished-individual choice
```

---

## fastsimcoal2

Best suited for:

```text
explicit demographic hypotheses
migration models
population-size parameters
divergence-time parameters
formal model comparison
```

Important consideration:

```text
candidate models must be explicitly specified
results are conditional on the models tested
```

---

# 102. Recommended order of analysis

For a population-genomics project, a practical order is:

```text
Step 1
Population structure / PCA / ADMIXTURE / phylogeny

Step 2
Population diversity and differentiation

Step 3
PSMC
Individual-level historical Ne

Step 4
MSMC2
Population Ne + relative cross-coalescence

Step 5
SMC++
Population Ne + split history using larger sample sets

Step 6
SFS construction with easySFS

Step 7
fastsimcoal2
Explicit demographic model testing

Step 8
Model comparison

Step 9
Bootstrap uncertainty

Step 10
Integrate results across methods
```

---

# 103. Consistent filtering across methods

Although the exact input format differs, preprocessing should be as consistent as possible.

Record:

```text
reference genome
mapping filters
base-quality filters
genotype-quality filters
depth filters
repeat masking
mappability masking
chromosomes included
sex chromosomes included/excluded
structural variants included/excluded
```

Different callable regions can create apparent differences among demographic trajectories.

---

# 104. Mutation rate

All time-scaled demographic methods depend strongly on the assumed mutation rate.

Record:

```text
μ = mutations per site per generation
```

Use the same biologically justified mutation rate across methods whenever possible.

If different mutation rates are used, comparisons must be rescaled before interpretation.

---

# 105. Generation time

Convert generations to years using:

```text
years = generations × generation time
```

Record:

```text
g = years per generation
```

For annual crops:

```text
g may often be approximately one year
```

but this should be justified from the biology and sampling context rather than assumed automatically.

---

# 106. Recombination rate

PSMC/MSMC/SMC++ infer demographic history from linkage/coalescence patterns and therefore interact with assumptions about recombination.

When a recombination rate is required or assumed, document:

```text
rate per base per generation
source of the estimate
whether a constant or variable rate was used
```

A genome-wide average can be a practical approximation, but recombination heterogeneity can affect local genealogical patterns.

---

# 107. Neutral-region filtering

For demographic inference, consider whether to restrict analysis to approximately neutral genomic regions.

Potential exclusions include:

```text
coding sequence
strongly selected regions
centromeres
large inversions
low-recombination regions
poorly mappable sequence
repetitive sequence
```

The choice depends on the biological question.

Document the retained genomic fraction and use compatible callable lengths when scaling SFS-based analyses.

---

# 108. Compare Ne trajectories across methods

After scaling all methods with the same mutation rate and generation time, compare:

```text
PSMC Ne
MSMC2 Ne
SMC++ Ne
```

Focus on broad features such as:

```text
major expansion
major contraction
bottleneck timing
long-term decline
relative differences among populations
```

Do not expect fine-scale curves to overlap exactly.

---

# 109. Compare divergence estimates

Population divergence can be summarized using:

```text
MSMC2 RCCR transition
SMC++ split estimate
fastsimcoal2 TDIV
```

These values measure related but not identical aspects of separation.

A useful final figure can show:

```text
MSMC2 RCCR curve
SMC++ split estimate
fastsimcoal2 divergence estimate + bootstrap CI
```

This makes methodological differences transparent.
