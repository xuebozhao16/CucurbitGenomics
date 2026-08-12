# GWAS and Genomic Selection Pipeline

## Overview

This document describes a reproducible workflow for genome-wide association studies (GWAS) and genomic selection (GS) using SNP and INDEL variation.

The workflow includes:

1. merging chromosome-level VCF files;
2. separating SNPs and INDELs;
3. converting VCF files to PLINK format;
4. principal component analysis;
5. genomic relationship matrix construction;
6. phenotype formatting;
7. SNP- and INDEL-based heritability estimation;
8. GWAS using GCTA-MLMA;
9. genomic selection using genomic best linear unbiased prediction (GBLUP);
10. cross-validation and prediction accuracy evaluation.

All file paths below are relative examples. Replace file names and sample lists as needed for the target dataset.

---

## Software requirements

The workflow requires:

- `bcftools`
- `PLINK`
- `GCTA`
- `R`
- R packages:
  - `data.table`
  - `rrBLUP`
  - `ggplot2`

Optional software for downstream visualization:

- `CMplot`
- `qqman`

---

# Part I. Variant preparation

## 1. Merge chromosome-level VCF files

Assume chromosome-level VCF files are named:

```text
chr1_filtered.vcf.gz
chr2_filtered.vcf.gz
...
chr12_filtered.vcf.gz
```

Index all input VCF files:

```bash
for file in chr*_filtered.vcf.gz; do
    bcftools index -f "$file"
done
```

Concatenate chromosomes in genomic order:

```bash
bcftools concat \
    -Oz \
    -o merged.vcf.gz \
    chr{1..12}_filtered.vcf.gz
```

Index the merged VCF:

```bash
bcftools index -f merged.vcf.gz
```

---

## 2. Separate SNPs and INDELs

Extract SNPs:

```bash
bcftools view \
    -v snps \
    merged.vcf.gz \
    -Oz \
    -o merged_snps.vcf.gz

bcftools index -f merged_snps.vcf.gz
```

Extract INDELs:

```bash
bcftools view \
    -v indels \
    merged.vcf.gz \
    -Oz \
    -o merged_indels.vcf.gz

bcftools index -f merged_indels.vcf.gz
```

---

## 3. Convert VCF files to PLINK format

Convert SNPs:

```bash
plink \
    --vcf merged_snps.vcf.gz \
    --make-bed \
    --out SNP
```

Convert INDELs:

```bash
plink \
    --vcf merged_indels.vcf.gz \
    --make-bed \
    --out INDEL
```

The output files are:

```text
SNP.bed
SNP.bim
SNP.fam

INDEL.bed
INDEL.bim
INDEL.fam
```

---

# Part II. Population structure correction

## 4. Principal component analysis

PCA can be calculated separately for SNPs and INDELs.

### SNP PCA

```bash
plink \
    --bfile SNP \
    --pca 20 header tabs \
    --out SNP_PCA
```

### INDEL PCA

```bash
plink \
    --bfile INDEL \
    --pca 20 header tabs \
    --out INDEL_PCA
```

The principal outputs are:

```text
SNP_PCA.eigenval
SNP_PCA.eigenvec

INDEL_PCA.eigenval
INDEL_PCA.eigenvec
```

---

## 5. Prepare PCA covariates for GCTA

The number of PCs included in GWAS should be determined based on population structure and model diagnostics.

The following R script extracts the first three PCs.

```r
pca <- read.table(
    "SNP_PCA.eigenvec",
    header = TRUE,
    stringsAsFactors = FALSE
)

# PLINK eigenvec files contain:
# FID, IID, PC1, PC2, ...

pca_result <- pca[, c(1, 2, 3, 4, 5)]

colnames(pca_result) <- c(
    "FID",
    "IID",
    "PC1",
    "PC2",
    "PC3"
)

write.table(
    pca_result,
    "SNP_PCA_covariates.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
```

If more PCs are required, simply retain additional PC columns.

---

# Part III. Genomic relationship matrices

## 6. Build genomic relationship matrices with GCTA

### SNP-based GRM

```bash
gcta64 \
    --bfile SNP \
    --make-grm \
    --out SNP_GRM \
    --thread-num 30
```

### INDEL-based GRM

```bash
gcta64 \
    --bfile INDEL \
    --make-grm \
    --out INDEL_GRM \
    --thread-num 30
```

The principal GRM files are:

```text
SNP_GRM.grm.bin
SNP_GRM.grm.N.bin
SNP_GRM.grm.id

INDEL_GRM.grm.bin
INDEL_GRM.grm.N.bin
INDEL_GRM.grm.id
```

These GRMs can be used for:

- SNP-based heritability estimation;
- INDEL-based heritability estimation;
- mixed-model GWAS;
- genomic prediction.

---

# Part IV. Phenotype preparation

## 7. Format phenotype files

Each phenotype file is assumed to contain:

```text
Sample  Trait
sample1 10.2
sample2 11.5
sample3  9.8
```

GCTA expects three columns:

```text
FID IID phenotype
```

If FID and IID are identical, convert phenotype files using:

```bash
mkdir -p phenotype_gcta

for file in phenotype/*.txt; do

    filename=$(basename "$file")
    trait="${filename%.*}"

    awk 'NR>1 {
        print $1, $1, $2
    }' OFS="\t" "$file" \
        > "phenotype_gcta/${trait}.txt"

done
```

Example output:

```text
sample1 sample1 10.2
sample2 sample2 11.5
sample3 sample3 9.8
```

Missing phenotype values should be coded consistently and excluded before analysis.

---

# Part V. Heritability estimation

## 8. Estimate SNP- and INDEL-based heritability

### SNP-based heritability

```bash
gcta64 \
    --grm SNP_GRM \
    --pheno phenotype_gcta/TRAIT.txt \
    --reml \
    --out TRAIT_SNP_heritability \
    --thread-num 30
```

### INDEL-based heritability

```bash
gcta64 \
    --grm INDEL_GRM \
    --pheno phenotype_gcta/TRAIT.txt \
    --reml \
    --out TRAIT_INDEL_heritability \
    --thread-num 30
```

The `.hsq` output contains estimates such as:

```text
V(G)
V(e)
Vp
V(G)/Vp
logL
```

`V(G)/Vp` represents the proportion of phenotypic variance explained by the corresponding marker set.

---

# Part VI. Genome-wide association study

## 9. SNP-based GWAS with GCTA-MLMA

GCTA-MLMA fits each marker as a fixed effect while accounting for the genomic relationship matrix.

```bash
gcta64 \
    --bfile SNP \
    --grm SNP_GRM \
    --pheno phenotype_gcta/TRAIT.txt \
    --qcovar SNP_PCA_covariates.txt \
    --mlma \
    --out TRAIT_SNP_GWAS \
    --thread-num 30
```

The principal output is:

```text
TRAIT_SNP_GWAS.mlma
```

Typical columns include:

```text
Chr
SNP
bp
A1
A2
Freq
b
se
p
```

---

## 10. INDEL-based GWAS

```bash
gcta64 \
    --bfile INDEL \
    --grm INDEL_GRM \
    --pheno phenotype_gcta/TRAIT.txt \
    --qcovar INDEL_PCA_covariates.txt \
    --mlma \
    --out TRAIT_INDEL_GWAS \
    --thread-num 30
```

---

## 11. Run GWAS for multiple traits

A simple Bash loop can automate the analysis.

```bash
mkdir -p GWAS_results
mkdir -p heritability

for pheno in phenotype_gcta/*.txt; do

    trait=$(basename "$pheno" .txt)

    # SNP-based heritability
    gcta64 \
        --grm SNP_GRM \
        --pheno "$pheno" \
        --reml \
        --out "heritability/${trait}_SNP" \
        --thread-num 30

    # INDEL-based heritability
    gcta64 \
        --grm INDEL_GRM \
        --pheno "$pheno" \
        --reml \
        --out "heritability/${trait}_INDEL" \
        --thread-num 30

    # SNP GWAS
    gcta64 \
        --bfile SNP \
        --grm SNP_GRM \
        --pheno "$pheno" \
        --qcovar SNP_PCA_covariates.txt \
        --mlma \
        --out "GWAS_results/${trait}_SNP" \
        --thread-num 30

    # INDEL GWAS
    gcta64 \
        --bfile INDEL \
        --grm INDEL_GRM \
        --pheno "$pheno" \
        --qcovar INDEL_PCA_covariates.txt \
        --mlma \
        --out "GWAS_results/${trait}_INDEL" \
        --thread-num 30

done
```

---

# Part VII. GWAS significance thresholds

## 12. Multiple-testing correction

A simple Bonferroni threshold can be calculated as:

```text
P_threshold = 0.05 / number_of_markers
```

For example, if 1,000,000 markers are tested:

```text
0.05 / 1,000,000 = 5 × 10^-8
```

A suggestive threshold can also be defined as:

```text
1 / number_of_markers
```

When LD among markers is strong, an effective-number-of-independent-markers threshold may be more appropriate than a strict Bonferroni threshold.

---

# Part VIII. Genomic selection

## 13. Overview of genomic selection

Genomic selection differs from GWAS in its objective.

GWAS aims to identify individual variants associated with a phenotype, whereas genomic selection uses genome-wide marker information simultaneously to predict the genetic value of individuals.

The basic GS workflow is:

```text
Genotype matrix
      +
Phenotype
      ↓
Training population
      ↓
Genomic prediction model
      ↓
Genomic estimated breeding values (GEBVs)
      ↓
Prediction accuracy
```

A standard baseline model is GBLUP.

---

# Part IX. Prepare genotype data for genomic selection

## 14. Export genotype matrix from PLINK

The additive genotype matrix can be exported using PLINK.

```bash
plink \
    --bfile SNP \
    --recode A \
    --out GS_genotype
```

This generates:

```text
GS_genotype.raw
```

Marker genotypes are generally coded as:

```text
0 = homozygous reference
1 = heterozygous
2 = homozygous alternative
```

The first columns contain sample metadata and should be removed before fitting the model.

---

## 15. Prepare phenotype data for GS

A simple phenotype table can be formatted as:

```text
Sample Trait
sample1 10.2
sample2 11.5
sample3  9.8
```

The sample order does not need to match the genotype file, but sample IDs must be identical so that genotype and phenotype records can be matched correctly.

---

# Part X. GBLUP genomic selection using rrBLUP

## 16. Read and prepare genotype data

```r
library(data.table)
library(rrBLUP)

geno_raw <- fread(
    "GS_genotype.raw",
    data.table = FALSE
)

# PLINK .raw files normally contain:
# FID IID PAT MAT SEX PHENOTYPE marker1 marker2 ...

sample_id <- geno_raw$IID

geno <- geno_raw[, 7:ncol(geno_raw)]

rownames(geno) <- sample_id

geno <- as.matrix(geno)
```

---

## 17. Marker quality control for genomic selection

Remove markers with no variation:

```r
keep <- apply(
    geno,
    2,
    function(x) var(x, na.rm = TRUE) > 0
)

geno <- geno[, keep]
```

Markers with excessive missing data can also be removed:

```r
missing_rate <- apply(
    geno,
    2,
    function(x) mean(is.na(x))
)

geno <- geno[, missing_rate <= 0.2]
```

---

## 18. Impute missing genotypes

A simple mean imputation can be used for a baseline GS analysis.

```r
for (j in seq_len(ncol(geno))) {

    marker_mean <- mean(
        geno[, j],
        na.rm = TRUE
    )

    geno[is.na(geno[, j]), j] <- marker_mean
}
```

For final analyses, more sophisticated genotype imputation may be preferable when missingness is substantial.

---

## 19. Build the additive genomic relationship matrix

The `A.mat()` function in `rrBLUP` constructs an additive genomic relationship matrix.

```r
G <- A.mat(
    geno,
    impute.method = "mean",
    return.imputed = FALSE
)
```

The matrix dimensions are:

```text
number_of_individuals × number_of_individuals
```

---

# Part XI. Match phenotype and genotype data

## 20. Read phenotype data

```r
pheno <- read.table(
    "phenotype.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
```

Assume the phenotype table contains:

```text
Sample Trait
```

Retain individuals present in both genotype and phenotype datasets.

```r
common <- intersect(
    rownames(G),
    pheno$Sample
)

G <- G[common, common]

pheno <- pheno[
    match(common, pheno$Sample),
]
```

Extract phenotype:

```r
y <- pheno$Trait
names(y) <- pheno$Sample
```

---

# Part XII. Cross-validation for genomic selection

## 21. Five-fold cross-validation

A standard approach is repeated k-fold cross-validation.

The following example performs five-fold cross-validation.

```r
library(rrBLUP)

set.seed(2026)

n <- length(y)
k <- 5

fold <- sample(
    rep(1:k, length.out = n)
)

prediction <- rep(NA, n)

for (f in 1:k) {

    test_index <- which(fold == f)

    y_train <- y
    y_train[test_index] <- NA

    fit <- kin.blup(
        data = data.frame(
            ID = names(y_train),
            phenotype = y_train
        ),
        geno = "ID",
        pheno = "phenotype",
        K = G
    )

    prediction[test_index] <- fit$g[
        match(
            names(y)[test_index],
            names(fit$g)
        )
    ]
}
```

---

## 22. Calculate prediction ability

Prediction ability is commonly calculated as the Pearson correlation between predicted genomic breeding values and observed phenotypes.

```r
prediction_ability <- cor(
    prediction,
    y,
    use = "complete.obs"
)

prediction_ability
```

A higher correlation indicates greater predictive performance.

---

## 23. Repeated cross-validation

Single cross-validation partitions may be sensitive to random sampling. Repeated cross-validation provides a more stable estimate.

```r
library(rrBLUP)

set.seed(2026)

n_rep <- 100
k <- 5

accuracy <- numeric(n_rep)

for (rep in 1:n_rep) {

    fold <- sample(
        rep(1:k, length.out = length(y))
    )

    prediction <- rep(
        NA,
        length(y)
    )

    for (f in 1:k) {

        test_index <- which(
            fold == f
        )

        y_train <- y
        y_train[test_index] <- NA

        fit <- kin.blup(
            data = data.frame(
                ID = names(y_train),
                phenotype = y_train
            ),
            geno = "ID",
            pheno = "phenotype",
            K = G
        )

        prediction[test_index] <- fit$g[
            match(
                names(y)[test_index],
                names(fit$g)
            )
        ]
    }

    accuracy[rep] <- cor(
        prediction,
        y,
        use = "complete.obs"
    )
}

mean(accuracy)
sd(accuracy)
```

Report both the mean and standard deviation of prediction ability across replicates.

---

# Part XIII. Genomic estimated breeding values

## 24. Fit GBLUP using all phenotyped individuals

After evaluating prediction performance, fit the final model using all available phenotyped individuals.

```r
fit_final <- kin.blup(
    data = data.frame(
        ID = names(y),
        phenotype = y
    ),
    geno = "ID",
    pheno = "phenotype",
    K = G
)
```

Extract genomic estimated breeding values:

```r
gebv <- data.frame(
    Sample = names(fit_final$g),
    GEBV = fit_final$g
)

gebv <- gebv[
    order(
        gebv$GEBV,
        decreasing = TRUE
    ),
]

write.table(
    gebv,
    "GEBV.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

---

# Part XIV. Training and prediction populations

## 25. Predict unphenotyped individuals

A typical breeding application includes:

- a **training population** with genotype and phenotype data;
- a **prediction population** with genotype data but no phenotype data.

The relationship matrix should include both populations.

Phenotypes of prediction individuals are coded as missing:

```r
phenotype_all <- data.frame(
    ID = rownames(G),
    phenotype = NA_real_
)

phenotype_all$phenotype[
    match(
        pheno$Sample,
        phenotype_all$ID
    )
] <- pheno$Trait
```

Fit GBLUP:

```r
fit_prediction <- kin.blup(
    data = phenotype_all,
    geno = "ID",
    pheno = "phenotype",
    K = G
)
```

Extract predicted breeding values for unphenotyped individuals:

```r
prediction_ids <- phenotype_all$ID[
    is.na(
        phenotype_all$phenotype
    )
]

predicted_GEBV <- data.frame(
    Sample = prediction_ids,
    GEBV = fit_prediction$g[
        match(
            prediction_ids,
            names(fit_prediction$g)
        )
    ]
)

write.table(
    predicted_GEBV,
    "Prediction_population_GEBV.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

---

# Part XV. SNP versus INDEL genomic prediction

## 26. Compare different marker sets

The same GS framework can be applied separately to:

```text
SNPs
INDELs
SNPs + INDELs
```

For example:

```text
Model 1: SNP-G
Model 2: INDEL-G
Model 3: SNP+INDEL-G
```

For each marker set:

1. create the genotype matrix;
2. construct the genomic relationship matrix;
3. use identical cross-validation folds;
4. estimate prediction ability;
5. compare mean prediction accuracy.

Using identical cross-validation partitions is important for fair comparisons among marker sets.

---

# Part XVI. Optional alternative GS model: ridge-regression BLUP

## 27. Marker-effect RR-BLUP

Instead of fitting GBLUP from a genomic relationship matrix, marker effects can be estimated directly using `mixed.solve()`.

```r
library(rrBLUP)

fit_rrblup <- mixed.solve(
    y = y,
    Z = geno
)

marker_effects <- fit_rrblup$u
```

Genomic breeding values can then be calculated as:

```r
gebv_rrblup <- geno %*% marker_effects
```

Because RR-BLUP and GBLUP are closely related under standard assumptions, GBLUP can be used as the primary baseline model and RR-BLUP as an alternative implementation.

---

# Part XVII. Recommended GS validation strategy

For manuscript-level analyses, prediction accuracy should ideally be evaluated using more than one validation strategy.

Recommended evaluations include:

### Random k-fold cross-validation

Randomly divide all individuals into training and validation sets.

This measures general prediction ability when training and validation individuals are genetically mixed.

### Population-stratified prediction

Train the model in one or several populations and predict another population.

This evaluates prediction transferability across population structure.

### Breeding-group prediction

Use elite breeding materials as the training population and predict candidate breeding materials.

This is often more relevant to practical selection.

### Environment-based prediction

For multi-environment phenotype datasets, train the model using phenotypes from selected environments and evaluate prediction in another environment.

This can be extended to genomic prediction models incorporating genotype-by-environment interactions.

---

# Part XVIII. Suggested output structure

A clean directory structure for the analysis is:

```text
GWAS_GS/
├── 01_VCF/
│   ├── merged.vcf.gz
│   ├── merged_snps.vcf.gz
│   └── merged_indels.vcf.gz
├── 02_PLINK/
│   ├── SNP.bed
│   ├── SNP.bim
│   ├── SNP.fam
│   ├── INDEL.bed
│   ├── INDEL.bim
│   └── INDEL.fam
├── 03_PCA/
├── 04_GRM/
├── 05_Phenotype/
├── 06_Heritability/
├── 07_GWAS/
└── 08_GS/
```

The GitHub repository does not need to contain large genotype files. Instead, it can contain:

- this workflow document;
- phenotype templates;
- sample lists;
- small example data;
- scripts used for data conversion and analysis.

---



# Part XX. Genomic Prediction Workflow

## 28. Overview

Genomic prediction extends genomic selection by using genome-wide marker information to predict phenotypic or breeding values for individuals whose phenotypes are unknown.

In this workflow, genomic prediction is treated as a complete prediction framework rather than only as a cross-validation step.

The general workflow is:

```text
Genotype data
     +
Phenotype data for training individuals
     ↓
Marker quality control
     ↓
Genotype matrix
     ↓
Training / validation / prediction populations
     ↓
Genomic prediction model
     ↓
Model validation
     ↓
Prediction of unphenotyped individuals
     ↓
Predicted phenotype or GEBV
     ↓
Ranking and candidate selection
```

The same framework can be applied to SNPs, INDELs, or combined marker sets.

---

## 29. Define training, validation, and prediction populations

Three groups can be distinguished:

```text
Training population
    Genotype + phenotype available

Validation population
    Genotype + phenotype available
    Phenotypes are hidden during model fitting and used only for evaluation

Prediction population
    Genotype available
    Phenotype unavailable
```

Example sample lists:

```text
training_samples.txt
validation_samples.txt
prediction_samples.txt
```

Each file contains one sample ID per line.

---

## 30. Prepare a genome-wide genotype matrix

Export additive marker genotypes from PLINK:

```bash
plink \
    --bfile SNP \
    --recode A \
    --out genomic_prediction
```

The resulting file is:

```text
genomic_prediction.raw
```

The genotype coding is generally:

```text
0 = homozygous reference
1 = heterozygous
2 = homozygous alternative
```

Read the genotype matrix in R:

```r
library(data.table)

geno_raw <- fread(
    "genomic_prediction.raw",
    data.table = FALSE
)

sample_id <- geno_raw$IID

geno <- geno_raw[, 7:ncol(geno_raw)]

rownames(geno) <- sample_id

geno <- as.matrix(geno)
```

---

## 31. Marker filtering before prediction

Remove markers with no variation:

```r
keep_var <- apply(
    geno,
    2,
    function(x) var(x, na.rm = TRUE) > 0
)

geno <- geno[, keep_var]
```

Remove markers with excessive missingness:

```r
marker_missing <- apply(
    geno,
    2,
    function(x) mean(is.na(x))
)

geno <- geno[, marker_missing <= 0.2]
```

Optional minor allele frequency filtering can be performed in PLINK before exporting the genotype matrix:

```bash
plink \
    --bfile SNP \
    --maf 0.01 \
    --geno 0.2 \
    --make-bed \
    --out SNP_GP_filtered
```

The exact filtering thresholds should be selected according to sample size and study design.

---

## 32. Genotype imputation

For a baseline analysis, missing genotypes can be imputed by marker mean:

```r
for (j in seq_len(ncol(geno))) {

    marker_mean <- mean(
        geno[, j],
        na.rm = TRUE
    )

    geno[is.na(geno[, j]), j] <- marker_mean
}
```

For large datasets or high missingness, a dedicated genotype-imputation method can be used before genomic prediction.

---

## 33. Prepare phenotype data

A phenotype file can be formatted as:

```text
Sample    Trait
sample1   10.2
sample2   11.5
sample3    9.8
```

Read phenotype data:

```r
pheno <- read.table(
    "phenotype.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
```

Match genotype and phenotype samples:

```r
common <- intersect(
    rownames(geno),
    pheno$Sample
)

geno_pheno <- geno[common, , drop = FALSE]

pheno <- pheno[
    match(common, pheno$Sample),
]
```

---

## 34. Construct the genomic relationship matrix

Using `rrBLUP`:

```r
library(rrBLUP)

G <- A.mat(
    geno_pheno,
    impute.method = "mean",
    return.imputed = FALSE
)
```

This genomic relationship matrix can be used by GBLUP for both model validation and prediction.

---

# Part XXI. Genomic Prediction Models

## 35. GBLUP

GBLUP predicts genomic breeding values using the genomic relationship matrix.

```r
y <- pheno$Trait
names(y) <- pheno$Sample

fit_gblup <- kin.blup(
    data = data.frame(
        ID = names(y),
        phenotype = y
    ),
    geno = "ID",
    pheno = "phenotype",
    K = G
)
```

Extract fitted genomic breeding values:

```r
gblup_gebv <- data.frame(
    Sample = names(fit_gblup$g),
    GEBV = fit_gblup$g
)
```

---

## 36. RR-BLUP

RR-BLUP estimates marker effects directly.

```r
fit_rrblup <- mixed.solve(
    y = y,
    Z = geno_pheno
)

marker_effects <- fit_rrblup$u
```

Calculate genomic values:

```r
pred_rrblup <- geno_pheno %*% marker_effects
```

---

## 37. Optional alternative prediction models

Depending on the genetic architecture of the trait, additional models can be compared.

Common genomic prediction models include:

```text
GBLUP
RR-BLUP
BayesA
BayesB
BayesC
Bayesian LASSO
Elastic Net
Random Forest
Gradient Boosting
Support Vector Regression
```

For a manuscript, GBLUP or RR-BLUP provides a useful baseline. Additional models can be added when prediction performance or genetic architecture is a major focus.

---

# Part XXII. Cross-Validation for Genomic Prediction

## 38. Repeated five-fold cross-validation

Repeated cross-validation provides a robust estimate of prediction performance.

```r
library(rrBLUP)

set.seed(2026)

n_rep <- 100
k <- 5

cv_results <- data.frame(
    Replicate = integer(),
    Correlation = numeric(),
    RMSE = numeric(),
    MAE = numeric()
)

for (rep in 1:n_rep) {

    fold <- sample(
        rep(
            1:k,
            length.out = length(y)
        )
    )

    prediction <- rep(
        NA_real_,
        length(y)
    )

    for (f in 1:k) {

        test_index <- which(
            fold == f
        )

        y_train <- y
        y_train[test_index] <- NA

        fit <- kin.blup(
            data = data.frame(
                ID = names(y_train),
                phenotype = y_train
            ),
            geno = "ID",
            pheno = "phenotype",
            K = G
        )

        prediction[test_index] <- fit$g[
            match(
                names(y)[test_index],
                names(fit$g)
            )
        ]
    }

    cor_value <- cor(
        prediction,
        y,
        use = "complete.obs"
    )

    rmse_value <- sqrt(
        mean(
            (prediction - y)^2,
            na.rm = TRUE
        )
    )

    mae_value <- mean(
        abs(prediction - y),
        na.rm = TRUE
    )

    cv_results <- rbind(
        cv_results,
        data.frame(
            Replicate = rep,
            Correlation = cor_value,
            RMSE = rmse_value,
            MAE = mae_value
        )
    )
}
```

Summarize predictive performance:

```r
summary_results <- data.frame(
    Mean_Correlation = mean(cv_results$Correlation),
    SD_Correlation = sd(cv_results$Correlation),
    Mean_RMSE = mean(cv_results$RMSE),
    SD_RMSE = sd(cv_results$RMSE),
    Mean_MAE = mean(cv_results$MAE),
    SD_MAE = sd(cv_results$MAE)
)

write.table(
    cv_results,
    "genomic_prediction_cross_validation.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

write.table(
    summary_results,
    "genomic_prediction_summary.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

Recommended statistics include:

```text
Pearson correlation
RMSE
MAE
```

The correlation between observed and predicted values is commonly referred to as **prediction ability**.

---

## 39. Prediction accuracy adjusted by heritability

When broad- or narrow-sense heritability is available, prediction accuracy can optionally be reported as:

```text
Prediction accuracy = correlation(observed, predicted) / sqrt(h²)
```

In R:

```r
h2 <- 0.50

prediction_accuracy <- mean(
    cv_results$Correlation
) / sqrt(h2)

prediction_accuracy
```

This adjustment should only be used when the heritability definition is compatible with the genomic prediction model and phenotype used.

---

# Part XXIII. Independent Validation

## 40. Train in one population and predict another

Random cross-validation may overestimate prediction performance when closely related individuals occur in both training and validation sets.

A stronger validation strategy is:

```text
Population A + Population B
          ↓
       Training
          ↓
     Prediction model
          ↓
     Population C
          ↓
   Independent validation
```

Example:

```r
train_id <- scan(
    "training_samples.txt",
    what = character()
)

test_id <- scan(
    "validation_samples.txt",
    what = character()
)

all_id <- c(
    train_id,
    test_id
)

G_sub <- G[
    all_id,
    all_id
]

pheno_model <- data.frame(
    ID = all_id,
    phenotype = NA_real_
)

pheno_model$phenotype[
    match(
        train_id,
        pheno_model$ID
    )
] <- y[
    match(
        train_id,
        names(y)
    )
]

fit <- kin.blup(
    data = pheno_model,
    geno = "ID",
    pheno = "phenotype",
    K = G_sub
)

pred_test <- fit$g[
    match(
        test_id,
        names(fit$g)
    )
]
```

Evaluate:

```r
obs_test <- y[
    match(
        test_id,
        names(y)
    )
]

cor(
    pred_test,
    obs_test,
    use = "complete.obs"
)
```

---

# Part XXIV. Predict Unphenotyped Individuals

## 41. Construct a combined relationship matrix

To predict individuals without phenotype data, the genomic relationship matrix must include both training and prediction individuals.

```r
training_id <- scan(
    "training_samples.txt",
    what = character()
)

prediction_id <- scan(
    "prediction_samples.txt",
    what = character()
)

all_id <- c(
    training_id,
    prediction_id
)

geno_all <- geno[
    all_id,
    ,
    drop = FALSE
]

G_all <- A.mat(
    geno_all,
    impute.method = "mean",
    return.imputed = FALSE
)
```

---

## 42. Fit the prediction model

Create a phenotype vector in which only training individuals have observed values:

```r
prediction_data <- data.frame(
    ID = all_id,
    phenotype = NA_real_
)

prediction_data$phenotype[
    match(
        training_id,
        prediction_data$ID
    )
] <- pheno$Trait[
    match(
        training_id,
        pheno$Sample
    )
]
```

Fit GBLUP:

```r
fit_prediction <- kin.blup(
    data = prediction_data,
    geno = "ID",
    pheno = "phenotype",
    K = G_all
)
```

---

## 43. Export genomic predictions

Extract predicted GEBVs for the prediction population:

```r
predicted_gebv <- data.frame(
    Sample = prediction_id,
    GEBV = fit_prediction$g[
        match(
            prediction_id,
            names(fit_prediction$g)
        )
    ]
)
```

Rank candidate individuals:

```r
predicted_gebv <- predicted_gebv[
    order(
        predicted_gebv$GEBV,
        decreasing = TRUE
    ),
]
```

Export:

```r
write.table(
    predicted_gebv,
    "prediction_population_GEBV.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

For traits in which lower values are favorable, reverse the ranking direction.

---

# Part XXV. Predict Phenotypic Values

## 44. From GEBV to predicted phenotype

For selection based on breeding value, GEBV is usually the principal output.

If predicted phenotypic values are required, they can be reconstructed as:

```text
Predicted phenotype = population mean + GEBV
```

Example:

```r
phenotype_mean <- mean(
    y,
    na.rm = TRUE
)

predicted_phenotype <- data.frame(
    Sample = predicted_gebv$Sample,
    Predicted_Phenotype =
        phenotype_mean + predicted_gebv$GEBV
)
```

This formulation assumes that no additional fixed effects need to be added back to the prediction.

---

# Part XXVI. Multi-Trait Genomic Prediction

## 45. Predict several traits independently

For multiple phenotypes, the basic model can be fitted trait by trait:

```r
trait_columns <- c(
    "Trait1",
    "Trait2",
    "Trait3"
)

for (trait in trait_columns) {

    y_trait <- pheno[[trait]]

    names(y_trait) <- pheno$Sample

    fit <- kin.blup(
        data = data.frame(
            ID = names(y_trait),
            phenotype = y_trait
        ),
        geno = "ID",
        pheno = "phenotype",
        K = G
    )

    output <- data.frame(
        Sample = names(fit$g),
        GEBV = fit$g
    )

    write.table(
        output,
        paste0(
            trait,
            "_GEBV.txt"
        ),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}
```

A dedicated multi-trait genomic prediction model can be considered when genetic correlations among traits are of interest.

---

# Part XXVII. Multi-Environment Genomic Prediction

## 46. Genotype-by-environment prediction

For multi-location or multi-year phenotypes, genomic prediction can be extended to model:

```text
Genotype effect
+
Environment effect
+
Genotype × Environment interaction
```

The conceptual model is:

```text
y = μ + E + G + G×E + e
```

Possible validation schemes include:

```text
CV1:
Predict newly genotyped individuals that have never been phenotyped.

CV2:
Predict missing genotype-by-environment combinations for individuals
that have phenotypes in some environments.

Leave-one-environment-out:
Train on all environments except one and predict the omitted environment.
```

This is particularly useful when phenotype data are collected across multiple years, locations, or treatments.

---

# Part XXVIII. Marker-Set Comparison for Genomic Prediction

## 47. Compare SNP, INDEL, and combined marker sets

Prepare three marker datasets:

```text
SNP
INDEL
SNP + INDEL
```

For each dataset:

```text
1. Apply identical sample filtering.
2. Apply comparable marker quality control.
3. Construct a genomic relationship matrix.
4. Use the same cross-validation folds.
5. Fit the same prediction model.
6. Calculate correlation, RMSE, and MAE.
7. Compare predictive performance.
```

Example summary table:

```text
Marker_set     Correlation     RMSE     MAE
SNP
INDEL
SNP_INDEL
```

Using identical validation partitions is essential for a fair comparison.

---

# Part XXIX. Genomic Prediction Outputs

## 48. Recommended output files

A genomic prediction analysis can generate:

```text
genomic_prediction_cross_validation.txt
genomic_prediction_summary.txt
prediction_population_GEBV.txt
predicted_phenotype.txt
SNP_prediction_summary.txt
INDEL_prediction_summary.txt
SNP_INDEL_prediction_summary.txt
```

For each trait, report:

```text
Number of training individuals
Number of validation individuals
Number of markers
Prediction model
Cross-validation design
Number of replicates
Mean prediction ability
Standard deviation
RMSE
MAE
Heritability
Prediction accuracy, if calculated
```

---

# Part XXX. Integrated GWAS, Genomic Selection, and Genomic Prediction Workflow

The complete integrated analysis is:

```text
                         Genotype data
                              │
                 ┌────────────┴────────────┐
                 │                         │
               GWAS                 Genomic prediction
                 │                         │
          PCA + GRM + MLMA           Marker QC + GRM
                 │                         │
       Significant loci             Training population
                 │                         │
      Candidate genes          Cross-validation / validation
                 │                         │
                 │                 Prediction model
                 │                         │
                 │                Predicted GEBV
                 │                         │
                 └────────────┬────────────┘
                              │
                        Candidate selection
                              │
                       Breeding decisions
```

GWAS, genomic selection, and genomic prediction provide complementary information:

- **GWAS** identifies loci statistically associated with phenotypic variation.
- **Genomic selection** uses genome-wide marker information to rank breeding candidates based on genomic breeding values.
- **Genomic prediction** evaluates and applies predictive models to estimate phenotypes or breeding values for individuals without phenotype records.

For breeding applications, prediction models should always be evaluated using an explicit validation design before being used to rank unphenotyped candidate individuals.


# Part XIX. Summary

The complete GWAS workflow is:

```text
VCF
 ↓
SNP / INDEL separation
 ↓
PLINK
 ↓
PCA
 ↓
GRM
 ↓
Phenotype
 ↓
REML heritability
 ↓
GCTA-MLMA GWAS
 ↓
Significant loci
```

The complete GS workflow is:

```text
Genotype
   +
Phenotype
   ↓
Marker QC
   ↓
Genomic relationship matrix
   ↓
Training / validation partition
   ↓
GBLUP or RR-BLUP
   ↓
Cross-validation
   ↓
Prediction accuracy
   ↓
Final GEBV
   ↓
Selection of candidate individuals
```
