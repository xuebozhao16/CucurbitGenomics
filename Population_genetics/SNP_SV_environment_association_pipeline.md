# Genotype–Environment Association Analysis for SNPs and SVs

This workflow identifies **SNPs and structural variants (SVs) associated with environmental variation** using three complementary methods: **LFMM**, **BayPass**, and **Bayenv2/pybayenv2**. The same environmental variables can be tested across SNPs and SVs, but genotype preparation differs slightly because SNPs are usually represented as nucleotide substitutions whereas SVs may be deletions, insertions, inversions, duplications, or other biallelic presence/absence variants.

For all analyses, use the same sample set, population assignments, environmental variables, and variant identifiers whenever possible so that results from different methods can be compared directly.

## Data preparation

Start from filtered SNP and SV VCF files:

```text
SNP.filtered.vcf.gz
SV.filtered.vcf.gz
```

The recommended minimum filtering steps are:

```text
biallelic variants
reliable genotype quality
reasonable missing-data rate
minimum minor allele count / minor allele frequency
removal of monomorphic variants
consistent sample IDs
consistent chromosome names
```

For SVs, retain only records with interpretable genotypes. A simple presence/absence coding is appropriate for biallelic SVs:

```text
0/0 -> 0 copies of the alternative SV allele
0/1 -> 1 copy
1/1 -> 2 copies
```

Multiallelic SVs should either be decomposed into biologically interpretable biallelic records or analyzed separately.

Index the VCF files:

```bash
bcftools index -f SNP.filtered.vcf.gz
bcftools index -f SV.filtered.vcf.gz
```

Check sample names:

```bash
bcftools query -l SNP.filtered.vcf.gz > samples.txt
bcftools query -l SV.filtered.vcf.gz > SV.samples.txt
```

Confirm that SNP and SV datasets contain the same individuals:

```bash
diff \
    <(sort samples.txt) \
    <(sort SV.samples.txt)
```

Create stable variant IDs. For SNPs:

```bash
bcftools annotate \
    --set-id '%CHROM_%POS_%REF_%FIRST_ALT' \
    SNP.filtered.vcf.gz \
    -Oz \
    -o SNP.withID.vcf.gz
```

For SVs, preserve an existing unique ID whenever possible. If the VCF contains `SVTYPE` and `END`, a useful external metadata table is:

```bash
bcftools query \
    -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%REF\t%ALT\n' \
    SV.filtered.vcf.gz \
    > SV.variant_info.tsv
```

Create a sample metadata file:

```text
Sample      Population
sample01    Pop1
sample02    Pop1
sample03    Pop2
sample04    Pop2
sample05    Pop3
sample06    Pop3
```

Environmental data should be matched to either individuals or populations depending on the method.

For individual-based LFMM:

```text
Sample      bio_3   bio_5   bio_8   bio_9   bio_12  bio_19
sample01    ...
sample02    ...
...
```

For population-based BayPass/Bayenv2:

```text
Population  bio_3   bio_5   bio_8   bio_9   bio_12  bio_19
Pop1        ...
Pop2        ...
Pop3        ...
```

Before association testing, remove strongly correlated environmental variables. A simple R workflow is:

```r
env <- read.table(
    "environment.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

env_numeric <- env[, sapply(env, is.numeric)]

cor_matrix <- cor(
    env_numeric,
    use = "pairwise.complete.obs",
    method = "pearson"
)

write.table(
    cor_matrix,
    "environment_correlation.tsv",
    sep = "\t",
    quote = FALSE
)
```

A practical approach is to retain one variable from highly correlated pairs, for example:

```text
|r| > 0.7
```

or to use environmental PCA and analyze a smaller set of orthogonal environmental PCs.

Standardize environmental variables before LFMM, BayPass, or Bayenv2:

```r
env_scaled <- scale(env_numeric)

write.table(
    env_scaled,
    "environment_scaled.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

To construct an individual-level genotype dosage matrix for SNPs or biallelic SVs, one convenient R approach is:

```r
library(vcfR)

vcf <- read.vcfR(
    "SNP.withID.vcf.gz",
    verbose = FALSE
)

gt <- extract.gt(
    vcf,
    element = "GT",
    as.numeric = FALSE
)

gt_to_dosage <- function(x) {

    x[x %in% c("./.", ".|.", ".")] <- NA
    x <- gsub("\\|", "/", x)

    out <- rep(NA_real_, length(x))

    out[x == "0/0"] <- 0
    out[x %in% c("0/1", "1/0")] <- 1
    out[x == "1/1"] <- 2

    out
}

geno <- apply(
    gt,
    2,
    gt_to_dosage
)

geno <- t(geno)
```

The final LFMM genotype matrix must be:

```text
rows    = individuals
columns = variants
values  = 0, 1, 2
```

The same code can be applied to a biallelic SV VCF:

```r
sv_vcf <- read.vcfR(
    "SV.filtered.vcf.gz",
    verbose = FALSE
)

sv_gt <- extract.gt(
    sv_vcf,
    element = "GT",
    as.numeric = FALSE
)

sv_geno <- apply(
    sv_gt,
    2,
    gt_to_dosage
)

sv_geno <- t(sv_geno)
```

Remove variants with no variation or excessive missingness:

```r
filter_genotype_matrix <- function(Y, max_missing = 0.2) {

    missing_rate <- apply(
        Y,
        2,
        function(x) mean(is.na(x))
    )

    variable <- apply(
        Y,
        2,
        function(x) {
            length(unique(x[!is.na(x)])) > 1
        }
    )

    Y[
        ,
        missing_rate <= max_missing & variable,
        drop = FALSE
    ]
}

geno <- filter_genotype_matrix(geno)
sv_geno <- filter_genotype_matrix(sv_geno)
```

Variant metadata must remain in the same column order as the genotype matrix.

## LFMM

LFMM tests associations between genotypes and environmental variables while accounting for unobserved population structure using latent factors.

The LFMM model can be written as:

$$
Y = X B^T + U V^T + E
$$

where \(Y\) is the individual-by-variant genotype matrix, \(X\) is the individual-by-environment matrix, \(B\) contains genotype–environment effect sizes, \(U\) contains individual latent-factor scores, \(V\) contains locus loadings, and \(E\) is residual variation.

For locus \(j\) and individual \(i\):

$$
Y_{ij}
=
X_i \beta_j
+
\sum_{k=1}^{K} U_{ik}V_{jk}
+
\varepsilon_{ij}
$$

Load the required packages:

```r
library(LEA)
library(lfmm)
library(qvalue)
library(ggplot2)
library(dplyr)
```

Convert the filtered VCF to LEA/LFMM format:

```r
vcf2lfmm(
    input.file = "SNP.withID.vcf.gz"
)
```

The same procedure can be used for the SV VCF if the records are biallelic and genotype-coded:

```r
vcf2lfmm(
    input.file = "SV.filtered.vcf.gz"
)
```

Use `sNMF` to explore the number of ancestry components:

```r
project <- snmf(
    "SNP.withID.lfmm",
    K = 1:10,
    entropy = TRUE,
    repetitions = 20,
    project = "new"
)

plot(
    project,
    pch = 19,
    cex = 1.1
)
```

Inspect cross-entropy across \(K\):

```r
K_values <- 1:10

ce <- sapply(
    K_values,
    function(k) {
        min(
            cross.entropy(
                project,
                K = k
            )
        )
    }
)

plot(
    K_values,
    ce,
    type = "b",
    xlab = "K",
    ylab = "Cross-entropy"
)
```

Choose \(K\) based on the minimum or plateau of cross-entropy together with PCA/ADMIXTURE/sNMF biological interpretation.

Prepare the environmental matrix in the exact same individual order as the genotype matrix:

```r
env <- read.table(
    "environment_individual.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

env_subset <- env[
    ,
    c(
        "bio_3",
        "bio_5",
        "bio_8",
        "bio_9",
        "bio_12",
        "bio_19"
    )
]

X <- scale(
    as.matrix(env_subset)
)
```

Use the SNP genotype matrix:

```r
Y <- geno
```

or the SV genotype matrix:

```r
Y <- sv_geno
```

Fit ridge LFMM:

```r
K <- 3

lfmm_fit <- lfmm_ridge(
    Y = Y,
    X = X,
    K = K
)
```

The ridge model minimizes a regularized least-squares criterion:

$$
\min_{B,U,V}
\left\|
Y-XB^T-UV^T
\right\|_F^2
+
\lambda
\left\|
B
\right\|_F^2
$$

Perform association testing:

```r
lfmm_test_result <- lfmm_test(
    Y = Y,
    X = X,
    lfmm = lfmm_fit,
    calibrate = "gif"
)
```

The genomic inflation factor is:

$$
\lambda_{GIF}
=
\frac{
\mathrm{median}(Z_j^2)
}{
\mathrm{median}(\chi_1^2)
}
$$

Inspect:

```r
lfmm_test_result$gif
```

Extract raw and calibrated P-values:

```r
raw_p <- as.data.frame(
    lfmm_test_result$pvalue
)

calibrated_p <- as.data.frame(
    lfmm_test_result$calibrated.pvalue
)

z_scores <- as.data.frame(
    lfmm_test_result$score
)

colnames(raw_p) <- colnames(X)
colnames(calibrated_p) <- colnames(X)
colnames(z_scores) <- colnames(X)
```

Inspect P-value distributions:

```r
par(
    mfrow = c(
        ncol(X),
        2
    )
)

for (i in seq_len(ncol(X))) {

    hist(
        raw_p[[i]],
        breaks = 50,
        main = paste("Raw:", colnames(X)[i]),
        xlab = "P-value"
    )

    hist(
        calibrated_p[[i]],
        breaks = 50,
        main = paste("GIF calibrated:", colnames(X)[i]),
        xlab = "P-value"
    )
}
```

Generate QQ plots:

```r
par(
    mfrow = c(
        2,
        3
    )
)

for (i in seq_len(ncol(X))) {

    observed <- sort(
        calibrated_p[[i]]
    )

    expected <- ppoints(
        length(observed)
    )

    plot(
        -log10(expected),
        -log10(observed),
        pch = 19,
        cex = 0.5,
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        main = colnames(X)[i]
    )

    abline(
        0,
        1
    )
}
```

Calculate q-values separately for each environmental variable:

```r
qvalue_matrix <- matrix(
    NA_real_,
    nrow = nrow(calibrated_p),
    ncol = ncol(calibrated_p)
)

colnames(qvalue_matrix) <- colnames(calibrated_p)

for (i in seq_len(ncol(calibrated_p))) {

    qvalue_matrix[, i] <- qvalue(
        calibrated_p[[i]]
    )$qvalues
}
```

Define candidates using:

$$
q < 0.05
$$

```r
LFMM_candidates <- lapply(
    seq_len(ncol(qvalue_matrix)),
    function(i) {
        which(
            qvalue_matrix[, i] < 0.05
        )
    }
)

names(LFMM_candidates) <- colnames(qvalue_matrix)
```

If `variant_info` contains the variant IDs:

```r
candidate_tables <- lapply(
    names(LFMM_candidates),
    function(env_name) {

        idx <- LFMM_candidates[[env_name]]

        data.frame(
            variant_info[idx, ],
            Environment = env_name,
            Z = z_scores[idx, env_name],
            P = calibrated_p[idx, env_name],
            Q = qvalue_matrix[idx, env_name]
        )
    }
)

LFMM_result <- bind_rows(
    candidate_tables
)
```

Save:

```r
write.table(
    LFMM_result,
    "LFMM_environment_associated_variants.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

The sign of the LFMM Z-score or effect estimate can be retained as the association direction:

```text
positive effect  -> alternative allele increases with the environmental variable
negative effect  -> alternative allele decreases with the environmental variable
```

Run the same LFMM pipeline separately for SNPs and SVs.

## BayPass

BayPass tests population-level allele-frequency differentiation and genotype–environment association while explicitly modeling the covariance structure among populations.

For each locus \(j\), let the population allele-frequency vector be:

$$
\boldsymbol{p}_j
=
(p_{1j},p_{2j},\ldots,p_{Kj})^T
$$

BayPass estimates a population covariance matrix:

$$
\Omega
$$

which describes shared demographic history and covariance of allele frequencies among populations.

The environmental association model can be written conceptually as:

$$
\boldsymbol{p}_j
=
\mu_j\mathbf{1}
+
\beta_j\boldsymbol{Z}
+
\boldsymbol{\varepsilon}_j
$$

with:

$$
\boldsymbol{\varepsilon}_j
\sim
N(
0,
\sigma_j^2\Omega
)
$$

BayPass uses **population allele counts**, not individual genotype rows.

Starting from a dosage matrix `Y`, convert genotypes to reference and alternative allele counts:

```r
sample_meta <- read.table(
    "sample_population.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

stopifnot(
    all(
        rownames(Y) ==
        sample_meta$Sample
    )
)

populations <- unique(
    sample_meta$Population
)
```

Construct BayPass allele counts:

```r
make_baypass_counts <- function(
    Y,
    populations,
    pop_vector
) {

    result <- matrix(
        NA_integer_,
        nrow = ncol(Y),
        ncol = length(populations) * 2
    )

    for (p in seq_along(populations)) {

        ids <- which(
            pop_vector == populations[p]
        )

        sub <- Y[
            ids,
            ,
            drop = FALSE
        ]

        alt_count <- apply(
            sub,
            2,
            function(x) {
                sum(
                    x,
                    na.rm = TRUE
                )
            }
        )

        n_called <- apply(
            sub,
            2,
            function(x) {
                sum(
                    !is.na(x)
                )
            }
        )

        ref_count <- (
            2 * n_called
            - alt_count
        )

        result[, 2 * p - 1] <- ref_count
        result[, 2 * p] <- alt_count
    }

    result
}

baypass_counts <- make_baypass_counts(
    Y = geno,
    populations = populations,
    pop_vector = sample_meta$Population
)
```

For SV analysis, replace `geno` with `sv_geno`.

Write BayPass genotype files without headers:

```r
write.table(
    baypass_counts,
    "SNP.baypass.geno",
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
```

For SVs:

```r
write.table(
    sv_baypass_counts,
    "SV.baypass.geno",
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
```

Environmental variables must be in the same population order as the allele-count columns.

```r
pop_env <- read.table(
    "population_environment.tsv",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

pop_env <- pop_env[
    match(
        populations,
        pop_env$Population
    ),
]

env_matrix <- scale(
    as.matrix(
        pop_env[
            ,
            c(
                "bio_3",
                "bio_5",
                "bio_8",
                "bio_9",
                "bio_12",
                "bio_19"
            )
        ]
    )
)
```

Write one environmental variable at a time:

```r
write.table(
    env_matrix[, "bio_3", drop = FALSE],
    "bio_3.env",
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)
```

Run the BayPass core model first:

```bash
g_baypass \
    -npop NPOP \
    -gfile SNP.baypass.geno \
    -coremodel \
    -outprefix SNP_core
```

Repeat independent core-model runs to verify covariance-matrix stability:

```bash
for run in 1 2 3; do

    g_baypass \
        -npop NPOP \
        -gfile SNP.baypass.geno \
        -coremodel \
        -outprefix SNP_core_run${run}

done
```

The core model estimates the covariance matrix \(\Omega\) and provides the \(X^TX\) differentiation statistic.

To test environmental association, run the auxiliary covariate model:

```bash
g_baypass \
    -npop NPOP \
    -gfile SNP.baypass.geno \
    -efile bio_3.env \
    -scalecov \
    -auxmodel \
    -outprefix SNP_bio3
```

Run each environmental variable separately:

```bash
for env in bio_3 bio_5 bio_8 bio_9 bio_12 bio_19; do

    g_baypass \
        -npop NPOP \
        -gfile SNP.baypass.geno \
        -efile ${env}.env \
        -scalecov \
        -auxmodel \
        -outprefix SNP_${env}

done
```

Apply the same workflow to SV allele counts:

```bash
for env in bio_3 bio_5 bio_8 bio_9 bio_12 bio_19; do

    g_baypass \
        -npop NPOP \
        -gfile SV.baypass.geno \
        -efile ${env}.env \
        -scalecov \
        -auxmodel \
        -outprefix SV_${env}

done
```

BayPass environmental association is summarized with a Bayes factor:

$$
BF_j
=
\frac{
P(D_j \mid M_1)
}{
P(D_j \mid M_0)
}
$$

where \(M_1\) includes the environmental effect and \(M_0\) does not.

Retain at least:

```text
variant ID
environmental variable
Bayes factor
posterior inclusion probability, if reported
effect direction / regression coefficient, if reported
XtX
```

For publication-level analyses, pseudo-observed datasets generated under the inferred covariance model can be used to establish empirical thresholds for \(X^TX\) and related statistics.

## Bayenv2 / pybayenv2

Bayenv2 also models population allele-frequency covariance and tests association between allele frequencies and environmental variables.

A simplified allele-frequency model is:

$$
\boldsymbol{p}_j
\sim
N(
\mu_j\mathbf{1},
\mu_j(1-\mu_j)\Omega
)
$$

where \(\Omega\) is estimated from a large set of approximately neutral loci.

For environmental association:

$$
\boldsymbol{p}_j
=
\mu_j\mathbf{1}
+
\beta_j\boldsymbol{E}
+
\boldsymbol{\varepsilon}_j
$$

with:

$$
\boldsymbol{\varepsilon}_j
\sim
N(
0,
\Omega
)
$$

A direct Bayenv2 workflow is:

```text
1. Select a large set of approximately neutral markers.
2. Estimate the covariance matrix Omega.
3. Run multiple independent covariance chains.
4. Compare or average covariance estimates.
5. Test each environmental variable.
6. Repeat association runs.
7. Summarize Bayes factors and rank correlations.
```

If using the `pybayenv2` wrapper, the input format is GENEPOP.

Load the wrapper configuration:

```bash
source pybayenv.conf
```

Display options:

```bash
pybayenv -h
```

A general run has the form:

```bash
pybayenv \
    -f genotype.gp \
    -e environment.env \
    -n 3 \
    -z 1 \
    -p 8 \
    -i 10000
```

Prepare separate genotype inputs for SNPs and SVs:

```text
SNP.genepop
SV.genepop
```

For biallelic SVs, each locus is treated as a diploid marker with reference and alternative SV alleles.

Environmental data should contain population-level values in the same population order as the genotype file.

Retain:

```text
Bayes factor
Spearman rho
Spearman P-value
run-to-run consistency
```

The Spearman rank correlation is:

$$
\rho
=
1
-
\frac{
6\sum_i d_i^2
}{
n(n^2-1)
}
$$

when there are no tied ranks.

An example candidate rule used in some workflows is:

```text
BF >= 10
and
Spearman P < 0.05
```

but this should be treated as a study-specific rule rather than a universal threshold.

Because Bayenv2 MCMC results can vary among runs, repeat the analysis:

```bash
for run in 1 2 3 4 5; do

    pybayenv \
        -f SNP.genepop \
        -e environment.env \
        -n 3 \
        -z 1 \
        -p 8 \
        -i 10000

done
```

Summarize independent runs in R:

```r
results <- read.table(
    "Bayenv_all_runs.tsv",
    header = TRUE,
    sep = "\t"
)

library(dplyr)

summary_result <- results %>%
    group_by(
        Variant,
        Environment
    ) %>%
    summarise(
        Median_BF = median(
            BF,
            na.rm = TRUE
        ),
        Median_rho = median(
            Spearman_rho,
            na.rm = TRUE
        ),
        Median_P = median(
            Spearman_P,
            na.rm = TRUE
        ),
        .groups = "drop"
    )
```

Extract candidates:

```r
Bayenv_candidates <- summary_result %>%
    filter(
        Median_BF >= 10,
        Median_P < 0.05
    )

write.table(
    Bayenv_candidates,
    "Bayenv_environment_associated_variants.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

Run SNPs and SVs independently but use the same environmental table.

## Integrating SNP and SV results

After LFMM, BayPass, and Bayenv2 are complete, combine all candidate statistics into one variant-level table:

```text
Variant
Variant_type
CHR
POS
END
SVTYPE
Environment
LFMM_Z
LFMM_P
LFMM_Q
BayPass_BF
BayPass_XtX
Bayenv_BF
Bayenv_rho
Bayenv_P
```

For SNPs:

```text
Variant_type = SNP
```

For structural variants:

```text
Variant_type = DEL
Variant_type = INS
Variant_type = INV
Variant_type = DUP
...
```

Merge results:

```r
lfmm <- read.table(
    "LFMM_environment_associated_variants.tsv",
    header = TRUE,
    sep = "\t"
)

baypass <- read.table(
    "BayPass_environment_associated_variants.tsv",
    header = TRUE,
    sep = "\t"
)

bayenv <- read.table(
    "Bayenv_environment_associated_variants.tsv",
    header = TRUE,
    sep = "\t"
)

combined <- full_join(
    lfmm,
    baypass,
    by = c(
        "Variant",
        "Environment"
    )
)

combined <- full_join(
    combined,
    bayenv,
    by = c(
        "Variant",
        "Environment"
    )
)
```

Add cross-method support:

```r
combined$LFMM_support <- (
    !is.na(combined$LFMM_Q) &
    combined$LFMM_Q < 0.05
)

combined$BayPass_support <- (
    !is.na(combined$BayPass_BF)
)

combined$Bayenv_support <- (
    !is.na(combined$Bayenv_BF) &
    combined$Bayenv_BF >= 10 &
    combined$Bayenv_P < 0.05
)

combined$Methods_supported <- (
    combined$LFMM_support +
    combined$BayPass_support +
    combined$Bayenv_support
)
```

A stringent cross-method candidate set can be defined as variants supported by at least two methods:

```r
robust_candidates <- combined[
    combined$Methods_supported >= 2,
]
```

Save:

```r
write.table(
    combined,
    "GEA_all_SNP_SV_results.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

write.table(
    robust_candidates,
    "GEA_cross_method_candidates.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

For downstream annotation, retain association direction as well as significance:

```text
LFMM Z-score
BayPass regression effect
Bayenv Spearman rho
```

For SNPs, use the SNP coordinate as a one-base interval.

For SVs, use the full SV interval where appropriate:

```text
DEL / INV / DUP -> POS to END
```

For insertions, use a consistent breakpoint interval or the representation used by the SV callset.

Intersect candidate variants with gene annotations:

```bash
bedtools intersect \
    -a GEA_candidates.bed \
    -b genes.gff3 \
    -wa \
    -wb \
    > GEA_candidate_gene_overlap.tsv
```

Identify nearby genes:

```bash
bedtools window \
    -a GEA_candidates.bed \
    -b genes.gff3 \
    -w 5000 \
    > GEA_candidate_nearby_genes.tsv
```