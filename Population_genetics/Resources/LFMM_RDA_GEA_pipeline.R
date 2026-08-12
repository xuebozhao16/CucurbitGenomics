#!/usr/bin/env Rscript

# ==============================================================================
# Genotype-Environment Association Analysis with LFMM and RDA
# ==============================================================================
# This script is adapted from the Landscape Genetic Data Analysis with R
# worked example and reorganized as a reusable analysis pipeline.
#
# INPUT:
#   1. Genotype matrix:
#        - rows = individuals
#        - columns = loci / variants
#        - values = 0, 1, 2
#        - first column = sample ID
#
#   2. Environmental table:
#        - rows = individuals
#        - first column = sample ID
#        - remaining selected columns = environmental predictors
#
# OUTPUT:
#   - environmental correlation matrix
#   - environmental PCA results
#   - genetic PCA results
#   - LFMM statistics and candidates
#   - RDA model statistics and candidates
#   - overlap between LFMM and RDA candidates
#
# NOTE:
#   LFMM and RDA require complete genotype matrices. Missing genotypes are
#   imputed here using the most common genotype at each locus, following the
#   workflow in the source tutorial.
# ==============================================================================


# ==============================================================================
# 1. User settings
# ==============================================================================

GENOTYPE_FILE <- "genotype_matrix.tsv"
ENVIRONMENT_FILE <- "environment.tsv"
OUTPUT_DIR <- "GEA_results"

# Column containing sample IDs in both files.
SAMPLE_COLUMN <- "Sample"

# Environmental variables to analyze.
# Replace these names with columns in ENVIRONMENT_FILE.
ENV_VARS <- c(
  "bio_3",
  "bio_5",
  "bio_8",
  "bio_9",
  "bio_12",
  "bio_19"
)

# Correlation threshold used to flag highly correlated predictors.
COR_THRESHOLD <- 0.70

# LFMM settings.
LFMM_K <- 3
LFMM_FDR <- 0.05

# If TRUE, LFMM uses environmental PC1 as one synthetic predictor.
# If FALSE, LFMM analyzes all ENV_VARS simultaneously.
USE_ENV_PC1_FOR_LFMM <- FALSE

# RDA settings.
RDA_N_AXES <- 3
RDA_Z_CUTOFF <- 3

# If TRUE, run partial RDA using selected genetic PCs as conditioning variables.
RUN_PARTIAL_RDA <- FALSE

# Genetic PCs used for partial RDA.
# Example: c(1, 2)
PARTIAL_RDA_PC_AXES <- c(1)

# Number of permutations for RDA significance tests.
RDA_PERMUTATIONS <- 999

# Random seed.
SEED <- 2026


# ==============================================================================
# 2. Required packages
# ==============================================================================

required_packages <- c(
  "vegan",
  "lfmm",
  "qvalue"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
}

suppressPackageStartupMessages({
  library(vegan)
  library(lfmm)
  library(qvalue)
})

set.seed(SEED)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Helper functions
# ==============================================================================

impute_mode <- function(x) {

  if (!anyNA(x)) {
    return(x)
  }

  observed <- x[!is.na(x)]

  if (length(observed) == 0) {
    return(x)
  }

  mode_value <- as.numeric(
    names(
      which.max(
        table(observed)
      )
    )
  )

  x[is.na(x)] <- mode_value

  x
}


remove_nonvariable_loci <- function(genotype_matrix) {

  keep <- apply(
    genotype_matrix,
    2,
    function(x) {
      length(unique(x[!is.na(x)])) > 1
    }
  )

  genotype_matrix[
    ,
    keep,
    drop = FALSE
  ]
}


loading_outliers <- function(x, z = 3) {

  limits <- mean(
    x,
    na.rm = TRUE
  ) +
    c(-1, 1) *
    z *
    sd(
      x,
      na.rm = TRUE
    )

  which(
    x < limits[1] |
      x > limits[2]
  )
}


write_tsv <- function(x, filename) {

  write.table(
    x,
    file = file.path(
      OUTPUT_DIR,
      filename
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


save_pdf <- function(filename, expr, width = 8, height = 6) {

  pdf(
    file.path(
      OUTPUT_DIR,
      filename
    ),
    width = width,
    height = height
  )

  on.exit(
    dev.off(),
    add = TRUE
  )

  force(expr)
}


# ==============================================================================
# 4. Import genotype data
# ==============================================================================

geno_df <- read.table(
  GENOTYPE_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

if (!SAMPLE_COLUMN %in% colnames(geno_df)) {
  stop(
    paste0(
      "Sample column '",
      SAMPLE_COLUMN,
      "' was not found in genotype file."
    )
  )
}

sample_ids <- geno_df[[SAMPLE_COLUMN]]

if (anyDuplicated(sample_ids)) {
  stop(
    "Duplicate sample IDs were found in the genotype file."
  )
}

geno <- geno_df[
  ,
  setdiff(
    colnames(geno_df),
    SAMPLE_COLUMN
  ),
  drop = FALSE
]

geno <- as.matrix(geno)
storage.mode(geno) <- "numeric"

rownames(geno) <- sample_ids

cat(
  "Genotype matrix:",
  nrow(geno),
  "individuals x",
  ncol(geno),
  "loci\n"
)

cat(
  "Missing genotype values:",
  sum(is.na(geno)),
  "\n"
)


# ==============================================================================
# 5. Remove loci without variation
# ==============================================================================

n_before <- ncol(geno)

geno <- remove_nonvariable_loci(geno)

cat(
  "Removed",
  n_before - ncol(geno),
  "non-variable loci.\n"
)


# ==============================================================================
# 6. Impute missing genotypes
# ==============================================================================

geno_imp <- apply(
  geno,
  2,
  impute_mode
)

if (is.null(dim(geno_imp))) {
  geno_imp <- matrix(
    geno_imp,
    ncol = 1
  )
}

rownames(geno_imp) <- rownames(geno)
colnames(geno_imp) <- colnames(geno)

if (anyNA(geno_imp)) {
  stop(
    "Missing values remain after genotype imputation."
  )
}

cat(
  "Missing genotype values after imputation:",
  sum(is.na(geno_imp)),
  "\n"
)


# ==============================================================================
# 7. Import environmental data and match sample order
# ==============================================================================

env_df <- read.table(
  ENVIRONMENT_FILE,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

if (!SAMPLE_COLUMN %in% colnames(env_df)) {
  stop(
    paste0(
      "Sample column '",
      SAMPLE_COLUMN,
      "' was not found in environmental file."
    )
  )
}

missing_env_vars <- setdiff(
  ENV_VARS,
  colnames(env_df)
)

if (length(missing_env_vars) > 0) {
  stop(
    paste0(
      "Environmental variables not found: ",
      paste(
        missing_env_vars,
        collapse = ", "
      )
    )
  )
}

if (anyDuplicated(env_df[[SAMPLE_COLUMN]])) {
  stop(
    "Duplicate sample IDs were found in the environmental file."
  )
}

match_index <- match(
  rownames(geno_imp),
  env_df[[SAMPLE_COLUMN]]
)

if (anyNA(match_index)) {
  stop(
    "Some genotype samples are missing from the environmental table."
  )
}

env_df <- env_df[
  match_index,
  ,
  drop = FALSE
]

if (!identical(
  rownames(geno_imp),
  env_df[[SAMPLE_COLUMN]]
)) {
  stop(
    "Genotype and environmental sample orders do not match."
  )
}

pred <- env_df[
  ,
  ENV_VARS,
  drop = FALSE
]

pred <- data.frame(
  lapply(
    pred,
    as.numeric
  ),
  check.names = FALSE
)

if (anyNA(pred)) {
  stop(
    "Missing values detected in selected environmental variables."
  )
}


# ==============================================================================
# 8. Environmental correlation analysis
# ==============================================================================

cor_matrix <- cor(
  pred,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.table(
  cor_matrix,
  file = file.path(
    OUTPUT_DIR,
    "environment_correlation_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

high_cor <- which(
  abs(cor_matrix) > COR_THRESHOLD &
    abs(cor_matrix) < 1,
  arr.ind = TRUE
)

if (nrow(high_cor) > 0) {

  high_cor_table <- data.frame(
    Variable1 = rownames(cor_matrix)[high_cor[, 1]],
    Variable2 = colnames(cor_matrix)[high_cor[, 2]],
    Correlation = cor_matrix[high_cor]
  )

  high_cor_table <- high_cor_table[
    !duplicated(
      t(
        apply(
          high_cor_table[, 1:2],
          1,
          sort
        )
      )
    ),
    ]

  write_tsv(
    high_cor_table,
    "highly_correlated_environmental_variables.tsv"
  )

  warning(
    paste0(
      "Environmental predictors with |r| > ",
      COR_THRESHOLD,
      " were detected. Review the output before final analysis."
    )
  )
}


# ==============================================================================
# 9. Environmental PCA
# ==============================================================================

env_pca <- rda(
  pred,
  scale = TRUE
)

env_pca_summary <- summary(env_pca)

env_importance <- as.data.frame(
  env_pca_summary$cont$importance
)

env_importance$Statistic <- rownames(env_importance)

write_tsv(
  env_importance,
  "environment_PCA_importance.tsv"
)

env_loadings <- scores(
  env_pca,
  choices = seq_len(
    min(
      ncol(pred),
      10
    )
  ),
  display = "species",
  scaling = 0
)

env_loadings_df <- data.frame(
  Variable = rownames(env_loadings),
  env_loadings,
  check.names = FALSE
)

write_tsv(
  env_loadings_df,
  "environment_PCA_loadings.tsv"
)

env_scores <- scores(
  env_pca,
  choices = seq_len(
    min(
      ncol(pred),
      10
    )
  ),
  display = "sites",
  scaling = 0
)

env_scores_df <- data.frame(
  Sample = rownames(geno_imp),
  env_scores,
  check.names = FALSE
)

write_tsv(
  env_scores_df,
  "environment_PCA_scores.tsv"
)

save_pdf(
  "environment_PCA_screeplot.pdf",
  screeplot(
    env_pca,
    bstick = TRUE,
    type = "barplot",
    main = "Environmental PCA"
  )
)

pred_PC1 <- env_scores[, 1, drop = FALSE]


# ==============================================================================
# 10. Genetic PCA
# ==============================================================================

gen_pca <- rda(
  geno_imp,
  scale = TRUE
)

save_pdf(
  "genetic_PCA_screeplot.pdf",
  screeplot(
    gen_pca,
    bstick = TRUE,
    type = "barplot",
    main = "Genetic PCA"
  )
)

gen_scores <- scores(
  gen_pca,
  display = "sites",
  scaling = 0
)

gen_scores_df <- data.frame(
  Sample = rownames(geno_imp),
  gen_scores,
  check.names = FALSE
)

write_tsv(
  gen_scores_df,
  "genetic_PCA_scores.tsv"
)


# ==============================================================================
# 11. Prepare LFMM environmental predictors
# ==============================================================================

if (USE_ENV_PC1_FOR_LFMM) {

  X_lfmm <- scale(
    pred_PC1
  )

  colnames(X_lfmm) <- "Env_PC1"

} else {

  X_lfmm <- scale(
    as.matrix(pred)
  )
}

Y_lfmm <- geno_imp


# ==============================================================================
# 12. Run LFMM
# ==============================================================================

cat(
  "Running LFMM with K =",
  LFMM_K,
  "\n"
)

lfmm_fit <- lfmm_ridge(
  Y = Y_lfmm,
  X = X_lfmm,
  K = LFMM_K
)

lfmm_test_result <- lfmm_test(
  Y = Y_lfmm,
  X = X_lfmm,
  lfmm = lfmm_fit,
  calibrate = "gif"
)

gif <- lfmm_test_result$gif

gif_table <- data.frame(
  Predictor = names(gif),
  GIF = as.numeric(gif)
)

write_tsv(
  gif_table,
  "LFMM_genomic_inflation_factor.tsv"
)


# ==============================================================================
# 13. LFMM P-value diagnostics
# ==============================================================================

raw_p <- as.matrix(
  lfmm_test_result$pvalue
)

calibrated_p <- as.matrix(
  lfmm_test_result$calibrated.pvalue
)

z_scores <- as.matrix(
  lfmm_test_result$score
)

if (is.null(dim(raw_p))) {
  raw_p <- matrix(raw_p, ncol = 1)
}

if (is.null(dim(calibrated_p))) {
  calibrated_p <- matrix(calibrated_p, ncol = 1)
}

if (is.null(dim(z_scores))) {
  z_scores <- matrix(z_scores, ncol = 1)
}

colnames(raw_p) <- colnames(X_lfmm)
colnames(calibrated_p) <- colnames(X_lfmm)
colnames(z_scores) <- colnames(X_lfmm)

rownames(raw_p) <- colnames(geno_imp)
rownames(calibrated_p) <- colnames(geno_imp)
rownames(z_scores) <- colnames(geno_imp)


pdf(
  file.path(
    OUTPUT_DIR,
    "LFMM_pvalue_histograms.pdf"
  ),
  width = 10,
  height = 5
)

for (i in seq_len(ncol(raw_p))) {

  par(mfrow = c(1, 2))

  hist(
    raw_p[, i],
    breaks = 50,
    main = paste(
      "Raw P-values:",
      colnames(raw_p)[i]
    ),
    xlab = "P-value"
  )

  hist(
    calibrated_p[, i],
    breaks = 50,
    main = paste(
      "GIF-calibrated P-values:",
      colnames(calibrated_p)[i]
    ),
    xlab = "P-value"
  )
}

dev.off()


pdf(
  file.path(
    OUTPUT_DIR,
    "LFMM_QQ_plots.pdf"
  ),
  width = 6,
  height = 6
)

for (i in seq_len(ncol(calibrated_p))) {

  observed <- sort(
    calibrated_p[, i]
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
    main = colnames(calibrated_p)[i]
  )

  abline(
    0,
    1
  )
}

dev.off()


# ==============================================================================
# 14. LFMM q-values and candidate loci
# ==============================================================================

qvalues <- matrix(
  NA_real_,
  nrow = nrow(calibrated_p),
  ncol = ncol(calibrated_p),
  dimnames = dimnames(calibrated_p)
)

for (i in seq_len(ncol(calibrated_p))) {

  qvalues[, i] <- qvalue(
    calibrated_p[, i]
  )$qvalues
}

lfmm_all_results <- list()
lfmm_candidate_results <- list()

for (i in seq_len(ncol(calibrated_p))) {

  predictor <- colnames(calibrated_p)[i]

  result_i <- data.frame(
    Variant = rownames(calibrated_p),
    Predictor = predictor,
    Z = z_scores[, i],
    Raw_P = raw_p[, i],
    Calibrated_P = calibrated_p[, i],
    Q = qvalues[, i],
    stringsAsFactors = FALSE
  )

  lfmm_all_results[[predictor]] <- result_i

  lfmm_candidate_results[[predictor]] <- result_i[
    result_i$Q < LFMM_FDR,
    ,
    drop = FALSE
  ]
}

lfmm_all_results <- do.call(
  rbind,
  lfmm_all_results
)

lfmm_candidates <- do.call(
  rbind,
  lfmm_candidate_results
)

rownames(lfmm_all_results) <- NULL
rownames(lfmm_candidates) <- NULL

write_tsv(
  lfmm_all_results,
  "LFMM_all_results.tsv"
)

write_tsv(
  lfmm_candidates,
  "LFMM_candidates.tsv"
)

cat(
  "LFMM candidates at FDR <",
  LFMM_FDR,
  ":",
  nrow(lfmm_candidates),
  "\n"
)

lfmm_unique_candidates <- unique(
  lfmm_candidates$Variant
)

writeLines(
  lfmm_unique_candidates,
  file.path(
    OUTPUT_DIR,
    "LFMM_unique_candidate_variants.txt"
  )
)


# ==============================================================================
# 15. Run RDA
# ==============================================================================

rda_data <- pred

rda_model <- rda(
  geno_imp ~ .,
  data = rda_data,
  scale = TRUE
)

capture.output(
  rda_model,
  file = file.path(
    OUTPUT_DIR,
    "RDA_model_summary.txt"
  )
)

r2 <- RsquareAdj(
  rda_model
)

r2_table <- data.frame(
  R2 = r2$r.squared,
  Adjusted_R2 = r2$adj.r.squared
)

write_tsv(
  r2_table,
  "RDA_R2.tsv"
)

vif_values <- vif.cca(
  rda_model
)

vif_table <- data.frame(
  Predictor = names(vif_values),
  VIF = as.numeric(vif_values)
)

write_tsv(
  vif_table,
  "RDA_VIF.tsv"
)

save_pdf(
  "RDA_screeplot.pdf",
  screeplot(
    rda_model,
    main = "RDA constrained axes"
  )
)


# ==============================================================================
# 16. RDA permutation tests
# ==============================================================================

rda_global_test <- anova.cca(
  rda_model,
  permutations = RDA_PERMUTATIONS
)

capture.output(
  rda_global_test,
  file = file.path(
    OUTPUT_DIR,
    "RDA_global_permutation_test.txt"
  )
)

rda_axis_test <- anova.cca(
  rda_model,
  by = "axis",
  permutations = RDA_PERMUTATIONS
)

capture.output(
  rda_axis_test,
  file = file.path(
    OUTPUT_DIR,
    "RDA_axis_permutation_test.txt"
  )
)


# ==============================================================================
# 17. Identify RDA candidate loci
# ==============================================================================

max_axes <- min(
  RDA_N_AXES,
  ncol(
    scores(
      rda_model,
      display = "species"
    )
  )
)

rda_loadings <- scores(
  rda_model,
  choices = seq_len(max_axes),
  display = "species",
  scaling = 0
)

if (is.null(dim(rda_loadings))) {
  rda_loadings <- matrix(
    rda_loadings,
    ncol = 1
  )
}

colnames(rda_loadings) <- paste0(
  "RDA",
  seq_len(ncol(rda_loadings))
)

rownames(rda_loadings) <- colnames(geno_imp)

rda_candidate_list <- list()

for (axis_i in seq_len(ncol(rda_loadings))) {

  idx <- loading_outliers(
    rda_loadings[, axis_i],
    z = RDA_Z_CUTOFF
  )

  axis_name <- colnames(
    rda_loadings
  )[axis_i]

  rda_candidate_list[[axis_name]] <- data.frame(
    Variant = rownames(rda_loadings)[idx],
    Axis = axis_name,
    Loading = rda_loadings[idx, axis_i],
    Z_cutoff = RDA_Z_CUTOFF,
    stringsAsFactors = FALSE
  )
}

rda_candidates <- do.call(
  rbind,
  rda_candidate_list
)

rownames(rda_candidates) <- NULL

write_tsv(
  rda_candidates,
  "RDA_candidates_by_axis.tsv"
)

rda_unique_candidates <- unique(
  rda_candidates$Variant
)

writeLines(
  rda_unique_candidates,
  file.path(
    OUTPUT_DIR,
    "RDA_unique_candidate_variants.txt"
  )
)

cat(
  "Unique RDA candidates:",
  length(rda_unique_candidates),
  "\n"
)


# ==============================================================================
# 18. Save RDA loadings
# ==============================================================================

rda_loading_table <- data.frame(
  Variant = rownames(rda_loadings),
  rda_loadings,
  check.names = FALSE
)

write_tsv(
  rda_loading_table,
  "RDA_variant_loadings.tsv"
)


# ==============================================================================
# 19. Correlations between predictors and RDA axes
# ==============================================================================

axis_cor <- intersetcor(
  rda_model
)

axis_cor <- axis_cor[
  ,
  seq_len(
    min(
      RDA_N_AXES,
      ncol(axis_cor)
    )
  ),
  drop = FALSE
]

axis_cor_table <- data.frame(
  Predictor = rownames(axis_cor),
  axis_cor,
  check.names = FALSE
)

write_tsv(
  axis_cor_table,
  "RDA_environment_axis_correlations.tsv"
)


# ==============================================================================
# 20. Plot RDA candidates
# ==============================================================================

candidate_flag <- colnames(geno_imp) %in% rda_unique_candidates

candidate_border <- ifelse(
  candidate_flag,
  "black",
  NA
)

candidate_fill <- ifelse(
  candidate_flag,
  "red",
  NA
)

pdf(
  file.path(
    OUTPUT_DIR,
    "RDA_candidate_plot_axes1_2.pdf"
  ),
  width = 8,
  height = 7
)

plot(
  rda_model,
  type = "n",
  scaling = 3,
  choices = c(1, 2),
  main = "RDA: axes 1 and 2"
)

points(
  rda_model,
  display = "species",
  pch = 21,
  cex = 0.6,
  col = "gray60",
  bg = "gray90",
  scaling = 3,
  choices = c(1, 2)
)

points(
  rda_model,
  display = "species",
  pch = 21,
  cex = 0.9,
  col = candidate_border,
  bg = candidate_fill,
  scaling = 3,
  choices = c(1, 2)
)

text(
  rda_model,
  scaling = 3,
  display = "bp",
  cex = 0.9,
  choices = c(1, 2)
)

dev.off()


# ==============================================================================
# 21. Optional partial RDA controlling for population structure
# ==============================================================================

if (RUN_PARTIAL_RDA) {

  if (
    max(PARTIAL_RDA_PC_AXES) >
      ncol(gen_scores)
  ) {
    stop(
      "Requested genetic PC does not exist."
    )
  }

  pc_covariates <- as.data.frame(
    gen_scores[
      ,
      PARTIAL_RDA_PC_AXES,
      drop = FALSE
    ]
  )

  colnames(pc_covariates) <- paste0(
    "PC",
    PARTIAL_RDA_PC_AXES
  )

  partial_data <- cbind(
    pred,
    pc_covariates
  )

  predictor_formula <- paste(
    ENV_VARS,
    collapse = " + "
  )

  condition_formula <- paste(
    paste0(
      "PC",
      PARTIAL_RDA_PC_AXES
    ),
    collapse = " + "
  )

  partial_formula <- as.formula(
    paste0(
      "geno_imp ~ ",
      predictor_formula,
      " + Condition(",
      condition_formula,
      ")"
    )
  )

  partial_rda <- rda(
    partial_formula,
    data = partial_data,
    scale = TRUE
  )

  capture.output(
    partial_rda,
    file = file.path(
      OUTPUT_DIR,
      "partial_RDA_model_summary.txt"
    )
  )

  partial_r2 <- RsquareAdj(
    partial_rda
  )

  partial_r2_table <- data.frame(
    R2 = partial_r2$r.squared,
    Adjusted_R2 = partial_r2$adj.r.squared
  )

  write_tsv(
    partial_r2_table,
    "partial_RDA_R2.tsv"
  )

  partial_global_test <- anova.cca(
    partial_rda,
    permutations = RDA_PERMUTATIONS
  )

  capture.output(
    partial_global_test,
    file = file.path(
      OUTPUT_DIR,
      "partial_RDA_global_permutation_test.txt"
    )
  )

  partial_axis_test <- anova.cca(
    partial_rda,
    by = "axis",
    permutations = RDA_PERMUTATIONS
  )

  capture.output(
    partial_axis_test,
    file = file.path(
      OUTPUT_DIR,
      "partial_RDA_axis_permutation_test.txt"
    )
  )
}


# ==============================================================================
# 22. Compare LFMM and RDA candidates
# ==============================================================================

candidate_overlap <- intersect(
  lfmm_unique_candidates,
  rda_unique_candidates
)

lfmm_only <- setdiff(
  lfmm_unique_candidates,
  rda_unique_candidates
)

rda_only <- setdiff(
  rda_unique_candidates,
  lfmm_unique_candidates
)

writeLines(
  candidate_overlap,
  file.path(
    OUTPUT_DIR,
    "LFMM_RDA_overlap.txt"
  )
)

writeLines(
  lfmm_only,
  file.path(
    OUTPUT_DIR,
    "LFMM_only_candidates.txt"
  )
)

writeLines(
  rda_only,
  file.path(
    OUTPUT_DIR,
    "RDA_only_candidates.txt"
  )
)

comparison_summary <- data.frame(
  Category = c(
    "LFMM_unique",
    "RDA_unique",
    "LFMM_RDA_overlap"
  ),
  Number_of_variants = c(
    length(lfmm_unique_candidates),
    length(rda_unique_candidates),
    length(candidate_overlap)
  )
)

write_tsv(
  comparison_summary,
  "candidate_comparison_summary.tsv"
)


# ==============================================================================
# 23. Save run information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

settings <- data.frame(
  Parameter = c(
    "LFMM_K",
    "LFMM_FDR",
    "USE_ENV_PC1_FOR_LFMM",
    "RDA_N_AXES",
    "RDA_Z_CUTOFF",
    "RUN_PARTIAL_RDA",
    "RDA_PERMUTATIONS",
    "COR_THRESHOLD",
    "SEED"
  ),
  Value = c(
    LFMM_K,
    LFMM_FDR,
    USE_ENV_PC1_FOR_LFMM,
    RDA_N_AXES,
    RDA_Z_CUTOFF,
    RUN_PARTIAL_RDA,
    RDA_PERMUTATIONS,
    COR_THRESHOLD,
    SEED
  )
)

write_tsv(
  settings,
  "analysis_settings.tsv"
)

cat(
  "\nAnalysis completed.\nResults written to:",
  OUTPUT_DIR,
  "\n"
)
