#!/usr/bin/env Rscript

# ==============================================================================
# Future Climate Adaptation and Genomic Offset with Gradient Forest
# ==============================================================================

# -----------------------------
# 1. User settings
# -----------------------------

GF_INPUT <- "gradient_forest_input.tsv"

PRESENT_POP_ENV <- "present_population_environment.tsv"
FUTURE_2050_POP_ENV <- "future_2050_population_environment.tsv"
FUTURE_2070_POP_ENV <- "future_2070_population_environment.tsv"

PRESENT_RASTER_DIR <- "present"
FUTURE_2050_RASTER_DIR <- "future_2050"
FUTURE_2070_RASTER_DIR <- "future_2070"

OUTPUT_DIR <- "Future_climate_adaptation_results"

POP_COLUMN <- "Population"
LON_COLUMN <- "Longitude"
LAT_COLUMN <- "Latitude"

ENV_VARS <- c(
  "bio_1",
  "bio_3",
  "bio_5",
  "bio_9",
  "bio_12"
)

# Variant prefixes in the GF input table.
SNP_PATTERN <- "^SNP_"
SV_PATTERN <- "^SV_"

# Which marker sets to run.
RUN_SNP <- TRUE
RUN_SV <- TRUE
RUN_COMBINED <- TRUE

NTREE <- 500
CORR_THRESHOLD <- 0.50

# Raster filename extension.
RASTER_PATTERN <- "\\.tif$"

set.seed(2026)


# -----------------------------
# 2. Packages
# -----------------------------

required_packages <- c(
  "gradientForest",
  "terra"
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
    paste(
      "Missing packages:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(gradientForest)
  library(terra)
})

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------
# 3. Helper functions
# -----------------------------

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


check_environment_columns <- function(dat, env_vars, object_name) {

  missing <- setdiff(
    env_vars,
    colnames(dat)
  )

  if (length(missing) > 0) {
    stop(
      paste0(
        object_name,
        " is missing environmental variables: ",
        paste(missing, collapse = ", ")
      )
    )
  }
}


filter_variant_columns <- function(dat, variants) {

  if (length(variants) == 0) {
    return(character(0))
  }

  keep <- vapply(
    dat[, variants, drop = FALSE],
    function(x) {

      x <- as.numeric(x)

      all(
        is.finite(x[!is.na(x)])
      ) &&
        sum(!is.na(x)) >= 3 &&
        sd(x, na.rm = TRUE) > 0
    },
    logical(1)
  )

  variants[keep]
}


complete_gf_rows <- function(dat, env_vars, response_vars) {

  cols <- c(
    env_vars,
    response_vars
  )

  dat[
    complete.cases(
      dat[, cols, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]
}


fit_gradient_forest <- function(
    dat,
    env_vars,
    response_vars,
    ntree = 500,
    corr_threshold = 0.50
) {

  response_vars <- filter_variant_columns(
    dat,
    response_vars
  )

  if (length(response_vars) == 0) {
    stop(
      "No variable response loci remain after filtering."
    )
  }

  dat2 <- complete_gf_rows(
    dat,
    env_vars,
    response_vars
  )

  if (nrow(dat2) < 5) {
    stop(
      "Too few complete populations remain for Gradient Forest."
    )
  }

  max_level <- log2(
    0.368 * nrow(dat2) / 2
  )

  gf_matrix <- cbind(
    dat2[, env_vars, drop = FALSE],
    dat2[, response_vars, drop = FALSE]
  )

  model <- gradientForest(
    gf_matrix,
    predictor.vars = env_vars,
    response.vars = response_vars,
    ntree = ntree,
    maxLevel = max_level,
    trace = TRUE,
    corr.threshold = corr_threshold
  )

  list(
    model = model,
    data = dat2,
    response_vars = response_vars,
    maxLevel = max_level
  )
}


euclidean_offset <- function(present, future) {

  if (!all(dim(present) == dim(future))) {
    stop(
      "Present and future transformed matrices have different dimensions."
    )
  }

  sqrt(
    rowSums(
      (future - present)^2
    )
  )
}


load_population_environment <- function(
    filename,
    pop_column,
    env_vars
) {

  x <- read.table(
    filename,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!pop_column %in% colnames(x)) {
    stop(
      paste0(
        "Population column not found in ",
        filename
      )
    )
  }

  check_environment_columns(
    x,
    env_vars,
    filename
  )

  rownames(x) <- x[[pop_column]]

  x[
    ,
    env_vars,
    drop = FALSE
  ]
}


match_future_to_present <- function(
    present,
    future
) {

  missing <- setdiff(
    rownames(present),
    rownames(future)
  )

  if (length(missing) > 0) {
    stop(
      paste(
        "Future environment is missing populations:",
        paste(missing, collapse = ", ")
      )
    )
  }

  future[
    rownames(present),
    colnames(present),
    drop = FALSE
  ]
}


read_climate_stack <- function(
    directory,
    env_vars,
    pattern = "\\.tif$"
) {

  files <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop(
      paste0(
        "No raster files found in: ",
        directory
      )
    )
  }

  r <- rast(files)

  # Remove extensions from layer names where needed.
  names(r) <- tools::file_path_sans_ext(
    basename(files)
  )

  missing <- setdiff(
    env_vars,
    names(r)
  )

  if (length(missing) > 0) {
    stop(
      paste0(
        "Raster directory ",
        directory,
        " is missing: ",
        paste(missing, collapse = ", ")
      )
    )
  }

  r[[env_vars]]
}


align_future_raster <- function(
    future,
    present
) {

  if (!same.crs(future, present)) {
    future <- project(
      future,
      present
    )
  }

  if (
    !all(
      res(future) == res(present)
    ) ||
      !all(
        ext(future) == ext(present)
      )
  ) {

    future <- resample(
      future,
      present,
      method = "bilinear"
    )
  }

  future
}


raster_to_env_df <- function(r, env_vars) {

  vals <- as.data.frame(
    r,
    cells = TRUE,
    xy = TRUE,
    na.rm = FALSE
  )

  keep <- complete.cases(
    vals[, env_vars, drop = FALSE]
  )

  vals[
    keep,
    c(
      "cell",
      "x",
      "y",
      env_vars
    ),
    drop = FALSE
  ]
}


project_raster_offset <- function(
    gf_model,
    present_raster,
    future_raster,
    env_vars,
    output_file
) {

  future_raster <- align_future_raster(
    future_raster,
    present_raster
  )

  present_df <- raster_to_env_df(
    present_raster,
    env_vars
  )

  future_vals <- as.data.frame(
    future_raster,
    cells = TRUE,
    xy = TRUE,
    na.rm = FALSE
  )

  future_df <- future_vals[
    match(
      present_df$cell,
      future_vals$cell
    ),
    ,
    drop = FALSE
  ]

  valid <- complete.cases(
    future_df[, env_vars, drop = FALSE]
  )

  present_df_valid <- present_df[
    valid,
    ,
    drop = FALSE
  ]

  future_df_valid <- future_df[
    valid,
    ,
    drop = FALSE
  ]

  t_present <- predict(
    gf_model,
    present_df_valid[, env_vars, drop = FALSE]
  )

  t_future <- predict(
    gf_model,
    future_df_valid[, env_vars, drop = FALSE]
  )

  offset <- euclidean_offset(
    t_present,
    t_future
  )

  out <- present_raster[[1]]
  values(out) <- NA_real_

  out[
    present_df_valid$cell
  ] <- offset

  names(out) <- "genomic_offset"

  writeRaster(
    out,
    file.path(
      OUTPUT_DIR,
      output_file
    ),
    overwrite = TRUE
  )

  out
}


run_marker_set <- function(
    label,
    response_vars,
    gf_data,
    present_pop_env,
    future_2050_pop_env,
    future_2070_pop_env,
    present_raster = NULL,
    future_2050_raster = NULL,
    future_2070_raster = NULL
) {

  message(
    "\nRunning marker set: ",
    label
  )

  fitted <- fit_gradient_forest(
    dat = gf_data,
    env_vars = ENV_VARS,
    response_vars = response_vars,
    ntree = NTREE,
    corr_threshold = CORR_THRESHOLD
  )

  gf_model <- fitted$model

  saveRDS(
    gf_model,
    file.path(
      OUTPUT_DIR,
      paste0(
        label,
        "_gradient_forest_model.rds"
      )
    )
  )

  importance <- sort(
    gf_model$overall.imp,
    decreasing = TRUE
  )

  importance_table <- data.frame(
    Environment = names(importance),
    Importance = as.numeric(importance)
  )

  write_tsv(
    importance_table,
    paste0(
      label,
      "_environment_importance.tsv"
    )
  )

  most_important_env <- importance_table$Environment[1]

  turnover_overall <- cumimp(
    gf_model,
    predictor = most_important_env,
    type = "Overall",
    standardize = TRUE
  )

  saveRDS(
    turnover_overall,
    file.path(
      OUTPUT_DIR,
      paste0(
        label,
        "_overall_turnover.rds"
      )
    )
  )

  present_pop <- present_pop_env
  future50_pop <- match_future_to_present(
    present_pop,
    future_2050_pop_env
  )

  future70_pop <- match_future_to_present(
    present_pop,
    future_2070_pop_env
  )

  turn_present <- predict(
    gf_model,
    present_pop
  )

  turn_2050 <- predict(
    gf_model,
    future50_pop
  )

  turn_2070 <- predict(
    gf_model,
    future70_pop
  )

  offset_2050 <- euclidean_offset(
    turn_present,
    turn_2050
  )

  offset_2070 <- euclidean_offset(
    turn_present,
    turn_2070
  )

  pop_result <- data.frame(
    Population = rownames(present_pop),
    Genomic_offset_2050 = offset_2050,
    Genomic_offset_2070 = offset_2070
  )

  write_tsv(
    pop_result,
    paste0(
      label,
      "_population_genomic_offset.tsv"
    )
  )

  # Optional descriptive adaptive-environment grouping.
  pop_turnover <- predict(
    gf_model,
    present_pop
  )

  threshold <- mean(
    pop_turnover[, most_important_env],
    na.rm = TRUE
  )

  adaptive_group <- ifelse(
    pop_turnover[, most_important_env] >= threshold,
    "high_environment",
    "low_environment"
  )

  adaptive_table <- data.frame(
    Population = rownames(pop_turnover),
    Most_important_environment = most_important_env,
    Turnover_score = pop_turnover[, most_important_env],
    Adaptive_group = adaptive_group
  )

  write_tsv(
    adaptive_table,
    paste0(
      label,
      "_adaptive_environment_groups.tsv"
    )
  )

  raster_2050 <- NULL
  raster_2070 <- NULL

  if (
    !is.null(present_raster) &&
      !is.null(future_2050_raster)
  ) {

    raster_2050 <- project_raster_offset(
      gf_model = gf_model,
      present_raster = present_raster,
      future_raster = future_2050_raster,
      env_vars = ENV_VARS,
      output_file = paste0(
        label,
        "_genomic_offset_2050.tif"
      )
    )
  }

  if (
    !is.null(present_raster) &&
      !is.null(future_2070_raster)
  ) {

    raster_2070 <- project_raster_offset(
      gf_model = gf_model,
      present_raster = present_raster,
      future_raster = future_2070_raster,
      env_vars = ENV_VARS,
      output_file = paste0(
        label,
        "_genomic_offset_2070.tif"
      )
    )
  }

  list(
    label = label,
    model = gf_model,
    importance = importance_table,
    most_important_env = most_important_env,
    population_offset = pop_result,
    adaptive_groups = adaptive_table,
    raster_2050 = raster_2050,
    raster_2070 = raster_2070
  )
}


# -----------------------------
# 4. Read GF input
# -----------------------------

gf_data <- read.table(
  GF_INPUT,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_meta <- c(
  POP_COLUMN,
  LON_COLUMN,
  LAT_COLUMN
)

missing_meta <- setdiff(
  required_meta,
  colnames(gf_data)
)

if (length(missing_meta) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_meta, collapse = ", ")
    )
  )
}

check_environment_columns(
  gf_data,
  ENV_VARS,
  GF_INPUT
)

if (anyDuplicated(gf_data[[POP_COLUMN]])) {
  stop(
    "Population IDs in GF input must be unique."
  )
}

rownames(gf_data) <- gf_data[[POP_COLUMN]]


# -----------------------------
# 5. Identify SNP and SV columns
# -----------------------------

snp_vars <- grep(
  SNP_PATTERN,
  colnames(gf_data),
  value = TRUE
)

sv_vars <- grep(
  SV_PATTERN,
  colnames(gf_data),
  value = TRUE
)

message(
  "SNP variants detected: ",
  length(snp_vars)
)

message(
  "SV variants detected: ",
  length(sv_vars)
)


# -----------------------------
# 6. Load population environments
# -----------------------------

present_pop_env <- load_population_environment(
  PRESENT_POP_ENV,
  POP_COLUMN,
  ENV_VARS
)

future_2050_pop_env <- load_population_environment(
  FUTURE_2050_POP_ENV,
  POP_COLUMN,
  ENV_VARS
)

future_2070_pop_env <- load_population_environment(
  FUTURE_2070_POP_ENV,
  POP_COLUMN,
  ENV_VARS
)


# -----------------------------
# 7. Load climate rasters if available
# -----------------------------

present_raster <- NULL
future_2050_raster <- NULL
future_2070_raster <- NULL

if (
  dir.exists(PRESENT_RASTER_DIR) &&
    dir.exists(FUTURE_2050_RASTER_DIR) &&
    dir.exists(FUTURE_2070_RASTER_DIR)
) {

  present_raster <- read_climate_stack(
    PRESENT_RASTER_DIR,
    ENV_VARS,
    RASTER_PATTERN
  )

  future_2050_raster <- read_climate_stack(
    FUTURE_2050_RASTER_DIR,
    ENV_VARS,
    RASTER_PATTERN
  )

  future_2070_raster <- read_climate_stack(
    FUTURE_2070_RASTER_DIR,
    ENV_VARS,
    RASTER_PATTERN
  )
}


# -----------------------------
# 8. Run analyses
# -----------------------------

results <- list()

if (RUN_SNP && length(snp_vars) > 0) {

  results$SNP <- run_marker_set(
    label = "SNP",
    response_vars = snp_vars,
    gf_data = gf_data,
    present_pop_env = present_pop_env,
    future_2050_pop_env = future_2050_pop_env,
    future_2070_pop_env = future_2070_pop_env,
    present_raster = present_raster,
    future_2050_raster = future_2050_raster,
    future_2070_raster = future_2070_raster
  )
}

if (RUN_SV && length(sv_vars) > 0) {

  results$SV <- run_marker_set(
    label = "SV",
    response_vars = sv_vars,
    gf_data = gf_data,
    present_pop_env = present_pop_env,
    future_2050_pop_env = future_2050_pop_env,
    future_2070_pop_env = future_2070_pop_env,
    present_raster = present_raster,
    future_2050_raster = future_2050_raster,
    future_2070_raster = future_2070_raster
  )
}

combined_vars <- c(
  snp_vars,
  sv_vars
)

if (
  RUN_COMBINED &&
    length(combined_vars) > 0
) {

  results$Combined <- run_marker_set(
    label = "SNP_SV_combined",
    response_vars = combined_vars,
    gf_data = gf_data,
    present_pop_env = present_pop_env,
    future_2050_pop_env = future_2050_pop_env,
    future_2070_pop_env = future_2070_pop_env,
    present_raster = present_raster,
    future_2050_raster = future_2050_raster,
    future_2070_raster = future_2070_raster
  )
}


# -----------------------------
# 9. Combine population offsets
# -----------------------------

offset_tables <- lapply(
  results,
  function(x) {

    y <- x$population_offset

    colnames(y)[-1] <- paste0(
      x$label,
      "_",
      colnames(y)[-1]
    )

    y
  }
)

if (length(offset_tables) > 0) {

  combined_offset <- Reduce(
    function(x, y) {
      merge(
        x,
        y,
        by = "Population",
        all = TRUE
      )
    },
    offset_tables
  )

  coords <- gf_data[
    ,
    c(
      POP_COLUMN,
      LON_COLUMN,
      LAT_COLUMN
    ),
    drop = FALSE
  ]

  colnames(coords)[1] <- "Population"

  combined_offset <- merge(
    coords,
    combined_offset,
    by = "Population",
    all.y = TRUE
  )

  write_tsv(
    combined_offset,
    "all_marker_sets_population_genomic_offset.tsv"
  )
}


# -----------------------------
# 10. Save complete R objects
# -----------------------------

saveRDS(
  results,
  file.path(
    OUTPUT_DIR,
    "future_climate_adaptation_results.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

message(
  "\nAnalysis completed. Results written to: ",
  OUTPUT_DIR
)
