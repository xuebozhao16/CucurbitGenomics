# Future Climate Adaptation and Genomic Offset with Gradient Forest

This workflow uses **Gradient Forest (GF)** to model the relationship between genomic variation and the present-day environment, project those relationships onto future climates, and estimate **genomic offset**.

For a location \(s\), let the Gradient Forest transformation of its present environment be

\[
\mathbf{T}_{present}(s)=\left(T_1(s),T_2(s),\ldots,T_p(s)\right)
\]

and the transformation of its future environment be

\[
\mathbf{T}_{future}(s)=\left(T'_1(s),T'_2(s),\ldots,T'_p(s)\right).
\]

Genomic offset is calculated as

\[
GO(s)=
\sqrt{
\sum_{j=1}^{p}
\left[T'_j(s)-T_j(s)\right]^2
}.
\]

Thus:

```text
low genomic offset
    -> future climate is relatively similar in adaptive genomic space
    -> comparatively little genomic change is predicted to be required

high genomic offset
    -> large mismatch between present adaptive genomic composition and future climate
    -> comparatively large genomic change is predicted to be required
```

Genomic offset is a measure of predicted **genotype–environment mismatch**, not a direct estimate of extinction probability.

## Data preparation

The main population-level input table should contain:

```text
Population
Longitude
Latitude
environmental variables
adaptive/candidate variant allele frequencies
optional reference/neutral variant allele frequencies
```

For example:

```text
Population  Longitude  Latitude  bio_1  bio_3  bio_5  bio_9  bio_12  SNP_001  SNP_002  SV_001
Pop01       ...        ...       ...    ...    ...    ...    ...     0.15     0.72     0.10
Pop02       ...        ...       ...    ...    ...    ...    ...     0.23     0.64     0.18
Pop03       ...        ...       ...    ...    ...    ...    ...     0.81     0.12     0.76
```

Allele-frequency responses should normally satisfy

\[
0 \leq p \leq 1.
\]

For diploid SNPs or biallelic SVs,

\[
p_{ALT}=\frac{N_{ALT}}{2N_{called}}.
\]

Candidate variants may come from previous GEA analyses such as:

```text
LFMM
RDA
BayPass
Bayenv2
```

Present and future climate data must contain the **same environmental variables**, with identical names and units.

For population-level prediction:

```text
Population  bio_1  bio_3  bio_5  bio_9  bio_12
Pop01       ...
Pop02       ...
```

For landscape projection, use aligned raster layers:

```text
present/
    bio_1.tif
    bio_3.tif
    bio_5.tif
    bio_9.tif
    bio_12.tif

future_2050/
    bio_1.tif
    bio_3.tif
    bio_5.tif
    bio_9.tif
    bio_12.tif

future_2070/
    bio_1.tif
    bio_3.tif
    bio_5.tif
    bio_9.tif
    bio_12.tif
```

All rasters should use the same CRS, extent, and resolution.

## Gradient Forest

Gradient Forest models nonlinear changes in allele frequency along environmental gradients. For variant \(l\),

\[
p_l=f_l(E_1,E_2,\ldots,E_p)+\varepsilon_l.
\]

Load the population table:

```r
gf_data <- read.table(
    "gradient_forest_input.tsv",
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
)
```

Define environmental predictors:

```r
ENV_VARS <- c(
    "bio_1",
    "bio_3",
    "bio_5",
    "bio_9",
    "bio_12"
)
```

Define candidate variants:

```r
candidate_vars <- grep(
    "^(SNP|SV)_",
    colnames(gf_data),
    value = TRUE
)
```

Remove invariant responses:

```r
candidate_vars <- candidate_vars[
    sapply(
        gf_data[, candidate_vars, drop = FALSE],
        function(x) sd(x, na.rm = TRUE) > 0
    )
]
```

Create the GF matrix:

```r
gf_candidate_data <- cbind(
    gf_data[, ENV_VARS, drop = FALSE],
    gf_data[, candidate_vars, drop = FALSE]
)
```

Calculate tree depth:

```r
maxLevel <- log2(
    0.368 * nrow(gf_data) / 2
)
```

Fit Gradient Forest:

```r
library(gradientForest)

gf_candidate <- gradientForest(
    gf_candidate_data,
    predictor.vars = ENV_VARS,
    response.vars = candidate_vars,
    ntree = 500,
    maxLevel = maxLevel,
    trace = TRUE,
    corr.threshold = 0.50
)
```

Environmental importance:

```r
env_importance <- sort(
    gf_candidate$overall.imp,
    decreasing = TRUE
)

most_important_env <- names(env_importance)[1]
```

## Allele-turnover functions

Overall turnover along the most important environmental gradient:

```r
overall_turnover <- cumimp(
    gf_candidate,
    predictor = most_important_env,
    type = "Overall",
    standardize = TRUE
)
```

Variant-specific turnover:

```r
variant_turnover <- cumimp(
    gf_candidate,
    predictor = most_important_env,
    type = "Species",
    standardize = TRUE
)
```

Transform the sampled populations:

```r
present_env_population <- gf_data[
    ,
    ENV_VARS,
    drop = FALSE
]

population_turnover <- predict(
    gf_candidate,
    present_env_population
)
```

A simple descriptive classification of contrasting adaptive environments is:

```r
turnover_threshold <- mean(
    population_turnover[, most_important_env],
    na.rm = TRUE
)

low_environment_group <- rownames(population_turnover)[
    population_turnover[, most_important_env] < turnover_threshold
]

high_environment_group <- rownames(population_turnover)[
    population_turnover[, most_important_env] >= turnover_threshold
]
```

## Future projection and genomic offset

Load population-level present and future environments:

```r
present_env <- read.table(
    "present_population_environment.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)

future_2050_env <- read.table(
    "future_2050_population_environment.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)

future_2070_env <- read.table(
    "future_2070_population_environment.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)
```

Match variables and population order:

```r
present_env <- present_env[, ENV_VARS, drop = FALSE]

future_2050_env <- future_2050_env[
    rownames(present_env),
    ENV_VARS,
    drop = FALSE
]

future_2070_env <- future_2070_env[
    rownames(present_env),
    ENV_VARS,
    drop = FALSE
]
```

Transform all environments using the **same fitted GF model**:

```r
turn_present <- predict(
    gf_candidate,
    present_env
)

turn_2050 <- predict(
    gf_candidate,
    future_2050_env
)

turn_2070 <- predict(
    gf_candidate,
    future_2070_env
)
```

Calculate genomic offset:

```r
euclidean_offset <- function(present, future) {
    sqrt(
        rowSums(
            (future - present)^2
        )
    )
}

offset_2050 <- euclidean_offset(
    turn_present,
    turn_2050
)

offset_2070 <- euclidean_offset(
    turn_present,
    turn_2070
)
```

Save:

```r
population_offset <- data.frame(
    Population = rownames(present_env),
    Genomic_offset_2050 = offset_2050,
    Genomic_offset_2070 = offset_2070
)

write.table(
    population_offset,
    "population_genomic_offset.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
```

For raster analysis:

```text
present climate rasters
        ↓
environment matrix for every cell
        ↓
GF transformation
        ↓
T_present

future climate rasters
        ↓
same GF transformation
        ↓
T_future

T_present vs T_future
        ↓
Euclidean distance
        ↓
genomic offset raster
```

The same fitted present-day GF model should be used for all future scenarios.

## SNP and SV analysis

Run marker classes separately and together:

```text
SNP-only GF
SV-only GF
SNP + SV GF
```

This allows comparison of:

```text
SNP genomic offset
SV genomic offset
combined adaptive genomic offset
```

When one marker class greatly outnumbers another, the combined model may be dominated by the larger class. Reporting all three analyses is therefore useful.

## Multiple future scenarios

The model can be projected to multiple scenarios:

```text
2050 SSP1-2.6
2050 SSP2-4.5
2050 SSP5-8.5
2070 SSP1-2.6
2070 SSP2-4.5
2070 SSP5-8.5
```

For scenario \(c\),

\[
GO_c(s)=
\sqrt{
\sum_j
\left[
T_{future,c,j}(s)-T_{present,j}(s)
\right]^2
}.
\]

A final population table can contain:

```text
Population
Longitude
Latitude
GO_2050_SSP126
GO_2050_SSP245
GO_2050_SSP585
GO_2070_SSP126
GO_2070_SSP245
GO_2070_SSP585
```
