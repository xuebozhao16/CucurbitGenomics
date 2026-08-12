# Introgression and Gene Flow Analysis Pipeline

This workflow summarizes the main analyses used to investigate **population admixture, introgression, and excess allele sharing** from population genomic SNP data. The main tools are **ADMIXTOOLS / admixr, ADMIXTOOLS 2, and TreeMix**.

The general workflow is:

```text
Filtered SNP VCF
      |
      +--> ADMIXTOOLS / ADMIXTOOLS 2
      |       |
      |       +--> D-statistic
      |       +--> f3
      |       +--> f4
      |       +--> f4-ratio
      |       +--> qpWave
      |       +--> qpAdm
      |       +--> qpGraph
      |
      +--> TreeMix
              |
              +--> population tree
              +--> migration edges
              +--> residuals
```

TreeMix is useful for exploring candidate migration relationships, while f-statistics provide formal tests of allele-sharing asymmetry and admixture.

---

## 1. Prepare SNP data

Start from a high-quality multisample VCF:

```text
population.filtered.vcf.gz
```

Index the VCF:

```bash
bcftools index population.filtered.vcf.gz
```

Check sample names:

```bash
bcftools query -l population.filtered.vcf.gz
```

Retain biallelic SNPs:

```bash
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    population.filtered.vcf.gz \
    -Oz \
    -o population.snps.vcf.gz

bcftools index population.snps.vcf.gz
```

Prepare a population assignment file:

```text
sample01    PopA
sample02    PopA
sample03    PopB
sample04    PopB
sample05    PopC
sample06    Outgroup
```

Population labels should be consistent across all downstream analyses.

---

## 2. ADMIXTOOLS / admixr

Classical ADMIXTOOLS commonly uses EIGENSTRAT input:

```text
population.geno
population.snp
population.ind
```

Convert the VCF using a tested VCF-to-EIGENSTRAT script:

```bash
sh convertVCFtoEigenstrat.sh population.snps.vcf.gz
```

Check the resulting `.ind` file carefully to ensure that individuals have the correct population labels.

Load the dataset in R:

```r
library(admixr)
library(dplyr)

data_prefix <- "population"

snps <- eigenstrat(data_prefix)

snp_info <- count_snps(snps)

snp_info
```

---

## 3. D-statistic / ABBA-BABA

For four populations arranged approximately as:

```text
((W, X), Y), Z
```

where `Z` is the outgroup, the two discordant genealogical patterns are:

```text
ABBA
W = A
X = B
Y = B
Z = A
```

and:

```text
BABA
W = B
X = A
Y = B
Z = A
```

For one haploid observation per population, the D-statistic can be written as:

$$
D = \frac{n_{\mathrm{BABA}}-n_{\mathrm{ABBA}}}
         {n_{\mathrm{BABA}}+n_{\mathrm{ABBA}}}
$$

For population allele frequencies, ADMIXTOOLS uses the allele-frequency formulation. At SNP \(i\), let the derived-allele frequencies in populations \(W,X,Y,Z\) be \(w_i,x_i,y_i,z_i\). Then:

$$
N_i=(w_i-x_i)(y_i-z_i)
$$

and

$$
D_i=(w_i+x_i-2w_ix_i)(y_i+z_i-2y_iz_i)
$$

The genome-wide statistic is:

$$
D(W,X;Y,Z)=
\frac{\sum_i (w_i-x_i)(y_i-z_i)}
{\sum_i (w_i+x_i-2w_ix_i)(y_i+z_i-2y_iz_i)}
$$

Under a simple null model with no asymmetric gene flow:

$$
D \approx 0
$$

ADMIXTOOLS estimates uncertainty using a block jackknife. The Z-score is:

$$
Z=\frac{D}{SE(D)}
$$

A commonly used descriptive threshold is:

$$
|Z|\geq 3
$$

Under the ADMIXTOOLS BABA-minus-ABBA convention:

```text
D > 0
```

indicates excess allele sharing between `W` and `Y`, whereas:

```text
D < 0
```

indicates excess allele sharing between `X` and `Y`.

Always verify the sign convention used by the specific software before biological interpretation.

Run one test:

```r
d1 <- d(
    W = "PopA",
    X = "PopB",
    Y = "PopC",
    Z = "Outgroup",
    data = snps
)

d1
```

Run multiple populations:

```r
test_pops <- c(
    "PopA1",
    "PopA2",
    "PopA3"
)

d_many <- d(
    W = test_pops,
    X = "PopB",
    Y = "PopC",
    Z = "Outgroup",
    data = snps
)

d_many
```

For interpretation, always record the exact quartet:

```text
D(W, X; Y, Z)
```

rather than reporting only the numerical D value.

---

## 4. f3 statistic

The basic f2 statistic between populations \(A\) and \(B\) is:

$$
f_2(A,B)=\frac{1}{M}\sum_{j=1}^{M}(a_j-b_j)^2
$$

where \(a_j\) and \(b_j\) are allele frequencies at SNP \(j\).

The f3 statistic can be expressed as:

$$
f_3(A;B,C)
=
\frac{1}{M}
\sum_{j=1}^{M}
(a_j-b_j)(a_j-c_j)
$$

or equivalently:

$$
f_3(A;B,C)
=
\frac{1}{2}
\left[
f_2(A,B)+f_2(A,C)-f_2(B,C)
\right]
$$

For an admixture f3 test:

```text
f3(Target; Source1, Source2)
```

a significantly negative value:

$$
f_3(Target;Source1,Source2)<0
$$

provides evidence that the target population contains ancestry related to the two source populations.

In ADMIXTOOLS 2:

```r
library(admixtools)
library(tidyverse)

prefix <- "population"
f2_dir <- "f2_data"

extract_f2(
    prefix,
    f2_dir
)

f2_blocks <- f2_from_precomp(
    f2_dir
)

count_snps(f2_blocks)
```

Run admixture f3:

```r
f3_result <- f3(
    f2_blocks,
    "Target",
    "Source1",
    "Source2"
)

f3_result
```

Outgroup f3 can be used to measure shared drift:

$$
f_3(O;A,B)
$$

where \(O\) is an outgroup.

Example:

```r
outgroup_f3 <- f3(
    f2_blocks,
    "Outgroup",
    "PopA",
    "PopB"
)

outgroup_f3
```

---

## 5. f4 statistic

The f4 statistic is:

$$
f_4(A,B;C,D)
=
\frac{1}{M}
\sum_{j=1}^{M}
(a_j-b_j)(c_j-d_j)
$$

It can also be expressed using f2 statistics:

$$
f_4(A,B;C,D)
=
\frac{1}{2}
\left[
f_2(A,D)+f_2(B,C)-f_2(A,C)-f_2(B,D)
\right]
$$

If populations \(A\) and \(B\) form a clade relative to \(C\) and \(D\), then:

$$
f_4(A,B;C,D)\approx0
$$

A significantly non-zero value indicates asymmetric shared drift.

Using `admixr`:

```r
f4_result <- f4(
    W = "PopA",
    X = "PopB",
    Y = "PopC",
    Z = "Outgroup",
    data = snps
)

f4_result
```

Using ADMIXTOOLS 2:

```r
f4_result <- f4(
    f2_blocks,
    "PopA",
    "PopB",
    "PopC",
    "Outgroup"
)

f4_result
```

The corresponding Z-score is:

$$
Z=\frac{f_4}{SE(f_4)}
$$

---

## 6. f4-ratio

The f4-ratio can estimate an admixture proportion under an appropriate population topology.

For a target population \(X\), the admixture proportion can be written conceptually as:

$$
\alpha=
\frac{f_4(A,O;X,C)}
     {f_4(A,O;B,C)}
$$

where the interpretation of \(A,B,C,O\) depends on the assumed admixture graph.

Using `admixr`:

```r
f4r <- f4ratio(
    X = "Target",
    A = "SourceA",
    B = "SourceB",
    C = "Sister",
    O = "Outgroup",
    data = snps
)

f4r
```

For multiple target populations:

```r
targets <- c(
    "Target1",
    "Target2",
    "Target3"
)

f4r_many <- f4ratio(
    X = targets,
    A = "SourceA",
    B = "SourceB",
    C = "Sister",
    O = "Outgroup",
    data = snps
)

f4r_many
```

The estimated \(\alpha\) should only be interpreted as an ancestry proportion when the assumed topology and source proxies are appropriate.

---

## 7. qpWave and qpAdm

qpWave evaluates whether a set of left populations can be explained by a given number of ancestry streams relative to a set of right populations.

Define:

```r
left <- c(
    "Source1",
    "Source2"
)

right <- c(
    "Outgroup1",
    "Outgroup2",
    "Outgroup3",
    "Outgroup4"
)
```

Run qpWave:

```r
qpwave_result <- qpwave(
    f2_blocks,
    left,
    right
)

qpwave_result$rankdrop
qpwave_result$f4
```

The method evaluates the rank of a matrix of f4 statistics. If the matrix has rank \(r\), the populations require at least:

$$
r+1
$$

independent ancestry streams relative to the selected right populations.

qpAdm uses the same f-statistic framework but models a target population as a mixture of candidate source populations.

Define:

```r
target <- "Target"
```

Run:

```r
qpadm_result <- qpadm(
    f2_blocks,
    left,
    right,
    target
)

qpadm_result$weights
qpadm_result$popdrop
qpadm_result$f4
```

For a two-source model:

$$
Target
=
\alpha Source_1
+
(1-\alpha)Source_2
$$

with:

$$
\alpha_1+\alpha_2=1
$$

For \(K\) sources:

$$
Target=\sum_{k=1}^{K}\alpha_k Source_k
$$

subject to:

$$
\sum_{k=1}^{K}\alpha_k=1
$$

The qpAdm output should be evaluated using the model fit, ancestry weights, standard errors, feasibility of the weights, and biological plausibility of the source and right populations.

To compare multiple possible sources:

```r
candidate_pops <- c(
    "Source1",
    "Source2",
    "Source3",
    "Source4"
)

fixed_outgroups <- c(
    "Outgroup1",
    "Outgroup2"
)

rotate_result <- qpadm_rotate(
    f2_blocks,
    leftright = candidate_pops,
    target = "Target",
    rightfix = fixed_outgroups
)

rotate_result
```

---

## 8. qpGraph

qpGraph fits an explicit admixture graph to the observed f-statistics.

The graph consists of population splits, drift branches, and admixture edges.

For an admixed population \(X\):

$$
p_X=\alpha p_A+(1-\alpha)p_B
$$

where \(p_A\) and \(p_B\) represent allele frequencies inherited from two ancestral branches.

Fit a graph:

```r
graph_result <- qpgraph(
    f2_blocks,
    my_graph,
    return_fstats = TRUE
)

graph_result$score
graph_result$edges
graph_result$worst_residual
```

Plot:

```r
plot_graph(
    graph_result$edges
)
```

or:

```r
plotly_graph(
    graph_result$edges
)
```

The important criterion is not simply whether the graph can be drawn, but whether the fitted f-statistics adequately reproduce the observed f-statistics.

For each fitted statistic:

$$
Z_{\mathrm{residual}}
=
\frac{f_{\mathrm{observed}}-f_{\mathrm{fitted}}}
     {SE}
$$

Large absolute residual Z-scores indicate parts of the population history that are not adequately explained by the proposed graph.

---

## 9. TreeMix

TreeMix uses population allele frequencies to infer a maximum-likelihood population tree and then adds migration edges to account for covariance not explained by the tree.

For populations \(i\) and \(j\), the model is based on covariance in allele-frequency drift. Conceptually:

$$
\hat{W}_{ij}
=
\frac{1}{M}
\sum_{k=1}^{M}
(\hat{p}_{ik}-\hat{p}_{Ak})
(\hat{p}_{jk}-\hat{p}_{Ak})
$$

where \(p_{ik}\) is the allele frequency at SNP \(k\), and \(p_A\) represents an ancestral/reference allele frequency used in the drift model.

TreeMix fits a graph that attempts to minimize the discrepancy between observed and model-predicted covariance.

First remove highly incomplete sites. A strict complete-data dataset can be generated with:

```bash
FILE="population.snps"

vcftools \
    --gzvcf "${FILE}.vcf.gz" \
    --max-missing 1 \
    --recode \
    --recode-INFO-all \
    --stdout | \
bgzip \
    > "${FILE}.no_missing.vcf.gz"

bcftools index \
    "${FILE}.no_missing.vcf.gz"
```

For LD pruning:

```bash
plink \
    --vcf "${FILE}.no_missing.vcf.gz" \
    --double-id \
    --make-bed \
    --out "${FILE}.no_missing"
```

Identify approximately independent SNPs:

```bash
plink \
    --bfile "${FILE}.no_missing" \
    --indep-pairwise 50 10 0.2 \
    --out "${FILE}.LD"
```

Create the pruned dataset:

```bash
plink \
    --bfile "${FILE}.no_missing" \
    --extract "${FILE}.LD.prune.in" \
    --make-bed \
    --out "${FILE}.LDpruned"
```

Export a VCF if required by the TreeMix conversion script:

```bash
plink \
    --bfile "${FILE}.LDpruned" \
    --recode vcf-iid bgz \
    --out "${FILE}.LDpruned"

bcftools index \
    "${FILE}.LDpruned.vcf.gz"
```

Prepare a cluster file:

```text
sample01    sample01    PopA
sample02    sample02    PopA
sample03    sample03    PopB
sample04    sample04    PopB
sample05    sample05    Outgroup
```

If the population is encoded in the sample name, the file can be generated automatically. For example:

```bash
bcftools query \
    -l "${FILE}.LDpruned.vcf.gz" | \
awk '{
    split($1, a, ".");
    print $1 "\t" $1 "\t" a[2]
}' \
    > population.clust
```

Only use this command when the sample naming convention actually supports it.

Convert the VCF to TreeMix format:

```bash
vcf2treemix.sh \
    "${FILE}.LDpruned.vcf.gz" \
    population.clust
```

The output is typically:

```text
population.treemix.frq.gz
```

Check:

```bash
zcat population.treemix.frq.gz | head
```

Set an outgroup:

```bash
ROOT="Outgroup"
```

Run the population tree without migration:

```bash
treemix \
    -i population.treemix.frq.gz \
    -m 0 \
    -root "$ROOT" \
    -o population.m0
```

Then test increasing numbers of migration edges:

```bash
for m in $(seq 0 5); do

    treemix \
        -i population.treemix.frq.gz \
        -m "$m" \
        -root "$ROOT" \
        -k 500 \
        -o "population.m${m}" \
        > "population.m${m}.log"

done
```

Here:

```text
-m    number of migration edges
-k    number of consecutive SNPs grouped into a block
-root population used to root the graph
```

For bootstrap analysis:

```bash
treemix \
    -i population.treemix.frq.gz \
    -m 2 \
    -root "$ROOT" \
    -bootstrap \
    -k 500 \
    -o population.m2.bootstrap
```

When the sampling design justifies it, TreeMix also provides:

```text
-noss
```

to disable the small-sample-size correction. This option should not be applied automatically.

Run multiple independent replicates for each migration number whenever possible because migration edges can vary among optimizations.

---

## 10. Plot TreeMix results

In R:

```r
library(RColorBrewer)
library(R.utils)

source("plotting_funcs.R")
```

Plot models with different numbers of migration edges:

```r
prefix <- "population"

par(mfrow = c(2, 3))

for (edge in 0:5) {

    plot_tree(
        cex = 0.8,
        paste0(
            prefix,
            ".m",
            edge
        )
    )

    title(
        paste(
            edge,
            "migration edges"
        )
    )
}
```

Prepare a population order file:

```text
PopA
PopB
PopC
PopD
Outgroup
```

Then plot residuals:

```r
for (edge in 0:5) {

    plot_resid(
        stem = paste0(
            prefix,
            ".m",
            edge
        ),
        pop_order = "population.list"
    )
}
```

For a population pair \(i,j\), the residual can be written conceptually as:

$$
R_{ij}
=
W_{ij}^{observed}
-
W_{ij}^{model}
$$

A large positive residual means that the two populations share more covariance than predicted by the fitted graph, suggesting that the current tree/migration model does not fully explain their relationship.

---

## 11. Combine TreeMix and f-statistics

TreeMix migration edges should generally be treated as hypotheses rather than final proof of introgression.

Suppose TreeMix suggests gene flow involving `PopA` and `PopC`.

Test:

```text
((PopA, PopB), PopC), Outgroup
```

with D:

```r
d_test <- d(
    W = "PopA",
    X = "PopB",
    Y = "PopC",
    Z = "Outgroup",
    data = snps
)

d_test
```

Then test f4:

```r
f4_test <- f4(
    W = "PopA",
    X = "PopB",
    Y = "PopC",
    Z = "Outgroup",
    data = snps
)

f4_test
```

If `PopA` is suspected to be admixed between populations related to `PopB` and `PopC`, test:

```r
f3_test <- f3(
    f2_blocks,
    "PopA",
    "PopB",
    "PopC"
)

f3_test
```

Then estimate ancestry proportions with qpAdm:

```r
qpadm_result <- qpadm(
    f2_blocks,
    left = c(
        "PopB",
        "PopC"
    ),
    right = c(
        "Outgroup1",
        "Outgroup2",
        "Outgroup3"
    ),
    target = "PopA"
)

qpadm_result$weights
qpadm_result$popdrop
```

If the genome-wide analyses support introgression and the aim is to identify introgressed genomic regions, continue with window-based statistics such as:

```text
windowed D
fd
fdM
Dxy
FST
local ancestry
```

For example, the \(f_d\) statistic is a normalized version of the ABBA-BABA numerator designed for genomic windows:

$$
f_d=
\frac{S(P_1,P_2,P_3,O)}
     {S(P_1,P_D,P_D,O)}
$$

where:

$$
S(P_1,P_2,P_3,O)
=
\sum_i
(p_{1i}-p_{2i})
(p_{3i}-p_{Oi})
$$

and \(P_D\) is chosen site-by-site from the population with the higher derived-allele frequency among the relevant donor/recipient populations according to the \(f_d\) formulation.

This allows genome-wide evidence of introgression to be followed by localization of candidate introgressed regions.

---

## 12. Important points for interpretation

For every D or f4 test, record the exact population order because changing the order changes the sign.

For D-statistics, always report:

```text
D(W,X;Y,Z)
D
SE
Z
number of informative SNPs
```
