# VACANT: Variant Annotation Clustering AssociatioN Test

VACANT is an R package for rare variant association testing in biobank-scale genomic data. It clusters rare variants into risk tiers using multi-dimensional annotation scores and a conservative **Pareto Staircase** boundary, then performs aggregated association testing via Firth logistic regression and ACAT.

---

## Overview

Standard burden tests assume all rare variants in a gene contribute equally to disease risk. VACANT relaxes this assumption by:

1. Using annotation scores (e.g., CADD, SpliceAI) to cluster variants into ordered risk tiers via GMM-initialized K-means
2. Identifying the Pareto frontier (Staircase boundary) of each risk cluster to enable interpretable clinical prediction
3. Testing each cluster against the binary phenotype using Firth logistic regression (with correction for separation), then combining across clusters via ACAT

The trained Pareto model can also be used to classify new variants into risk clusters without re-running the full association pipeline.

---

## Installation

```r
install.packages(c(
  "remotes", "logistf", "mclust", "stringi",
  "data.table", "dplyr", "optparse", "parallel"
))
remotes::install_github("ShuHsienCho/VACANT")
```

---

## Input File Formats

### Genotype Matrix (`--matrix`)

Gzipped XPAT 2.0 format (`.matrix.gz`):

- **Cols 1-11**: variant metadata; col 7 = gene name (default)
- **Header cols 12+**: sample IDs
- **Data col 12**: one concatenated genotype string per variant, where each character encodes one sample's allele count
  - `0` = homozygous reference
  - `1` = heterozygous
  - `2` = homozygous alternate
  - `5` = missing (recoded to `0` automatically)

### Phenotype File (`--ped`)

Standard PLINK PED format with header:

```
#fid  iid  pid  mid  sex  aff  platform  ignore
```

The `aff` column uses PLINK coding (`1` = control, `2` = case), which is automatically recoded to `0`/`1`.

### Score File (`--score`)

Tab-delimited, no header:

- **Cols 1-8**: genomic coordinate metadata (must match matrix cols 1-8 exactly for join)
- **Col 9+**: numeric annotation scores (e.g., CADD, aPC, SpliceAI)

### Covariate File (`--covariates`, optional)

No header; col 1 = IID, cols 2+ = numeric covariates (principal components, age, etc.)

---

## Quick Start

### R Interface

```r
library(VACANT)

# --- Analysis: file-path interface ---
result <- vacant(
  matrix.file      = "XPAT.region.ATM_QC.matrix.gz",
  ped.file         = "UKB.BREAST.unrelated.ped",
  score.file       = "casm_avg_spliceAI.txt",
  cov.file         = "UKB.BREAST.unrelated.pca",
  score.cols       = c("CADD", "aPC"),
  maf.threshold    = 0.005,
  transform.method = "none",
  test             = "multi"
)

# --- Multi-gene matrix with parallel execution (Linux/macOS) ---
result <- vacant(
  matrix.file = "XPAT.chr17.matrix.gz",
  ped.file    = "UKB.BREAST.unrelated.ped",
  score.file  = "casm_avg_spliceAI.txt",
  n.cores     = 8
)

# --- Advanced: R-object interface (data already in memory) ---
result.obj <- vacant_core(
  geno             = prepared$geno,
  score            = prepared$score,
  phenotype        = prepared$phenotype,
  covariates       = prepared$covariates,
  test             = "multi",
  transform.method = "none"
)

# --- Clinical prediction: classify new variants into risk clusters ---
model.obj     <- readRDS("results/ATM_ATM.rds")
new.scores    <- read.csv("new_variants_scores.csv")
pred.clusters <- predict_vacant_cluster(model.obj$model, new.scores)
```

### CLI (Bash)

#### One-time setup

```bash
# Make the CLI executable
chmod +x $(Rscript -e "cat(system.file('bin', 'vacant', package='VACANT'))")

# Optional: copy to PATH so 'vacant' works from anywhere
sudo cp $(Rscript -e "cat(system.file('bin', 'vacant', package='VACANT'))") /usr/local/bin/
```

#### Analysis mode

```bash
vacant \
  --matrix     "XPAT.region.ATM_QC.matrix.gz" \
  --ped        "UKB.BREAST.unrelated.ped" \
  --score      "casm_avg_spliceAI.txt" \
  --covariates "UKB.BREAST.unrelated.pca" \
  --output     "results/ATM.csv" \
  --score_cols "CADD,aPC" \
  --maf        0.005 \
  --test       multi \
  --transform  none
```

#### Prediction mode

```bash
vacant \
  --model      "results/ATM_ATM.rds" \
  --input      "new_variants_scores.csv" \
  --score_cols "CADD,aPC" \
  --output     "predictions.csv"
```

#### Without PATH setup

```bash
Rscript $(Rscript -e "cat(system.file('bin', 'vacant', package='VACANT'))") \
  --matrix "XPAT.region.ATM_QC.matrix.gz" \
  --ped    "UKB.BREAST.unrelated.ped" \
  --score  "casm_avg_spliceAI.txt" \
  --output "results/ATM.csv"
```

---

## Parameters

### `vacant()` -- Main R Function

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `matrix.file` | character | required | Path to gzipped XPAT 2.0 matrix (`.gz`) |
| `ped.file` | character | required | Path to PED phenotype file |
| `score.file` | character | required | Path to annotation score file (no header) |
| `cov.file` | character | `NULL` | Path to covariate/PCA file (no header) |
| `score.cols` | character vector | `NULL` | Names for score columns (col 9+ of score file); `NULL` uses `"score"` |
| `maf.threshold` | numeric | `0.01` | Upper MAF bound for rare variant filtering |
| `size.threshold` | integer | `10` | Minimum variant count per K-means cluster |
| `transform.method` | character | `"none"` | Score transformation before clustering; see below |
| `test` | character | `"multi"` | `"multi"` (joint multivariate Firth) or `"uni"` (sequential univariate + ACAT) |
| `acat.weight` | character | `"score"` | ACAT weights: `"score"` (cluster score magnitude) or `"equal"` |
| `gene.col` | integer | `7` | Column index of gene name in matrix |
| `meta.ncols` | integer | `11` | Number of metadata columns before the genotype string column |
| `n.cores` | integer | `1` | CPU cores for parallel gene processing (Linux/macOS only) |

**Tabix auto-detection**: if a pre-built `.bgz` + `.tbi` index exists alongside the `.gz` file, `vacant()` automatically uses tabix for fast per-gene random access. Otherwise it falls back to `zcat` streaming. Recommended for chromosome-level or region-level matrices.

### Score Transformations (`transform.method`)

| Value | Description | Recommended for |
|-------|-------------|-----------------|
| `"none"` | No transformation | Pre-normalized scores |
| `"raw_squared"` | Truncate negatives to 0, then square | CADD raw scores, aPC |
| `"phred_to_chisq"` | Convert PHRED-scaled scores to chi-squared (df=1) | CADD PHRED, REVEL |
| `"log"` | Shift to non-negative then apply `log1p` | Skewed positive scores |
| `"sigmoid"` | Logistic transform | Log-odds scale scores |

### CLI Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--matrix` | Yes (analysis) | | Gzipped XPAT 2.0 matrix |
| `--ped` | Yes (analysis) | | PED phenotype file |
| `--score` | Yes (analysis) | | Annotation score file |
| `--output` | Yes | | Output file path |
| `--covariates` | No | NULL | Covariate/PCA file |
| `--score_cols` | No | NULL | Comma-separated score column names |
| `--test` | No | `multi` | `multi` or `uni` |
| `--weight` | No | `score` | `score` or `equal` |
| `--maf` | No | `0.01` | MAF threshold |
| `--size_threshold` | No | `10` | Minimum cluster size |
| `--transform` | No | `none` | `none`, `raw_squared`, `phred_to_chisq`, `log`, `sigmoid` |
| `--gene_col` | No | `7` | Column index of gene name |
| `--meta_ncols` | No | `11` | Number of metadata columns |
| `--n_cores` | No | `1` | CPU cores for parallel execution |
| `--no_rds` | No | `FALSE` | Skip saving per-gene `.rds` model files. Recommended for null/inflation simulations |
| `--rds_dir` | No | NULL | Directory to save `.rds` files separately from the CSV. `NULL` saves alongside the CSV |
| `--model` | Yes (prediction) | | Path to trained `.rds` model |
| `--input` | Yes (prediction) | | CSV/TXT file with scores for new variants |

---

## Output

Each run produces:

- **`{output}.csv`**: one row per gene with p-values, effect estimates, and cluster statistics

| Column | Description |
|--------|-------------|
| `gene` | Gene name |
| `p.firth` | One-sided Firth p-value for aggregated burden (all variants) |
| `clusters` | Number of clusters identified |
| `tx.score` | Weighted average cluster score (L2 magnitude) |
| `beta.1`, `p.1`, ... | Per-cluster effect estimate and one-sided p-value (`multi` mode) |
| `acat.multi` | ACAT-combined p-value across clusters (`multi` mode) |
| `p.lrt.multi` | Profile likelihood ratio test p-value (`multi` mode) |
| `p.score.uni1`, ... | Per-cluster cumulative burden p-value (`uni` mode) |
| `acat.uni` | ACAT-combined p-value (`uni` mode) |
| `n_variants` | Number of variants passing MAF filter |
| `n_samples` | Number of samples with non-missing phenotype |

- **`{stem}_{gene}.rds`**: Pareto staircase model per gene for use with `predict_vacant_cluster()`. Saved alongside the CSV by default, or in a separate directory if `--rds_dir` is specified. Can be suppressed entirely with `--no_rds TRUE`.

The `.rds` object contains:

```r
list(
  results = list(...),   # same statistics as the CSV row
  model   = list(...),   # Pareto staircase model for prediction
  gene    = "GENE_NAME"
)
```

### RDS output control (CLI only)

R users receive the model directly in the returned list and manage persistence themselves. The CLI provides two flags to control automatic `.rds` saving:

| Flag | Use case |
|------|----------|
| `--no_rds TRUE` | Null/inflation simulations where Pareto models are never needed. Avoids generating thousands of unused files. |
| `--rds_dir <path>` | Keep model files in a dedicated `models/` subdirectory, separate from result CSVs. |

---

## HPC Usage (LSF / MD Anderson)

`generate_lsf.py` generates LSF job scripts for three scenarios:

```bash
# Real-data analysis: one job per genomic region
python3 generate_lsf.py whole --cancer breast

# Inflation simulation: one job per region
python3 generate_lsf.py inflation --size_threshold 10

# Power simulation: one job per gene, from a curated gene list
python3 generate_lsf.py power --cancer breast --reps 200
```

---

## Function Reference

| Function | Description |
|----------|-------------|
| `vacant()` | Main entry point; reads files, detects genes, runs full pipeline |
| `vacant_core()` | Core engine; accepts R objects directly (for advanced use) |
| `internal_prepare_inputs()` | Per-gene data preparation from XPAT matrix block |
| `cluster_score()` | GMM + K-means clustering with Pareto staircase boundary |
| `analyze_set()` | Firth regression + ACAT statistical tests |
| `predict_vacant_cluster()` | Classify new variants into risk clusters via saved Pareto model |
| `safe_logistf()` | Firth logistic regression with auto-retry on singular Fisher matrix |
| `find_pareto_anchors_optimized()` | O(N log N) Pareto frontier identification |
| `perform_kmeans()` | K-means with GMM-derived initial centers and error retry |
| `acat_t()` / `acat_p()` | ACAT statistic computation and p-value conversion |
| `extract_sub()` | Per-sample allele extraction from XPAT genotype strings |

---

## File Structure

```
R/
├── vacant.R                   # Main entry point (file-path interface)
├── vacant_core.R              # Core analysis engine (R-object interface)
├── vacant_helpers.R           # Helper functions + internal_prepare_inputs
├── import.R                   # Centralized @importFrom declarations
├── cluster_score.R            # GMM + K-means clustering + Pareto anchors
├── analyze_set.R              # Firth regression + ACAT tests
└── predict_vacant_cluster.R   # Clinical prediction via Pareto model
inst/
└── bin/
    └── vacant                 # CLI executable (Rscript + optparse)
inst/extdata/
├── geno.txt                   # Example genotype data
├── score.txt                  # Example annotation scores
├── pheno.txt                  # Example phenotype data
└── cov.txt                    # Example covariates
```

---

## Dependencies

| Package | Role |
|---------|------|
| `logistf` | Firth-corrected logistic regression |
| `mclust` | Gaussian mixture model for cluster initialization |
| `stringi` | Fast string operations on genotype strings |
| `data.table` | File I/O and result aggregation |
| `dplyr` | Score table alignment via left join |
| `parallel` | Multi-core gene processing |
| `optparse` | CLI argument parsing |