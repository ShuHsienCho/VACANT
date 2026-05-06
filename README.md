# VACANT: Variant Annotation Clustering AssociatioN Test

**VACANT** is a robust R framework for rare variant association testing. It leverages **multi-dimensional annotation scores** to cluster variants into risk tiers using a **GMM-based Pareto Staircase** approach, followed by an aggregated association test (ACAT) or a joint multivariate Firth regression.

---

## Important: Data Requirements

1. **Rare Variants Only**: Filter your variants (e.g., MAF < 0.01) before creating input files.
2. **XPAT 2.0 Format**: Genotype matrix must be in gzipped XPAT 2.0 format (`.matrix.gz`).
3. **File Alignment**: PED, PCA, and score files must share the same sample IDs.

---

## Installation

```r
install.packages(c("remotes", "logistf", "mclust", "stringi",
                   "data.table", "dplyr", "optparse", "parallel"))
remotes::install_github("ShuHsienCho/VACANT")
```

---

## Input File Formats

### 1. Genotype Matrix (`--matrix`)

Gzipped XPAT 2.0 format (`.matrix.gz`):
- Cols 1-11: variant metadata; col 7 = gene name
- Header cols 12+: sample IDs
- Data col 12: one concatenated genotype string per variant
  (each character = one sample's allele count: 0 = ref, 1 = het, 2 = hom-alt, 5 = missing)

### 2. Phenotype File (`--ped`)

Standard PLINK PED format:
- Header: `#fid iid pid mid sex aff platform ignore`
- `aff`: 1 = control, 2 = case (recoded to 0/1 automatically)

### 3. Score File (`--score`)

Tab-delimited, no header:
- Cols 1-8: genomic coordinate metadata (must match matrix cols 1-8)
- Col 9+: numeric annotation scores (e.g., CASM, CADD, aPC)

### 4. Covariate File (`--covariates`)

No header; col 1 = IID, cols 2+ = numeric covariates (PCs, age, etc.)

---

## Usage

### R Library

```r
library(VACANT)

# Single-score analysis (e.g., CASM)
result <- vacant(
  matrix.file   = "XPAT.region.ATM_QC.matrix.gz",
  ped.file      = "UKB.BREAST.unrelated.ped",
  score.file    = "casm_spliceAI.txt",
  cov.file      = "UKB.BREAST.unrelated.pca",
  maf.threshold = 0.01,
  test          = "multi"
)

# Multi-score analysis (e.g., CADD + aPC)
result <- vacant(
  matrix.file   = "XPAT.region.ATM_QC.matrix.gz",
  ped.file      = "UKB.BREAST.unrelated.ped",
  score.file    = "casm_avg_spliceAI.txt",
  cov.file      = "UKB.BREAST.unrelated.pca",
  score.cols    = c("CADD", "aPC"),
  maf.threshold = 0.005,
  test          = "multi"
)

# Chromosome-level matrix with parallel execution
result <- vacant(
  matrix.file = "XPAT.chr17.matrix.gz",
  ped.file    = "UKB.BREAST.unrelated.ped",
  score.file  = "casm_spliceAI.txt",
  n.cores     = 8
)

# HPC: route temp files to scratch to avoid /tmp exhaustion
result <- vacant(
  matrix.file = "XPAT.chr17.matrix.gz",
  ped.file    = "UKB.BREAST.unrelated.ped",
  score.file  = "casm_spliceAI.txt",
  n.cores     = 8,
  tmpdir      = "/scratch/your_id/tmp"
)

# Clinical prediction
model_obj     <- readRDS("results/ATM.rds")
new_scores    <- read.csv("data/new_variants_scores.csv")
pred_clusters <- predict_vacant_cluster(model_obj$model, new_scores)
```

### Tier Summary (Post-Analysis Characterization)

```r
# Run single-gene analysis via vacant_core()
result <- vacant_core(geno, score, phenotype, covariates)

# Base summary: score distribution, carrier counts, effect estimates per tier
ts <- tier_summary(result, score, geno, phenotype,
                   variant.info = data.frame(chr, pos, ref, alt))

# Per-tier table
ts$tier
#   tier  n.variants  n.carriers.case  carrier.freq.case  OR  OR.lower  OR.upper  ...

# Per-variant table (ready for ggplot2)
ts$variant
#   chr  pos  ref  alt  tier  score  carrier.count  carrier.count.case  ...

# Optional enrichment with external annotations
ts <- enrich_vep(ts, vep.table)       # adds consequence breakdown
ts <- enrich_clinvar(ts, clinvar.table) # adds P/LP/VUS/LB/B counts

# Visualization
library(ggplot2)
ggplot(ts$variant, aes(x = factor(tier), y = score)) +
  geom_boxplot() +
  labs(x = "Tier", y = "Annotation Score")
```

### CLI (Bash)

#### Setup (one-time)

```bash
chmod +x $(Rscript -e "cat(system.file('bin', 'vacant', package='VACANT'))")
sudo cp $(Rscript -e "cat(system.file('bin', 'vacant', package='VACANT'))") /usr/local/bin/
```

#### Run

```bash
vacant \
  --matrix     "XPAT.region.ATM_QC.matrix.gz" \
  --ped        "UKB.BREAST.unrelated.ped" \
  --score      "casm_spliceAI.txt" \
  --covariates "UKB.BREAST.unrelated.pca" \
  --output     "results/ATM.csv" \
  --maf        0.01 \
  --test       multi
```

#### CLI Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--matrix` | Yes | | Gzipped XPAT 2.0 matrix file. |
| `--ped` | Yes | | PED phenotype file. |
| `--score` | Yes | | Annotation score file (no header). |
| `--output` | Yes | | Output CSV path (`.rds` also saved). |
| `--covariates` | No | NULL | Covariate/PCA file (no header). |
| `--score_cols` | No | NULL | Comma-separated score column names. |
| `--test` | No | `multi` | `multi` or `uni`. |
| `--weight` | No | `score` | ACAT weighting strategy: `score` (weight by cumulative tier center) or `equal`. |
| `--maf` | No | `0.01` | MAF threshold. |
| `--size_threshold` | No | `10` | Minimum variant count per tier (triggers merging). |
| `--gene_col` | No | `7` | Column index of gene name in matrix. |
| `--meta_ncols` | No | `11` | Number of metadata columns before genotype string. |
| `--n_cores` | No | `1` | CPU cores (set > 1 for chromosome-level matrices). |
| `--tmpdir` | No | `tempdir()` | Temp directory for `fread` intermediate files. |

---

## HPC Note: Parallel Jobs and /tmp Exhaustion

When submitting hundreds of concurrent LSF/SLURM jobs, `data.table::fread` with `cmd=` writes intermediate files to `/tmp`. On shared compute nodes, this can exhaust `/tmp`:

```
External command failed with exit code 2. This can happen when the disk is
full in the temporary directory ('/tmp/Rtmp...'). See ?fread for the tmpdir argument.
```

**Fix**: Redirect temp files to a scratch partition in your job script:

```bash
export TMPDIR=/scratch/your_id/tmp
mkdir -p $TMPDIR
```

---

## Output

Each run produces two files per gene:

- **`{output}.csv`**: Statistical results (p-values, per-tier betas, ORs, 95% CIs, cluster statistics).
- **`{output}_{gene}.rds`**: Pareto staircase model for use with `predict_vacant_cluster()`.

When using the R library interface, `vacant()` returns a `data.table` with one row per gene. List columns `.model` (Pareto staircase model) and `.cluster` (per-variant tier assignments) are included for downstream use with `predict_vacant_cluster()` and `tier_summary()`.

---

## Function Reference

| Function | Role |
|---|---|
| `vacant()` | Main entry point; accepts file paths |
| `vacant_core()` | Core engine; accepts R objects directly |
| `internal_prepare_inputs()` | Per-gene data preparation from XPAT matrix |
| `cluster_score()` | GMM-based clustering (E/EEI model) with Pareto anchors |
| `analyze_set()` | Firth regression (per-tier beta/OR/CI + gene-level penalized LRT) |
| `tier_summary()` | Post-analysis tier characterization (score, carrier, effect summary) |
| `enrich_vep()` | Add VEP functional consequence breakdown to tier summary |
| `enrich_clinvar()` | Add ClinVar classification counts to tier summary |
| `predict_vacant_cluster()` | Predict risk tier for new variants |
| `safe_logistf()` | Firth regression with auto-retry |
| `find_pareto_anchors_optimized()` | Pareto frontier identification |
| `acat_t()` / `acat_p()` | ACAT statistic and p-value conversion |
| `extract_sub()` | Genotype string extraction |
| ~~`perform_kmeans()`~~ | Deprecated. Retained for backward compatibility. |

---

## File Structure

```text
R/
├── vacant.R                   # Main entry point (file path interface)
├── vacant_core.R              # Core analysis engine (R object interface)
├── vacant_helpers.R           # Helper functions + internal_prepare_inputs
├── imports.R                  # Centralized @importFrom declarations
├── cluster_score.R            # GMM-based clustering (E/EEI)
├── analyze_set.R              # Firth regression + ACAT tests
├── tier_summary.R             # Post-analysis tier characterization + enrichment
└── predict_vacant_cluster.R   # Clinical prediction function
inst/
└── bin/
    └── vacant                 # CLI executable (Rscript with optparse)
inst/extdata/
├── geno.txt                   # Example genotype data
├── score.txt                  # Example annotation scores
├── pheno.txt                  # Example phenotype data
└── cov.txt                    # Example covariates
```
