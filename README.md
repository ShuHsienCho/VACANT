# VACANT: Variant Annotation Clustering AssociatioN Test

**VACANT** is a robust R framework for rare variant association testing. It leverages **multi-dimensional annotation scores** to cluster variants into risk tiers using a conservative **Pareto Staircase** approach, followed by an aggregated association test (ACAT) or a joint multivariate Firth regression.

---

## ⚠️ Important: Data Requirements

Before preparing your input files, please adhere to these critical rules:

1.  **Rare Variants Only**: You must filter your variants (e.g., MAF < 0.01) **BEFORE** creating the input files.
2.  **Row Alignment**: All input files (`geno`, `score`, `pheno`, `cov`) must be **row-aligned**. The $i$-th row in every file must correspond to the same sample/variant.
3.  **Genotype Format**: Genotypes must be stored as **strings** (e.g., "010010..."), never as numbers.

---

## ⚡️ Quick Start with Example Data

VACANT comes with built-in example datasets in `inst/extdata`.

### Run in R

```r
library(VACANT)
library(data.table)

# 1. Locate example files
geno_file  <- system.file("extdata", "geno.txt", package = "VACANT")
score_file <- system.file("extdata", "score.txt", package = "VACANT")
pheno_file <- system.file("extdata", "pheno.txt", package = "VACANT")
cov_file   <- system.file("extdata", "cov.txt", package = "VACANT")

# 2. Load Data
# [CRITICAL] Read genotype as "character" to preserve string format
geno_df  <- data.table::fread(geno_file, header = TRUE, colClasses = "character", data.table = FALSE)
score_df <- data.table::fread(score_file, header = TRUE, data.table = FALSE)
pheno_df <- data.table::fread(pheno_file, header = TRUE, data.table = FALSE)
cov_df   <- data.table::fread(cov_file, header = FALSE, data.table = FALSE)

# 3. Run Analysis
result <- vacant(
  geno         = geno_df[[1]],      # Extract the genotype column
  score        = as.matrix(score_df), 
  phenotype    = data.frame(phenotype = pheno_df[[1]]), 
  covariates   = cov_df,
  test         = "multi"
)

print(result$results)

```

---

## 📄 Input File Formats

VACANT accepts standard text files (CSV/TXT).

### 1. Genotype File (`--geno`)

* **Format**: A file with a header and a column of genotype strings.
* **Crucial**: Must be read as `character` in R.

**Example (`geno.txt`):**

```text
genotype
01000000100000...
00000000000000...

```

### 2. Score File (`--score`)

* **Format**: A matrix of numeric scores with header names.

**Example (`score.txt`):**

```text
CADD,REVEL
25.4,0.8
10.2,0.1

```

### 3. Phenotype File (`--pheno`)

* **Format**: A file containing phenotype (0/1).

### 4. Covariate File (`--cov`)

* **Format**: A file containing covariates (e.g. PC1, Age).

---

## 🛠 Installation & Setup

VACANT can be used as a **Command Line Tool** (like SAIGE/PLINK) or as a standard **R Library**.

### 1. Install the R Package

```r
# In R console
install.packages(c("devtools", "optparse", "data.table", "stringi", "logistf", "mclust", "dplyr"))
devtools::install_local(".") 

```

### 2. Set up the Command Line Tool (CLI)

```bash
cd inst/bin/
chmod +x vacant
sudo cp vacant /usr/local/bin/

```

---

## 🚀 CLI Usage (Terminal)

### Mode 1: Analysis (Training)

```bash
vacant \
  --geno data/geno.txt \
  --score data/score.txt \
  --pheno data/pheno.txt \
  --cov data/cov.txt \
  --score_cols "CADD,REVEL" \
  --test multi \
  --output results/TP53_result

```

| Argument | Description |
| --- | --- |
| `--geno` | **[Required]** Path to genotype string file. |
| `--score` | **[Required]** Path to annotation scores. |
| `--pheno` | **[Required]** Path to phenotype file. |
| `--cov` | **[Optional]** Path to covariate file. |
| `--score_cols` | **[Required]** Score column names (must match header). |
| `--output` | **[Required]** Output prefix. |

### Mode 2: Clinical Prediction

Use a previously trained model (`.rds`) to predict risk tiers for new variants.

```bash
vacant \
  --input new_variants_scores.csv \
  --model results/TP53_result.rds \
  --output results/predictions.csv

```

*(Note: `--input` should contain the raw scores for new variants. The CLI handles log-transformation automatically.)*

---

## 💻 R Library Usage

How to load your data in R.

### 1. Analysis (Training)

```r
library(VACANT)
library(data.table)

# 1. Load Data
# [CRITICAL] Always force colClasses = "character" for genotypes
geno_df <- fread("data/geno.txt", header = TRUE, colClasses = "character", data.table = FALSE)
geno_vec <- geno_df[[1]] 

# Load scores & pheno
score_mat <- as.matrix(fread("data/score.txt", header = TRUE))
pheno_df  <- fread("data/pheno.txt", header = TRUE, data.table = FALSE)
cov_df    <- fread("data/cov.txt", header = TRUE, data.table = FALSE)

# 2. Run Analysis
result_obj <- vacant(
  geno         = geno_vec,
  score        = score_mat,
  phenotype    = pheno_df,
  covariates   = cov_df,
  test         = "multi", 
  acat.weight  = "score"
)

# View Results
print(result_obj$results)

```

### 2. Clinical Prediction

You can predict risk clusters for new variants using the model object.

```r
library(VACANT)

# Step A: Load a trained model
model_obj <- readRDS("results/TP53_result.rds")

# Step B: Prepare new data (Raw scores)
# Ensure columns match the training data
new_scores <- read.csv("data/new_variants.csv")

# Step C: Predict Risk Clusters
# Returns integer IDs (1 = Lowest Risk, K = Highest Risk)
pred_clusters <- predict_vacant_cluster(model_obj$model, new_scores)

# View
print(pred_clusters)

```

---

## 📂 File Structure

```text
R/
├── vacant.R                  # Main analysis engine
├── cluster_score.R           # Log-transform & Pareto Clustering logic
├── analyze_set.R             # Statistical tests (Firth / ACAT)
├── predict_vacant_cluster.R  # Prediction function
inst/bin/
└── vacant                    # CLI Executable Script
inst/extdata/
└── geno.txt                  # Example Data
└── score.txt                 # Example Data
└── pheno.txt                 # Example Data
└── cov.txt                   # Example Data
