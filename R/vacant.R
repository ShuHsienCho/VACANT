#' Run VACANT Rare Variant Association Analysis
#'
#' Main entry point for the VACANT package. Accepts file paths to an XPAT 2.0
#' genotype matrix, phenotype (PED), annotation scores, and optionally
#' covariates. Automatically detects all genes in the matrix, prepares inputs,
#' and runs the full VACANT pipeline (clustering + Firth association test) for
#' each gene. Supports parallel execution for multi-gene matrices.
#'
#' @param matrix.file Character. Path to gzipped XPAT 2.0 matrix file (.gz).
#'   May contain one or multiple genes.
#' @param ped.file Character. Path to PED file. Expected header:
#'   \code{#fid iid pid mid sex aff platform ignore}. Phenotype column
#'   (\code{aff}) uses PLINK 1/2 coding (1 = control, 2 = case).
#' @param score.file Character. Path to annotation score file (no header).
#'   Cols 1-8 must be genomic coordinate metadata matching the matrix;
#'   cols 9+ are numeric scores.
#' @param cov.file Character or NULL. Path to covariate file (no header).
#'   Col 1 = IID, cols 2+ = numeric covariates (PCs, age, etc.).
#'   (default: NULL).
#' @param score.cols Character vector or NULL. Names to assign to score
#'   columns (col 9+ of \code{score.file}). NULL uses a single column
#'   named \code{"score"}.
#' @param maf.threshold Numeric. Upper MAF bound for rare variant filtering
#'   (default: 0.01).
#' @param size.threshold Integer. Minimum cluster size for K-means
#'   (default: 10).
#' @param transform.method Character. Score transformation applied before
#'   clustering: \code{"none"} (default), \code{"raw_squared"},
#'   \code{"phred_to_chisq"}, \code{"log"}, or \code{"sigmoid"}.
#' @param test Character. \code{"multi"} (joint multivariate Firth; default)
#'   or \code{"uni"} (sequential univariate + ACAT).
#' @param acat.weight Character. ACAT weighting: \code{"score"} (default)
#'   or \code{"equal"}.
#' @param gene.col Integer. Column index of the gene name field in the matrix
#'   (default: 7, standard XPAT 2.0).
#' @param meta.ncols Integer. Number of metadata columns before the genotype
#'   string column (default: 11, standard XPAT 2.0).
#' @param n.cores Integer. Number of CPU cores for parallel gene processing.
#'   Use 1 (default) for single-gene matrices. Set higher for
#'   chromosome-level matrices with many genes (Linux/macOS only).
#' @param use.tabix Logical. If TRUE, converts the matrix to bgzip format and
#'   builds a tabix index for fast per-gene random access. Recommended for
#'   large region/chromosome-level matrices. If FALSE (default), uses direct
#'   zcat streaming, which is sufficient for small gene-level matrices.
#'
#' @return A \code{data.table} with one row per gene containing p-values,
#'   effect estimates, cluster statistics, and metadata. Returns NULL if
#'   all genes fail.
#'
#' @examples
#' \dontrun{
#' result <- vacant(
#'   matrix.file      = "XPAT.region.ATM_QC.matrix.gz",
#'   ped.file         = "UKB.BREAST.unrelated.ped",
#'   score.file       = "casm_avg_spliceAI.txt",
#'   cov.file         = "UKB.BREAST.unrelated.pca",
#'   score.cols       = c("CADD", "aPC"),
#'   maf.threshold    = 0.005,
#'   transform.method = "none",
#'   test             = "multi"
#' )
#' }
#'
#' @seealso \code{\link{vacant_core}}, \code{\link{predict_vacant_cluster}}
#' @export
vacant <- function(matrix.file,
                   ped.file,
                   score.file,
                   cov.file         = NULL,
                   score.cols       = NULL,
                   maf.threshold    = 0.01,
                   size.threshold   = 10L,
                   transform.method = c("none", "raw_squared",
                                        "phred_to_chisq", "log", "sigmoid"),
                   test             = c("multi", "uni"),
                   acat.weight      = c("score", "equal"),
                   gene.col         = 7L,
                   meta.ncols       = 11L,
                   n.cores          = 1L,
                   use.tabix        = FALSE) {

  transform.method <- match.arg(transform.method)
  test             <- match.arg(test)
  acat.weight      <- match.arg(acat.weight)

  # ---- 1. Validate file paths ----
  for (fp in c(matrix.file, ped.file, score.file)) {
    if (!file.exists(fp)) stop("File not found: ", fp)
  }
  if (!is.null(cov.file) && !file.exists(cov.file)) {
    stop("Covariate file not found: ", cov.file)
  }

  # ---- 2. Read header, PED, scores, covariates into memory ----
  # matrix is NOT read here - bgz index handles per-gene access.
  message(sprintf("[%s] Reading matrix header...",
                  format(Sys.time(), "%H:%M:%S")))

  col.names <- strsplit(
    system(paste("zcat", shQuote(matrix.file), "| head -1"), intern = TRUE),
    "\t"
  )[[1]]

  message(sprintf("[%s] %d samples in matrix.",
                  format(Sys.time(), "%H:%M:%S"),
                  length(col.names) - meta.ncols))

  message(sprintf("[%s] Reading ped...", format(Sys.time(), "%H:%M:%S")))
  ped.dt <- data.table::fread(ped.file, header = TRUE, data.table = FALSE)
  colnames(ped.dt) <- tolower(gsub("^#", "", colnames(ped.dt)))

  message(sprintf("[%s] Reading score file...", format(Sys.time(), "%H:%M:%S")))
  score.dt <- data.table::fread(score.file, header = FALSE, data.table = FALSE)

  cov.dt <- NULL
  if (!is.null(cov.file)) {
    message(sprintf("[%s] Reading covariates...", format(Sys.time(), "%H:%M:%S")))
    cov.dt <- data.table::fread(cov.file, header = FALSE, data.table = FALSE)
    colnames(cov.dt) <- c("iid", paste0("PC", seq_len(ncol(cov.dt) - 1L)))
  }

  # ---- 3. Prepare matrix access (bgz or gz) ----
  if (use.tabix) {
    mat.handle <- ensure_bgz_index(matrix.file, gene.col = gene.col)

    # Detect genes via tabix --list-chroms
    message(sprintf("[%s] Detecting genes...", format(Sys.time(), "%H:%M:%S")))
    raw.genes <- system(paste("tabix --list-chroms", shQuote(mat.handle)), intern = TRUE)
  } else {
    mat.handle <- matrix.file

    # Detect genes via zcat | awk (suitable for small/gene-level matrices)
    message(sprintf("[%s] Detecting genes...", format(Sys.time(), "%H:%M:%S")))
    raw.genes <- system(
      paste("zcat", shQuote(matrix.file),
            "| awk 'NR>1{print $", gene.col, "}' | sort -u"),
      intern = TRUE
    )
  }

  split.genes  <- trimws(unlist(strsplit(raw.genes, "[,;]")))
  unique.genes <- unique(
    split.genes[split.genes != "" & split.genes != "." & !is.na(split.genes)]
  )
  n.genes <- length(unique.genes)

  if (n.genes == 0L) stop("No valid genes detected in matrix column ", gene.col)

  message(sprintf("[%s] Detected %d gene(s).", format(Sys.time(), "%H:%M:%S"), n.genes))

  # ---- 4. Worker: per-gene matrix access ----
  process_gene <- function(target.gene) {

    result <- tryCatch({

      if (use.tabix) {
        # tabix: instant random access by gene name
        sub.mat <- data.table::fread(
          cmd        = paste("tabix", shQuote(mat.handle), shQuote(target.gene)),
          header     = FALSE,
          colClasses = "character"
        )
      } else {
        # zcat | awk: stream and filter (fast for small matrices)
        sub.mat <- data.table::fread(
          cmd        = paste0("zcat ", shQuote(mat.handle),
                              " | awk 'NR>1 && $", gene.col,
                              "==\"", target.gene, "\"'"),
          header     = FALSE,
          colClasses = "character"
        )
      }

      if (nrow(sub.mat) == 0L) {
        return(list(status = "fail", gene = target.gene,
                    reason = "no variants found in matrix"))
      }

      prepared <- internal_prepare_inputs(
        sub.matrix    = sub.mat,
        col.names     = col.names,
        ped.dt        = ped.dt,
        cov.dt        = cov.dt,
        score.dt      = score.dt,
        score.cols    = score.cols,
        maf.threshold = maf.threshold,
        meta.ncols    = meta.ncols
      )

      if (prepared$status != "success") {
        return(list(status = "fail", gene = target.gene,
                    reason = prepared$message))
      }

      vacant.obj <- tryCatch(
        vacant_core(
          geno             = prepared$geno,
          score            = prepared$score,
          phenotype        = prepared$phenotype,
          covariates       = prepared$covariates,
          test             = test,
          acat.weight      = acat.weight,
          size.threshold   = size.threshold,
          transform.method = transform.method
        ),
        error = function(e) {
          list(status = "fail", gene = target.gene,
               reason = paste("vacant_core() error:", e$message))
        }
      )

      # vacant_core returns NULL when there are no carriers
      if (is.null(vacant.obj)) {
        return(list(status = "fail", gene = target.gene,
                    reason = "vacant_core() returned NULL (no carriers)"))
      }
      # vacant_core returned a fail list (from inner tryCatch above)
      if (is.list(vacant.obj) && identical(vacant.obj$status, "fail")) {
        return(vacant.obj)
      }

      res.df            <- as.data.frame(t(vacant.obj$results))
      res.df$gene       <- target.gene
      res.df$n_variants <- length(prepared$geno)
      res.df$n_samples  <- nrow(prepared$phenotype)
      res.df$.model     <- list(vacant.obj$model)
      res.df

    }, error = function(e) {
      # Catch any unexpected error (e.g. from grepl, as.data.frame, etc.)
      list(status = "fail", gene = target.gene,
           reason = paste("unexpected error:", e$message))
    })

    result
  }

  # ---- 5. Execute (parallel or sequential) ----
  message(sprintf("[%s] Starting analysis...", format(Sys.time(), "%H:%M:%S")))

  if (n.cores > 1L) {
    results.list <- parallel::mclapply(unique.genes, process_gene,
                                       mc.cores = min(n.cores, n.genes))
  } else {
    results.list <- lapply(unique.genes, process_gene)
  }

  # ---- 6. Separate successes from failures and report ----
  # mclapply error objects (inherits "error") are caught by the outer tryCatch
  # inside process_gene, so they appear as fail lists, not raw error objects.
  is.success <- function(x) is.data.frame(x)
  is.fail    <- function(x) is.list(x) && identical(x$status, "fail")

  failed.list  <- Filter(is.fail,    results.list)
  results.list <- Filter(is.success, results.list)

  # Report all failures with reasons (visible in .e log even with mclapply)
  if (length(failed.list) > 0L) {
    fail.summary <- vapply(failed.list, function(x) {
      sprintf("  SKIP [%s]: %s", x$gene, x$reason)
    }, character(1L))
    message(paste(fail.summary, collapse = "\n"))
  }

  if (length(results.list) == 0L) {
    warning("All genes failed or returned no results.")
    return(NULL)
  }

  final.dt <- data.table::rbindlist(results.list, fill = TRUE)

  message(sprintf("[%s] Done. %d / %d genes returned results.",
                  format(Sys.time(), "%H:%M:%S"), nrow(final.dt), n.genes))
  final.dt
}
