# ==============================================================================
# vacant_helpers.R
# Helper functions for the VACANT package.
# All @importFrom declarations are centralized in imports.R.
# ==============================================================================


# ---- ACAT Helpers ------------------------------------------------------------

#' Calculate ACAT Test Statistic
#'
#' Computes the aggregated Cauchy association test statistic from a vector
#' of p-values and corresponding weights.
#'
#' @param p Numeric vector of p-values (each in (0, 1)).
#' @param w Numeric vector of positive weights of the same length as \code{p}.
#' @return Numeric scalar: the ACAT test statistic.
#' @export
acat_t <- function(p, w) {
  sum(w * tan((0.5 - p) * pi))
}

#' Convert ACAT Statistic to a P-value
#'
#' Given an ACAT test statistic and its weights, returns the combined p-value.
#'
#' @param t Numeric scalar: ACAT test statistic from \code{acat_t()}.
#' @param w Numeric vector of weights used to compute \code{t}.
#' @return Numeric scalar: the combined p-value.
#' @export
acat_p <- function(t, w) {
  0.5 - atan(t / sum(w)) / pi
}


# ---- Clustering Helper -------------------------------------------------------

#' Perform K-means Clustering with Retry
#'
#' Attempts k-means clustering on \code{data} using GMM-derived initial
#' centers. If an error occurs, waits one second and retries recursively.
#'
#' @param data Numeric matrix or data frame.
#' @param mod A Mclust model object containing a \code{classification} vector.
#' @return An object of class \code{\link[stats]{kmeans}}.
#' @export
perform_kmeans <- function(data, mod) {
  tryCatch({
    k <- length(unique(mod$classification))
    stats::kmeans(data, centers = k, nstart = 25)
  }, error = function(e) {
    message("Error in perform_kmeans: ", e$message, "; retrying in 1 second...")
    Sys.sleep(1)
    perform_kmeans(data, mod)
  })
}


# ---- Genotype Helper ---------------------------------------------------------

#' Extract Per-Sample Alleles from Genotype Strings
#'
#' Given a vector of XPAT-format genotype strings (one per variant) and a
#' vector of sample positions, extracts a single concatenated string of
#' allele counts for those samples.
#'
#' @param x Character vector of genotype strings.
#' @param idx Integer vector. Sample positions to extract.
#' @return Character scalar of extracted alleles (one character per sample).
#' @export
extract_sub <- function(x, idx) {
  stringi::stri_c(stringi::stri_sub(x, from = idx, to = idx), collapse = "")
}


# ---- Regression Helper -------------------------------------------------------

#' Safe Firth Logistic Regression with Retry
#'
#' Wraps \code{logistf::logistf} to perform Firth-corrected logistic
#' regression, automatically retrying if the Fisher information matrix
#' becomes singular.
#'
#' @param formula A \code{formula} specifying the logistic regression model.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family Character or function. Model family (default: \code{"binomial"}).
#' @param control A control object from \code{logistf::logistf.control()}.
#' @param max_retries Integer. Maximum retry attempts (default: 5).
#' @return An object of class \code{logistf}.
#' @export
safe_logistf <- function(formula, data, family = "binomial",
                         control, max_retries = 5) {
  attempt <- 1
  repeat {
    warn_msg <- NULL
    fit <- withCallingHandlers(
      try(logistf::logistf(formula, data = data,
                           family = family, control = control),
          silent = TRUE),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("Determinant of Fisher information matrix was numerically 0",
                  msg)) {
          warn_msg <<- msg
          invokeRestart("muffleWarning")
        }
      }
    )
    if (is.null(warn_msg) && !inherits(fit, "try-error")) return(fit)
    message(sprintf("Warning in logistf (attempt %d/%d): %s. Retrying...",
                    attempt, max_retries, warn_msg))
    attempt <- attempt + 1
    if (attempt > max_retries) {
      stop("logistf failed after ", max_retries,
           " attempts due to singular Fisher information.")
    }
  }
}


# ---- Pareto Helper -----------------------------------------------------------

#' Find Pareto Anchors (Lower-Left Frontier)
#'
#' Implements a Sort-and-Scan algorithm to identify non-dominated points
#' (Pareto frontier) in O(N log N) for 2D or O(N x k) for higher dimensions.
#'
#' @param points Numeric matrix of points (e.g., annotation scores of variants
#'   in a risk cluster).
#' @return A matrix of anchor points defining the lower-left boundary,
#'   subset from the input \code{points}.
#' @export
find_pareto_anchors_optimized <- function(points) {
  if (is.null(points) || nrow(points) == 0) return(NULL)
  if (nrow(points) == 1) return(points)

  pts        <- unique(points)
  ord        <- do.call(order, as.data.frame(pts))
  pts.sorted <- pts[ord, , drop = FALSE]
  anchors    <- matrix(NA_real_, nrow = nrow(pts.sorted), ncol = ncol(pts.sorted))
  n.anchors  <- 0

  if (ncol(pts) == 2) {
    min.y.so.far <- Inf
    for (i in seq_len(nrow(pts.sorted))) {
      curr.y <- pts.sorted[i, 2]
      if (curr.y < min.y.so.far) {
        n.anchors            <- n.anchors + 1
        anchors[n.anchors, ] <- pts.sorted[i, ]
        min.y.so.far         <- curr.y
      }
    }
  } else {
    for (i in seq_len(nrow(pts.sorted))) {
      pt           <- pts.sorted[i, ]
      is.dominated <- FALSE
      if (n.anchors > 0) {
        curr.anchors <- anchors[1:n.anchors, , drop = FALSE]
        diffs        <- sweep(curr.anchors, 2, pt, "-")
        if (any(rowSums(diffs <= 1e-9) == ncol(pts))) is.dominated <- TRUE
      }
      if (!is.dominated) {
        n.anchors            <- n.anchors + 1
        anchors[n.anchors, ] <- pt
      }
    }
  }
  anchors[1:n.anchors, , drop = FALSE]
}


# ---- Data Preparation --------------------------------------------------------

#' Prepare VACANT Inputs from an XPAT Matrix Block
#'
#' Takes a single-gene subset of an XPAT 2.0 genotype matrix (already loaded
#' in memory) and prepares all inputs required by \code{vacant_core()}.
#' Handles sample alignment, phenotype recoding, covariate subsetting,
#' score joining, genotype string extraction, and MAF filtering.
#'
#' @param sub.matrix Data frame. Rows = variants for one gene; column layout
#'   follows XPAT 2.0: cols 1-\code{meta.ncols} are metadata, col
#'   \code{meta.ncols + 1} is the concatenated genotype string.
#' @param col.names Character vector. Full column names from the matrix header
#'   (length = \code{meta.ncols} + n_samples). Used to recover sample IDs.
#' @param ped.dt Data frame. Phenotype table with lowercase column names
#'   (\code{iid} and \code{aff} required). PLINK 1/2 coding (1 = control,
#'   2 = case) is recoded to 0/1 automatically.
#' @param cov.dt Data frame or NULL. Covariate table; col 1 = IID,
#'   cols 2+ = numeric covariates (PCs, age, etc.).
#' @param score.dt Data frame. Annotation score table; cols 1-8 are genomic
#'   coordinate metadata (matching the matrix), cols 9+ are scores.
#' @param score.cols Character vector or NULL. Names to assign to score columns
#'   (col 9+ of \code{score.dt}). NULL uses a single column named
#'   \code{"score"}.
#' @param maf.threshold Numeric. Upper MAF bound; variants with MAF >=
#'   \code{maf.threshold} or MAF == 0 are excluded (default: 0.01).
#' @param meta.ncols Integer. Number of metadata columns before the genotype
#'   string column in the XPAT matrix (default: 11).
#'
#' @return A named list. On success: \code{status = "success"} plus
#'   \code{geno} (character vector), \code{score} (data frame),
#'   \code{phenotype} (data frame), \code{covariates} (data frame or NULL).
#'   On failure: \code{status = "fail"} plus \code{message}.
#' @export
internal_prepare_inputs <- function(sub.matrix,
                                    col.names,
                                    ped.dt,
                                    cov.dt        = NULL,
                                    score.dt,
                                    score.cols    = NULL,
                                    maf.threshold = 0.01,
                                    meta.ncols    = 11L) {

  # ---- 1. Sample alignment ----
  matrix.samples <- col.names[(meta.ncols + 1L):length(col.names)]

  if (!"iid" %in% colnames(ped.dt)) {
    return(list(status = "fail", message = "'iid' column not found in ped"))
  }

  common.samples <- intersect(matrix.samples, ped.dt[["iid"]])
  if (length(common.samples) == 0L) {
    return(list(status = "fail",
                message = "No sample overlap between matrix and ped"))
  }

  sample.idx <- match(common.samples, matrix.samples)

  # ---- 2. Phenotype ----
  pheno.rows <- match(common.samples, ped.dt[["iid"]])
  aff.vals   <- ped.dt[["aff"]][pheno.rows]

  if (all(na.omit(aff.vals) %in% c(1L, 2L))) {
    aff.vals <- aff.vals - 1L
  }

  pheno.final    <- data.frame(phenotype = aff.vals, row.names = common.samples)
  pheno.final    <- stats::na.omit(pheno.final)
  common.samples <- rownames(pheno.final)
  sample.idx     <- match(common.samples, matrix.samples)

  if (length(common.samples) == 0L) {
    return(list(status = "fail",
                message = "No samples remaining after phenotype NA removal"))
  }

  # ---- 3. Covariates ----
  cov.final <- NULL
  if (!is.null(cov.dt)) {
    cov.iid.col <- colnames(cov.dt)[1]
    cov.rows    <- match(common.samples, cov.dt[[cov.iid.col]])
    cov.sub     <- cov.dt[cov.rows, -1L, drop = FALSE]
    rownames(cov.sub) <- common.samples
    colnames(cov.sub) <- paste0("PC", seq_len(ncol(cov.sub)))
    cov.final <- cov.sub
  }

  # ---- 4. Score alignment ----
  # Convert only the small 8-col slice to data.frame for dplyr::left_join.
  # Avoids converting the full matrix (potentially thousands of variants).
  var.info   <- as.data.frame(sub.matrix[, 1:8, drop = FALSE])
  colnames(var.info) <- col.names[1:8]

  score.work <- score.dt
  colnames(score.work)[1:8] <- col.names[1:8]
  # var.info is character (fread colClasses="character"); match types for join
  score.work[, 1:8] <- lapply(score.work[, 1:8], as.character)

  if (is.null(score.cols)) {
    colnames(score.work)[9] <- "score"
    target.cols <- "score"
  } else {
    n.sc <- length(score.cols)
    colnames(score.work)[9:(9 + n.sc - 1L)] <- score.cols
    target.cols <- score.cols
  }

  joined    <- suppressMessages(
    dplyr::left_join(var.info, score.work, by = col.names[1:8])
  )
  score.mat <- as.data.frame(joined)[, target.cols, drop = FALSE]

  keep.score <- which(rowSums(is.na(score.mat)) == 0L)
  if (length(keep.score) == 0L) {
    return(list(status = "fail",
                message = "All variants missing scores after join"))
  }

  sub.matrix <- sub.matrix[keep.score, ]
  score.mat  <- score.mat[keep.score, , drop = FALSE]

  # ---- 5. Genotype extraction ----
  geno.strings <- sub.matrix[[meta.ncols + 1L]]

  geno.raw <- vapply(
    geno.strings,
    extract_sub,
    FUN.VALUE = character(1L),
    idx = sample.idx
  )
  geno.raw <- chartr("5", "0", geno.raw)

  # ---- 6. MAF filtering ----
  n.chroms <- 2L * length(common.samples)
  ac       <- stringi::stri_count_fixed(geno.raw, "1") +
    2L * stringi::stri_count_fixed(geno.raw, "2")
  maf.val  <- ac / n.chroms
  maf.pass <- which(maf.val > 0 & maf.val < maf.threshold)

  if (length(maf.pass) == 0L) {
    return(list(status = "fail", message = "No variants passed MAF filter"))
  }

  list(
    status     = "success",
    geno       = geno.raw[maf.pass],
    score      = score.mat[maf.pass, , drop = FALSE],
    phenotype  = pheno.final,
    covariates = cov.final
  )
}

#' Ensure bgzip + tabix index exists for a matrix file
#'
#' Converts a .gz matrix file to bgzip format and builds a tabix index
#' keyed on the gene name column. The .bgz and .tbi files are written
#' next to the original .gz file. This is a one-time setup; subsequent
#' calls return immediately if the index already exists.
#'
#' @param gz.file Character. Path to the gzipped matrix file.
#' @param gene.col Integer. Column index of the gene name field (default 7).
#' @return Character. Path to the .bgz file.
#' @export
ensure_bgz_index <- function(gz.file, gene.col = 7L) {

  bgz.file <- sub("\\.gz$", ".bgz", gz.file)
  tbi.file <- paste0(bgz.file, ".tbi")

  if (!file.exists(bgz.file) || !file.exists(tbi.file)) {
    message(sprintf(
      "[%s] Building bgzip+tabix index (one-time setup, may take a few minutes)...",
      format(Sys.time(), "%H:%M:%S")
    ))

    ret <- system(paste("zcat", shQuote(gz.file), "| bgzip >", shQuote(bgz.file)))
    if (ret != 0L) stop("bgzip conversion failed for: ", gz.file)

    ret <- system(paste(
      "tabix -s", gene.col, "-b 2 -e 3 -S 1", shQuote(bgz.file)
    ))
    if (ret != 0L) stop("tabix indexing failed for: ", bgz.file)

    message(sprintf("[%s] Index ready: %s", format(Sys.time(), "%H:%M:%S"), bgz.file))
  } else {
    message(sprintf("[%s] Using existing bgz index: %s",
                    format(Sys.time(), "%H:%M:%S"), bgz.file))
  }

  bgz.file
}
