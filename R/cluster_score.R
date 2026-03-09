#' Cluster annotation scores with Multi-Layer Pareto Staircase Boundaries
#'
#' Performs K-means clustering on standardized multi-dimensional scores and
#' identifies the Pareto frontier (Staircase) for each risk cluster.
#'
#' @param score Numeric matrix or vector. Rows are variants, columns are score types.
#' @param geno Character vector. Genotype strings for allele count weighting.
#' @param size.threshold Integer. Minimum cluster size threshold (default 10).
#' @param transform.method Character. "none", "raw_squared", "phred_to_chisq", or "log".
#' @param cond.number.warn Numeric. If the condition number of the feature
#'   covariance matrix exceeds this threshold, a warning is printed and
#'   near-collinear columns are reported. High condition numbers indicate
#'   near-singular covariance, which destabilizes GMM BIC selection and
#'   K-means initialization. Default: 100.
#'
#' @return A list containing:
#'   \item{group.assignments}{Integer vector of cluster IDs.}
#'   \item{score.centers}{Matrix of cluster centers (Standardized & Sorted).}
#'   \item{ac.weights}{Numeric vector of scalar weights (Un-scaled sum) for ACAT.}
#'   \item{cluster.sizes}{Integer vector of cluster sizes.}
#'   \item{prediction.model}{List containing \code{anchors.list} for clinical prediction.}
#'   \item{cond.number}{Condition number of the standardized feature covariance.}
#' @importFrom mclust Mclust mclustBIC mclustBICupdate hcRandomPairs
#' @importFrom stats kmeans sd scale
#' @importFrom stringi stri_extract_all_regex
#' @export
cluster_score <- function(score,
                          geno,
                          size.threshold   = 10,
                          transform.method = c("none", "raw_squared",
                                               "phred_to_chisq", "log", "sigmoid"),
                          cond.number.warn = 100) {

  transform.method <- match.arg(transform.method)

  # ---- 1. Pre-processing (Matrix Conversion) ----
  if (is.vector(score)) {
    score.mat <- matrix(score, ncol = 1)
    colnames(score.mat) <- "Score"
  } else {
    score.mat <- as.matrix(score)
  }

  if (is.null(colnames(score.mat))) {
    colnames(score.mat) <- paste0("Score_", seq_len(ncol(score.mat)))
  }

  if (any(is.na(score.mat))) stop("NA values found in score matrix.")

  # ---- 2. Transformation ----
  shift.vals <- numeric(ncol(score.mat))

  if (transform.method == "raw_squared") {
    score.mat <- pmax(score.mat, 0)^2

  } else if (transform.method == "phred_to_chisq") {
    for (i in seq_len(ncol(score.mat))) {
      score.mat[, i] <- pmax(score.mat[, i], 0)
      curr.vals <- score.mat[, i]
      curr.vals[curr.vals == 0] <- 1e-6
      log.p.vals <- -(curr.vals / 10) * log(10)
      score.mat[, i] <- qchisq(log.p.vals, df = 1, lower.tail = FALSE, log.p = TRUE)
      if (any(is.infinite(score.mat[, i]))) {
        max.val <- max(score.mat[!is.infinite(score.mat[, i]), i], na.rm = TRUE)
        score.mat[is.infinite(score.mat[, i]), i] <- max.val * 1.1
      }
    }
  } else if (transform.method == "log") {
    for (i in seq_len(ncol(score.mat))) {
      col.min <- min(score.mat[, i])
      if (col.min < 0) {
        shift.vals[i] <- abs(col.min)
        score.mat[, i] <- score.mat[, i] + shift.vals[i]
      }
    }
    score.mat <- log1p(score.mat)
  } else if (transform.method == "sigmoid") {
    score.mat <- 1 / (1 + exp(-score.mat))
  }

  # ---- 3. Scaling (Standardization) ----
  col.means <- colMeans(score.mat)
  col.sds   <- apply(score.mat, 2, sd)
  col.sds[col.sds == 0] <- 1

  score.scaled <- scale(score.mat, center = col.means, scale = col.sds)

  # ---- 3c. Jitter for near-zero variance columns ----
  # When all indels in a gene receive an identical imputed score on one
  # dimension (e.g. apc_protein_function: all frameshifts get the same
  # gene-local median), that column has variance = 0 in score.scaled.
  # A zero-variance column creates a point mass in the feature space,
  # which causes the GMM covariance matrix to be singular and can induce
  # spurious clusters along other dimensions.
  # Fix: add a small N(0, jitter.sd^2) perturbation to any column whose
  # empirical SD in the scaled matrix is below jitter.threshold.
  # The perturbation is negligible relative to the inter-variant spread
  # on other dimensions but restores numerical rank to the covariance matrix.
  jitter.threshold <- 1e-6
  jitter.sd        <- 1e-4
  for (j in seq_len(ncol(score.scaled))) {
    if (sd(score.scaled[, j]) < jitter.threshold) {
      warning(sprintf(
        paste0("cluster_score: column '%s' has near-zero variance (sd < %.0e) ",
               "after standardization. This typically occurs when all variants ",
               "share an identical imputed value (e.g. gene-local median for ",
               "apc_protein_function). Adding N(0, %.0e) jitter to restore ",
               "covariance matrix rank."),
        colnames(score.scaled)[j], jitter.threshold, jitter.sd^2
      ))
      set.seed(42L + j)
      score.scaled[, j] <- score.scaled[, j] +
        stats::rnorm(nrow(score.scaled), mean = 0, sd = jitter.sd)
    }
  }

  # ---- 3b. Condition number check (Challenge C defense) ----
  # A near-singular covariance matrix (high condition number) indicates that
  # two or more feature columns are nearly collinear. This destabilizes GMM
  # covariance estimation and makes BIC-based K selection unreliable.
  # The check is performed on the complete-case standardized matrix BEFORE
  # allele-count expansion, so it reflects the user's input feature set.
  # Recommended pre-screening: exclude columns with pairwise |r| > 0.90.
  cond.number <- NA_real_
  if (ncol(score.scaled) > 1) {
    cov.mat     <- crossprod(score.scaled) / max(nrow(score.scaled) - 1, 1)
    evals       <- tryCatch(eigen(cov.mat, symmetric = TRUE, only.values = TRUE)$values,
                            error = function(e) NULL)
    if (!is.null(evals) && min(abs(evals)) > 0) {
      cond.number <- max(abs(evals)) / min(abs(evals))
    }

    if (!is.na(cond.number) && cond.number > cond.number.warn) {
      # Identify the most collinear column pair
      cor.mat  <- stats::cor(score.scaled)
      diag(cor.mat) <- 0
      max.abs.r <- max(abs(cor.mat), na.rm = TRUE)
      idx       <- which(abs(cor.mat) == max.abs.r, arr.ind = TRUE)[1, ]
      warning(sprintf(
        paste0("cluster_score: feature covariance condition number = %.1f ",
               "(threshold = %.0f). The most correlated pair is columns ",
               "'%s' and '%s' (|r| = %.4f). ",
               "Consider excluding one of these columns before clustering to ",
               "stabilize GMM BIC selection. ",
               "If using apc_*_raw inputs, verify that no pair with |r| > 0.90 ",
               "was included (see inter-block correlation in Phase 1 output)."),
        cond.number, cond.number.warn,
        colnames(score.scaled)[idx[1]],
        colnames(score.scaled)[idx[2]],
        max.abs.r
      ))
    }
  }

  # ---- 4. Expand Scores (Allele Count Weighting) ----
  ac <- vapply(stringi::stri_extract_all_regex(geno, "\\d"),
               function(chars) sum(as.integer(chars)), integer(1))

  if (sum(ac) == 0) return(NULL)

  idx.expanded   <- rep(seq_len(nrow(score.scaled)), ac)
  score.expanded <- score.scaled[idx.expanded, , drop = FALSE]

  # ---- 5. Clustering (GMM + K-means) ----
  n.unique.points  <- nrow(unique(score.expanded))
  max.clust.search <- min(20, n.unique.points - 1)

  if (max.clust.search < 1) {
    km <- stats::kmeans(score.expanded,
                        centers = score.expanded[1, , drop = FALSE])
  } else {
    bic.all <- NULL
    for (i in seq_len(max.clust.search)) {
      suppressMessages({
        bic.all <- mclust::mclustBICupdate(
          bic.all,
          mclust::mclustBIC(score.expanded, verbose = FALSE, G = i,
                            initialization = list(
                              hcPairs = mclust::hcRandomPairs(score.expanded)))
        )
      })
    }
    mod <- mclust::Mclust(score.expanded, x = bic.all, verbose = FALSE)

    if (is.null(mod) || max(mod$classification) == 1) {
      km <- stats::kmeans(score.expanded, centers = 1)
    } else {
      km <- perform_kmeans(score.expanded, mod)
    }
  }

  # ---- 6. Merge Small Clusters ----
  sizes           <- km$size
  centers         <- km$centers
  cluster.indices <- km$cluster

  while (any(sizes < size.threshold) && length(sizes) > 1) {
    i.min <- which.min(sizes)
    dists <- apply(centers, 1, function(x) sum((x - centers[i.min, ])^2))
    dists[i.min] <- Inf
    merge.with <- which.min(dists)

    new.size   <- sizes[i.min] + sizes[merge.with]
    new.center <- (sizes[i.min] * centers[i.min, ] +
                     sizes[merge.with] * centers[merge.with, ]) / new.size

    sizes[merge.with]    <- new.size
    centers[merge.with, ] <- new.center
    sizes   <- sizes[-i.min]
    centers <- centers[-i.min, , drop = FALSE]

    old.cluster.vec <- cluster.indices
    cluster.indices[old.cluster.vec == i.min] <- merge.with
    current.ids     <- sort(unique(cluster.indices))
    cluster.indices <- match(cluster.indices, current.ids)
  }

  # ---- 7. Sort Clusters & Calculate Weights ----
  center.vals    <- rowSums(centers)
  ord            <- order(center.vals)
  sorted.centers <- centers[ord, , drop = FALSE]
  rownames(sorted.centers) <- seq_len(nrow(sorted.centers))
  sorted.sizes   <- sizes[ord]

  map.vec <- integer(nrow(centers))
  map.vec[ord] <- seq_len(nrow(centers))
  sorted.assignments.expanded <- map.vec[cluster.indices]

  final.assignments <- as.integer(
    tapply(sorted.assignments.expanded, idx.expanded, function(x) x[1])
  )
  final.assignments[is.na(final.assignments)] <- 1

  unscaled.centers <- t(t(sorted.centers) * col.sds + col.means)
  raw.burden       <- rowSums(unscaled.centers)
  min.burden       <- min(raw.burden)
  final.weights    <- if (min.burden <= 0) {raw.burden - min.burden + 0.1 }  else {raw.burden}

  # ---- 8. Pareto Anchors ----
  score.expanded.raw <- score.mat[idx.expanded, , drop = FALSE]
  K            <- nrow(sorted.centers)
  anchors.list <- vector("list", K)

  for (k in seq_len(K)) {
    k.indices <- which(sorted.assignments.expanded == k)
    if (length(k.indices) > 0) {
      k.points  <- score.expanded.raw[k.indices, , drop = FALSE]
      anchors.k <- find_pareto_anchors_optimized(k.points)
      if (!is.null(anchors.k)) {
        anchors.list[[k]] <- anchors.k[order(anchors.k[, 1]), , drop = FALSE]
      }
    }
  }

  list(
    group.assignments = final.assignments,
    score.centers     = sorted.centers,
    ac.weights        = final.weights,
    cluster.sizes     = sorted.sizes,
    cond.number       = cond.number,
    prediction.model  = list(
      type             = "layered_pareto",
      anchors.list     = anchors.list,
      K                = K,
      transform.method = transform.method,
      shift.vals       = shift.vals,
      scale.mean       = col.means,
      scale.sd         = col.sds
    )
  )
}
