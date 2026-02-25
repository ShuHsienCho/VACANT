#' Cluster annotation scores with Multi-Layer Pareto Staircase Boundaries
#'
#' Performs K-means clustering on standardized multi-dimensional scores and
#' identifies the Pareto frontier (Staircase) for each risk cluster.
#'
#' @param score Numeric matrix or vector. Rows are variants, columns are score types.
#' @param geno Character vector. Genotype strings for allele count weighting.
#' @param size.threshold Integer. Minimum cluster size threshold (default 10).
#' @param transform.method Character. "none", "raw_squared", "phred_to_chisq", or "log".
#' @return A list containing:
#'   \item{group.assignments}{Integer vector of cluster IDs.}
#'   \item{score.centers}{Matrix of cluster centers (Standardized & Sorted).}
#'   \item{ac.weights}{Numeric vector of scalar weights (Un-scaled sum) for ACAT.}
#'   \item{cluster.sizes}{Integer vector of cluster sizes.}
#'   \item{prediction.model}{List containing \code{anchors.list} for clinical prediction.}
#' @importFrom mclust Mclust mclustBIC mclustBICupdate hcRandomPairs
#' @importFrom stats kmeans sd scale
#' @importFrom stringi stri_extract_all_regex
#' @export
cluster_score <- function(score,
                          geno,
                          size.threshold = 10,
                          transform.method = c("none", "raw_squared", "phred_to_chisq", "log", "sigmoid")) {

  transform.method <- match.arg(transform.method)

  # ---- 1. Pre-processing (Matrix Conversion) ----
  if (is.vector(score)) {
    score.mat <- matrix(score, ncol = 1)
    colnames(score.mat) <- "Score"
  } else {
    score.mat <- as.matrix(score)
  }

  if (is.null(colnames(score.mat))) {
    colnames(score.mat) <- paste0("Score_", 1:ncol(score.mat))
  }

  if (any(is.na(score.mat))) stop("NA values found in score matrix.")

  # ---- 2. Transformation (NEW: Added raw_squared) ----
  shift.vals <- numeric(ncol(score.mat))

  if (transform.method == "raw_squared") {
    # [NEW] CADD/aPC Raw Score Logic:
    # 1. Truncate negatives (benign) to 0
    # 2. Square to expand dynamic range for K-means
    score.mat <- pmax(score.mat, 0)^2

  } else if (transform.method == "phred_to_chisq") {
    for (i in seq_len(ncol(score.mat))) {
      # PHRED scores are theoretically >= 0
      score.mat[, i] <- pmax(score.mat[, i], 0)

      # Replace exact zeros with a small floor before log transform
      # to avoid log(0) = -Inf, which maps to qchisq(Inf) = Inf.
      # This floor matches the handling in predict_vacant_cluster() so
      # training and prediction use the same numerical scale.
      curr.vals <- score.mat[, i]
      curr.vals[curr.vals == 0] <- 1e-6

      # Compute ln(P) numerically stably: ln(P) = -(PHRED/10) * ln(10)
      log.p.vals <- -(curr.vals / 10) * log(10)

      # Convert log-probability to chi-squared (df = 1) using log.p = TRUE
      # for numerical stability at extreme PHRED values
      score.mat[, i] <- qchisq(log.p.vals, df = 1, lower.tail = FALSE, log.p = TRUE)

      # Replace any residual Inf (can occur at very high PHRED) with
      # 110% of the largest finite value
      if (any(is.infinite(score.mat[, i]))) {
        max.val <- max(score.mat[!is.infinite(score.mat[, i]), i], na.rm = TRUE)
        score.mat[is.infinite(score.mat[, i]), i] <- max.val * 1.1
      }
    }
  } else if (transform.method == "log") {
    # Shift columns with negative values to non-negative before log1p.
    # shift.vals records the per-column shift so predict_vacant_cluster()
    # can apply the same offset to new data.
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
  # If a column is constant (sd = 0), set sd = 1 to avoid division by zero
  col.sds[col.sds == 0] <- 1

  score.scaled <- scale(score.mat, center = col.means, scale = col.sds)

  # ---- 4. Expand Scores (Allele Count Weighting) ----
  # Parse the genotype string to derive per-variant allele counts.
  # Heterozygous (AC=1) contributes one row; homozygous alt (AC=2) contributes
  # two rows. This weights the clustering toward higher-AC variants.
  ac <- vapply(stringi::stri_extract_all_regex(geno, "\\d"),
               function(chars) sum(as.integer(chars)), integer(1))

  if (sum(ac) == 0) return(NULL)

  # Expand the score matrix according to allele count
  idx.expanded   <- rep(seq_len(nrow(score.scaled)), ac)
  score.expanded <- score.scaled[idx.expanded, , drop = FALSE]

  # ---- 5. Clustering (GMM + K-means) ----
  n.unique.points  <- nrow(unique(score.expanded))
  max.clust.search <- min(20, n.unique.points - 1)

  if (max.clust.search < 1) {
    # Only one unique point; degenerate case
    km <- stats::kmeans(score.expanded, centers = score.expanded[1, , drop = FALSE])
  } else {
    bic.all <- NULL
    for (i in 1:max.clust.search) {
      suppressMessages({
        bic.all <- mclust::mclustBICupdate(
          bic.all,
          mclust::mclustBIC(score.expanded, verbose = FALSE, G = i,
                            initialization = list(hcPairs = mclust::hcRandomPairs(score.expanded)))
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
  sizes   <- km$size
  centers <- km$centers
  cluster.indices <- km$cluster

  while (any(sizes < size.threshold) && length(sizes) > 1) {
    i.min <- which.min(sizes)
    # Calculate Euclidean distance to all other centers
    dists <- apply(centers, 1, function(x) sum((x - centers[i.min,])^2))
    dists[i.min] <- Inf
    merge.with <- which.min(dists)

    # Merge logic
    new.size <- sizes[i.min] + sizes[merge.with]
    new.center <- (sizes[i.min] * centers[i.min,] + sizes[merge.with] * centers[merge.with,]) / new.size

    sizes[merge.with]   <- new.size
    centers[merge.with,] <- new.center

    # Remove the merged cluster
    sizes   <- sizes[-i.min]
    centers <- centers[-i.min, , drop = FALSE]

    # Re-assign indices
    old.cluster.vec <- cluster.indices
    cluster.indices[old.cluster.vec == i.min] <- merge.with

    # Compact indices to 1..K
    current.ids <- sort(unique(cluster.indices))
    cluster.indices <- match(cluster.indices, current.ids)
  }

  # ---- 7. Sort Clusters & Calculate Weights ----
  # Sort by magnitude of centers (risk)
  center.vals <- rowSums(centers)
  ord <- order(center.vals)

  sorted.centers <- centers[ord, , drop = FALSE]
  rownames(sorted.centers) <- seq_len(nrow(sorted.centers)) # Reset names

  sorted.sizes <- sizes[ord]

  # Remap assignments based on sort order
  map.vec <- integer(nrow(centers))
  map.vec[ord] <- 1:nrow(centers)
  sorted.assignments.expanded <- map.vec[cluster.indices]

  # Contract back to variant level (taking the first assignment for each variant)
  final.assignments <- as.integer(
    tapply(sorted.assignments.expanded, idx.expanded, function(x) x[1])
  )

  # Handle NAs if any variant was lost (unlikely with this logic)
  final.assignments[is.na(final.assignments)] <- 1

  # 2. Calculate Weights based on "Original" (Un-scaled) Magnitude
  # This ensures weights reflect the raw signal strength (e.g. Chi-sq magnitude)
  unscaled.centers <- t(t(sorted.centers) * col.sds + col.means)

  # Sum across dimensions to get a scalar burden
  raw.burden <- rowSums(unscaled.centers)

  # Ensure Weights are Positive for ACAT
  min.burden <- min(raw.burden)
  if (min.burden <= 0) {
    final.weights <- raw.burden - min.burden + 0.1
  } else {
    final.weights <- raw.burden
  }

  # ---- 8. Identify "Layered" Pareto Anchors ----
  # Use the pre-transform (but pre-scale) scores so anchor thresholds are
  # expressed in the same units as the raw input to predict_vacant_cluster().
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

  # ---- Output ----
  list(
    group.assignments = final.assignments,
    score.centers     = sorted.centers,
    ac.weights        = final.weights,
    cluster.sizes     = sorted.sizes,
    prediction.model  = list(
      type = "layered_pareto",
      anchors.list = anchors.list,
      K = K,
      transform.method = transform.method,
      shift.vals = shift.vals,
      scale.mean = col.means,
      scale.sd = col.sds
    )
  )
}
