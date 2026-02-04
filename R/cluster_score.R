#' Cluster annotation scores with Multi-Layer Pareto Staircase Boundaries
#'
#' Performs K-means clustering on standardized multi-dimensional scores and
#' identifies the Pareto frontier (Staircase) for each risk cluster.
#'
#' @param score Numeric matrix or vector. Rows are variants, columns are score types.
#' @param geno Character vector. Genotype strings for allele count weighting.
#' @param size.threshold Integer. Minimum cluster size threshold (default 10).
#' @param transform.method Character. "log" (default, applies log1p) or "none".
#' @return A list containing:
#'   \item{group.assignments}{Integer vector of cluster IDs.}
#'   \item{score.centers}{Matrix of cluster centers (Standardized & Sorted).}
#'   \item{ac.weights}{Numeric vector of scalar weights (Un-scaled sum) for ACAT.}
#'   \item{cluster.sizes}{Integer vector of cluster sizes.}
#'   \item{prediction.model}{List containing \code{anchors.list} for clinical prediction.}
#' @import mclust stats stringi
#' @export
cluster_score <- function(score,
                          geno,
                          size.threshold = 10,
                          transform.method = c("log", "none")) {

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

  # ---- 2. Transformation (Fixing Skewness) ----
  shift.vals <- numeric(ncol(score.mat))

  if (transform.method == "log") {
    for (i in seq_len(ncol(score.mat))) {
      col.min <- min(score.mat[, i])
      if (col.min < 0) {
        shift.vals[i] <- abs(col.min)
        score.mat[, i] <- score.mat[, i] + shift.vals[i]
      }
    }
    score.mat <- log1p(score.mat)
  }

  # ---- 3. Scaling (Standardization) ----
  col.means <- colMeans(score.mat)
  col.sds   <- apply(score.mat, 2, sd)
  col.sds[col.sds == 0] <- 1

  score.scaled <- scale(score.mat, center = col.means, scale = col.sds)

  # ---- 4. Expand Scores (Allele Count Weighting) ----
  ac <- vapply(stringi::stri_extract_all_regex(geno, "\\d"),
               function(chars) sum(as.integer(chars)), integer(1))

  if (sum(ac) == 0) return(NULL)

  idx.expanded <- rep(seq_len(nrow(score.scaled)), ac)
  score.expanded <- score.scaled[idx.expanded, , drop = FALSE]

  # ---- 5. Clustering (GMM + K-means) ----
  n.unique.points <- nrow(unique(score.expanded))
  max.clust.search <- min(20, n.unique.points - 1)

  if (max.clust.search < 1) {
    km <- stats::kmeans(score.expanded, centers = score.expanded[1,,drop=FALSE])
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
      # Use helper function defined in vacant_helpers.R
      km <- perform_kmeans(score.expanded, mod)
    }
  }

  # ---- 6. Merge Small Clusters ----
  sizes   <- km$size
  centers <- km$centers
  cluster.indices <- km$cluster

  while (any(sizes < size.threshold) && length(sizes) > 1) {
    i.min <- which.min(sizes)
    dists <- apply(centers, 1, function(x) sum((x - centers[i.min,])^2))
    dists[i.min] <- Inf
    merge.with <- which.min(dists)

    new.size <- sizes[i.min] + sizes[merge.with]
    new.center <- (sizes[i.min] * centers[i.min,] + sizes[merge.with] * centers[merge.with,]) / new.size

    sizes[merge.with]   <- new.size
    centers[merge.with,] <- new.center
    sizes   <- sizes[-i.min]
    centers <- centers[-i.min, , drop = FALSE]

    old.cluster.vec <- cluster.indices
    cluster.indices[old.cluster.vec == i.min] <- merge.with
    current.ids <- sort(unique(cluster.indices))
    cluster.indices <- match(cluster.indices, current.ids)
  }

  # ---- 7. Sort Clusters & Calculate Weights ----
  # [Corrected Logic]
  # 1. Sort by Algebraic Value (Risk): Low (Benign) -> High (Pathogenic)
  center.vals <- rowSums(centers)
  ord <- order(center.vals)

  sorted.centers <- centers[ord, , drop = FALSE]
  rownames(sorted.centers) <- seq_len(nrow(sorted.centers)) # Reset names

  sorted.sizes <- sizes[ord]

  # Remap assignments
  map.vec <- integer(nrow(centers))
  map.vec[ord] <- 1:nrow(centers)
  sorted.assignments.expanded <- map.vec[cluster.indices]

  final.assignments <- as.integer(
    tapply(sorted.assignments.expanded, idx.expanded, function(x) x[1])
  )

  # 2. Calculate Weights based on "Original" (Un-scaled) Magnitude
  # Un-scale the centers: x_raw = x_z * sd + mean
  # We use the transformed space (Log space), which is proportional to risk.
  unscaled.centers <- t(t(sorted.centers) * col.sds + col.means)

  # Sum across dimensions to get a scalar burden
  raw.burden <- rowSums(unscaled.centers)

  # Ensure Weights are Positive for ACAT
  min.burden <- min(raw.burden)
  if (min.burden <= 0) {
    # Shift so the lowest risk cluster has a small positive weight (e.g. 0.1)
    final.weights <- raw.burden - min.burden + 0.1
  } else {
    final.weights <- raw.burden
  }

  # ---- 8. Identify "Layered" Pareto Anchors ----
  score.expanded.raw <- score.mat[idx.expanded, , drop = FALSE]

  K <- nrow(sorted.centers)
  anchors.list <- vector("list", K)

  if (K > 1) {
    for (k in 2:K) {
      k.indices <- which(sorted.assignments.expanded == k)
      if (length(k.indices) > 0) {
        k.points <- score.expanded.raw[k.indices, , drop = FALSE]
        anchors.k <- find_pareto_anchors_optimized(k.points)
        anchors.list[[k]] <- anchors.k[order(anchors.k[,1]), , drop = FALSE]
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
