# ==============================================================================
# cluster_score.R
# GMM-based variant tier identification for the VACANT package.
#
# Pipeline:
#   1. Pre-processing      : matrix conversion
#   2. Transformation      : none / raw_squared / phred_to_chisq / log / sigmoid
#   3. Scaling             : standardize columns to zero mean, unit variance
#   4. GMM (E/EEI)         : BIC-based K selection + posterior classification
#                            p >= 2: EEI (diagonal, equal volume and shape)
#                            p  = 1: E   (univariate, equal variance)
#   4b. Empty component drop: remove components with 0 assigned variants
#   5. Merge small tiers   : pooled-covariance Mahalanobis distance on centers
#   6. Pareto anchors      : identify Pareto frontier per tier
#
# Design notes:
#   - AC-weighting removed (2026-05). GMM fits directly on annotation score
#     vectors without allele-count row replication. Rationale: AC expansion
#     biases GMM toward higher-frequency variants and inflates the effective
#     sample size used by BIC, distorting K selection. The GMM should cluster
#     variants by annotation profile alone; allele frequency information enters
#     only through the association test (Firth regression), not through tier
#     identification.
#
#   - E/EEI covariance constraint (2026-05). V/VVI allows a single component's
#     variance to expand until it absorbs the entire dataset, causing K=1
#     collapse (empirically confirmed on BRCA1, CHEK2, PALB2 without AC
#     expansion). E/EEI's equal-variance constraint acts as regularization:
#     all components share the same variance, so no single component can
#     dominate by inflating its variance. This forces BIC to use multiple
#     components when the score distribution is genuinely multimodal.
#
#   - Empty component drop (step 4b). mclust may fit G components but assign
#     zero variants to some via MAP classification. These are dropped before
#     merging so that effective K = number of non-empty components. This
#     prevents phantom tiers from inflating the LRT degrees of freedom.
#
#   - K-means step removed (2026-03). GMM mod$classification is the primary
#     assignment. Rationale: GMM BIC selects K under the chosen covariance
#     model; using the same model's posterior classification is internally
#     consistent. K-means re-assignment used a weaker spherical assumption
#     and introduced a model-selection / assignment inconsistency.
#
#   - Merge uses pooled within-cluster Mahalanobis distance on cluster centers,
#     consistent with the GMM geometry used for primary assignment.
#
# Note on perform_kmeans() in vacant_helpers.R:
#   perform_kmeans() is no longer called by cluster_score(). It is retained
#   in vacant_helpers.R for backward compatibility but is deprecated.
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helper: pooled within-cluster covariance
# ------------------------------------------------------------------------------

compute_pooled_cov <- function(data, classification) {
  d    <- ncol(data)
  n    <- nrow(data)
  k    <- max(classification)
  sigma.pooled <- matrix(0, d, d)

  for (ck in seq_len(k)) {
    idx <- which(classification == ck)
    nk  <- length(idx)
    wk  <- nk / n
    if (nk > d) {
      sigma.pooled <- sigma.pooled + wk * cov(data[idx, , drop = FALSE])
    } else {
      sigma.pooled <- sigma.pooled + wk * diag(d)
    }
  }
  sigma.pooled
}


# ------------------------------------------------------------------------------
# Internal helper: safe inverse of a symmetric positive-definite matrix
# ------------------------------------------------------------------------------

safe_solve <- function(sigma) {
  tryCatch(
    chol2inv(chol(sigma)),
    error = function(e) {
      sv  <- svd(sigma)
      tol <- max(dim(sigma)) * .Machine$double.eps * sv$d[1]
      sv$v %*% diag(ifelse(sv$d > tol, 1 / sv$d, 0)) %*% t(sv$u)
    }
  )
}


# ------------------------------------------------------------------------------
# Internal helper: Mahalanobis distance between two vectors given sigma.inv
# ------------------------------------------------------------------------------

mahal_dist <- function(c1, c2, sigma.inv) {
  diff <- as.numeric(c1) - as.numeric(c2)
  sqrt(max(0, as.numeric(t(diff) %*% sigma.inv %*% diff)))
}


# ------------------------------------------------------------------------------
# Main clustering function
# ------------------------------------------------------------------------------

#' Cluster annotation scores with Multi-Layer Pareto Staircase Boundaries
#'
#' Performs GMM-based clustering on standardized multi-dimensional annotation
#' scores and identifies the Pareto frontier (Staircase) for each risk tier.
#'
#' The pipeline uses a Gaussian mixture model (GMM) with E covariance
#' parameterization (equal variance; EEI for multivariate input) for both
#' K selection via BIC and primary tier assignment via posterior classification.
#' Empty components (fitted by mclust but assigned zero variants via MAP) are
#' dropped before merging.  Small tiers below \code{size.threshold} variants
#' are then merged to the nearest tier by pooled within-cluster Mahalanobis
#' distance.
#'
#' @param score Numeric matrix or vector. Rows are variants, columns are score
#'   dimensions.
#' @param size.threshold Integer. Minimum number of variants per tier before
#'   merging is triggered (default 10).
#' @param transform.method Character. Pre-clustering score transformation:
#'   \code{"none"}, \code{"raw_squared"}, \code{"phred_to_chisq"},
#'   \code{"log"}, or \code{"sigmoid"}.
#'
#' @return A list with elements:
#'   \item{group.assignments}{Integer vector of tier IDs, one per variant.}
#'   \item{score.centers}{Matrix of tier centers on the standardized scale,
#'     sorted by ascending mean score.}
#'   \item{ac.weights}{Numeric vector of ACAT weights per tier.}
#'   \item{cluster.sizes}{Integer vector of tier sizes (variant counts).}
#'   \item{prediction.model}{List for clinical prediction via
#'     \code{predict_vacant_cluster()}.}
#'
#' @importFrom mclust Mclust mclustBIC hcRandomPairs
#' @importFrom stats cov sd scale
#' @export
cluster_score <- function(score,
                          size.threshold = 10,
                          transform.method = c("none", "raw_squared",
                                               "phred_to_chisq", "log",
                                               "sigmoid")) {

  transform.method <- match.arg(transform.method)

  # ---- 1. Pre-processing (matrix conversion) --------------------------------
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

  n.variants <- nrow(score.mat)
  p.dim      <- ncol(score.mat)

  # ---- 2. Transformation ----------------------------------------------------
  shift.vals <- numeric(p.dim)

  if (transform.method == "raw_squared") {
    score.mat <- pmax(score.mat, 0)^2

  } else if (transform.method == "phred_to_chisq") {
    for (i in seq_len(p.dim)) {
      score.mat[, i] <- pmax(score.mat[, i], 0)
      curr.vals <- score.mat[, i]
      curr.vals[curr.vals == 0] <- 1e-6
      log.p.vals <- -(curr.vals / 10) * log(10)
      score.mat[, i] <- qchisq(log.p.vals, df = 1,
                               lower.tail = FALSE, log.p = TRUE)
      if (any(is.infinite(score.mat[, i]))) {
        max.val <- max(score.mat[!is.infinite(score.mat[, i]), i], na.rm = TRUE)
        score.mat[is.infinite(score.mat[, i]), i] <- max.val * 1.1
      }
    }

  } else if (transform.method == "log") {
    for (i in seq_len(p.dim)) {
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

  # ---- 3. Scaling (standardize to zero mean, unit variance) ----------------
  col.means <- colMeans(score.mat)
  col.sds   <- apply(score.mat, 2, sd)
  col.sds[col.sds == 0] <- 1
  score.scaled <- scale(score.mat, center = col.means, scale = col.sds)

  # ---- 4. GMM: BIC-based K selection + posterior tier assignment -----------
  # Covariance model: E for p = 1 (equal variance, univariate),
  #                   EEI for p >= 2 (diagonal, equal volume and shape).
  n.unique.points  <- nrow(unique(score.scaled))
  max.clust.search <- min(20L, n.unique.points - 1L)

  gmm.model.name <- if (p.dim == 1L) "E" else "EEI"

  if (max.clust.search < 1L) {
    tier.indices <- rep(1L, n.variants)
    gmm.mod      <- NULL
  } else {
    suppressWarnings(suppressMessages({
      bic.obj <- mclust::mclustBIC(
        score.scaled, G = 1:max.clust.search, verbose = FALSE,
        modelNames = gmm.model.name,
        initialization = list(hcPairs = mclust::hcRandomPairs(score.scaled))
      )
    }))
    suppressWarnings(suppressMessages({
      gmm.mod <- mclust::Mclust(score.scaled, x = bic.obj, verbose = FALSE)
    }))

    if (is.null(gmm.mod) || max(gmm.mod$classification) == 1L) {
      tier.indices <- rep(1L, n.variants)
      gmm.mod      <- NULL
    } else {
      tier.indices <- as.integer(gmm.mod$classification)
    }
  }

  # ---- 4b. Drop empty components -------------------------------------------
  # mclust may fit G components but MAP assigns zero variants to some.
  # Remap to contiguous 1:K_effective before merging.
  occupied     <- sort(unique(tier.indices))
  if (length(occupied) < max(tier.indices)) {
    remap        <- integer(max(tier.indices))
    remap[occupied] <- seq_along(occupied)
    tier.indices <- remap[tier.indices]
  }

  # ---- 5. Merge small tiers (pooled-covariance Mahalanobis) ---------------
  k.current <- max(tier.indices)

  compute_centers <- function(data, assignment, k) {
    raw <- sapply(seq_len(k), function(ck) {
      idx <- which(assignment == ck)
      if (length(idx) == 1L) as.numeric(data[idx, , drop = FALSE])
      else as.numeric(colMeans(data[idx, , drop = FALSE]))
    })
    # sapply: p=1 -> vector of length k; p>1 -> p x k matrix
    if (is.matrix(raw)) t(raw) else matrix(raw, ncol = 1)
  }

  centers <- compute_centers(score.scaled, tier.indices, k.current)
  sizes   <- tabulate(tier.indices, nbins = k.current)

  while (any(sizes < size.threshold) && length(sizes) > 1L) {

    sigma.pooled <- compute_pooled_cov(score.scaled, tier.indices)
    sigma.inv    <- safe_solve(sigma.pooled)

    i.min <- which.min(sizes)

    dists <- vapply(seq_len(length(sizes)), function(j) {
      if (j == i.min) return(Inf)
      mahal_dist(centers[i.min, ], centers[j, ], sigma.inv)
    }, numeric(1))

    merge.with <- which.min(dists)

    new.size   <- sizes[i.min] + sizes[merge.with]
    new.center <- (sizes[i.min]    * centers[i.min, ] +
                     sizes[merge.with] * centers[merge.with, ]) / new.size

    sizes[merge.with]    <- new.size
    centers[merge.with, ] <- new.center
    sizes   <- sizes[-i.min]
    centers <- centers[-i.min, , drop = FALSE]

    tier.indices[tier.indices == i.min] <- merge.with
    current.ids  <- sort(unique(tier.indices))
    tier.indices <- match(tier.indices, current.ids)
    k.current    <- length(sizes)
  }

  # ---- 6. Sort tiers by ascending mean score (tier 1 = lowest risk) --------
  group.assignments <- tier.indices

  tier.mean.score <- vapply(sort(unique(group.assignments)), function(ck) {
    idx <- which(group.assignments == ck)
    mean(score.mat[idx, 1])
  }, numeric(1))
  tier.order      <- order(tier.mean.score)
  tier.map        <- integer(max(group.assignments))
  tier.map[sort(unique(group.assignments))[tier.order]] <- seq_along(tier.order)
  group.assignments <- tier.map[group.assignments]

  # Recompute centers and sizes in sorted order
  k.final <- max(group.assignments)
  centers.sorted <- compute_centers(score.scaled,
                                    group.assignments,
                                    k.final)
  sizes.sorted   <- tabulate(group.assignments, nbins = k.final)

  # ---- 7. ACAT weights for nested test design --------------------------------
  # Weight w_k = cumulative mean center (unscaled) from tier k to K,
  # approximating the NCP of nested burden test k under the rare variant
  # assumption.

  # Back-transform standardized centers to original score space
  unscaled.centers <- t(t(centers.sorted) * col.sds + col.means)
  center.vals      <- rowSums(unscaled.centers)

  cum.sizes   <- rev(cumsum(rev(sizes.sorted)))
  cum.sc.sz   <- rev(cumsum(rev(center.vals * sizes.sorted)))
  cum.centers <- cum.sc.sz / cum.sizes

  min.cc <- min(cum.centers)
  ac.weights <- if (min.cc <= 0) cum.centers - min.cc + 0.1 else cum.centers

  # ---- 8. Pareto anchor identification (per-tier) --------------------------
  anchors.list <- vector("list", k.final)
  anchors.list[[1]] <- NULL   # tier 1 = default (lowest risk), no anchors
  if (k.final > 1L) {
    for (k in 2:k.final) {
      tier.scores <- score.mat[group.assignments == k, , drop = FALSE]
      if (nrow(tier.scores) > 0L) {
        anchors.list[[k]] <- find_pareto_anchors_optimized(tier.scores)
      }
    }
  }

  # ---- Return ---------------------------------------------------------------
  list(
    group.assignments = group.assignments,
    score.centers     = centers.sorted,
    ac.weights        = ac.weights,
    cluster.sizes     = sizes.sorted,
    prediction.model  = list(
      type             = "layered_pareto",
      K                = k.final,
      anchors.list     = anchors.list,
      scale.mean       = col.means,
      scale.sd         = col.sds,
      shift.vals       = shift.vals,
      transform.method = transform.method
    )
  )
}
