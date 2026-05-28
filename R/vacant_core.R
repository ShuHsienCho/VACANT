#' Core VACANT Analysis Engine
#'
#' Runs the VACANT statistical analysis on pre-processed R objects.
#' This is the internal engine called by \code{vacant()}. Advanced users
#' may call this directly when genotype and score data are already loaded
#' into memory as R objects.
#'
#' @param geno Character vector. One string per variant; each character
#'   in the string is one sample's allele count (0/1/2).
#' @param score Numeric matrix or vector. Annotation scores; rows = variants.
#' @param phenotype Data frame or numeric vector. Binary phenotype (0/1).
#'   If a data frame, must contain a column named \code{phenotype}.
#' @param covariates Numeric matrix or data frame. Covariates (default: NULL).
#' @param test Character. \code{"multi"} (joint multivariate Firth; default)
#'   or \code{"uni"} (sequential univariate + ACAT).
#' @param acat.weight Character. ACAT weighting: \code{"score"} (default)
#'   or \code{"equal"}.
#' @param size.threshold Integer. Minimum variant count per tier (default: 10).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{results}}{Named numeric vector of p-values and effect estimates.}
#'     \item{\code{group.assignments}}{Integer vector of tier IDs per variant.}
#'     \item{\code{cluster.sizes}}{Integer vector of tier sizes.}
#'     \item{\code{score.centers}}{Matrix of tier centers (standardized scale).}
#'     \item{\code{model}}{Pareto staircase prediction model (for use with
#'       \code{predict_vacant_cluster()}).}
#'     \item{\code{tier.summary}}{Data frame with per-tier score ranges,
#'       carrier counts, MAF, and effect estimates.}
#'     \item{\code{info}}{List of metadata (test type, weights, n variants).}
#'   }
#'   Returns NULL if no carriers are found or the model fails.
#'
#' @seealso \code{\link{vacant}}, \code{\link{predict_vacant_cluster}}
#' @export
vacant_core <- function(geno,
                        score,
                        phenotype,
                        covariates     = NULL,
                        test           = c("multi", "uni"),
                        acat.weight    = c("score", "equal"),
                        size.threshold = 10) {

  test        <- match.arg(test)
  acat.weight <- match.arg(acat.weight)

  # ---- 1. Initial Firth P-value (Aggregated Burden) ----
  gt <- if (length(geno) == 1) {
    utf8ToInt(geno) - utf8ToInt("0")
  } else {
    Reduce(`+`, lapply(geno, function(s) utf8ToInt(s) - utf8ToInt("0")))
  }
  gt.bin <- ifelse(gt > 0, 1, 0)

  if (is.null(covariates)) {
    df.firth <- cbind(phenotype = phenotype, gt.bin)
  } else {
    df.firth <- cbind(phenotype = phenotype, gt.bin, covariates)
  }

  if (sum(gt.bin) == 0) return(NULL)

  m1 <- safe_logistf(phenotype ~ ., data = df.firth,
                     family      = "binomial",
                     control     = logistf::logistf.control(maxit = 100),
                     max_retries = 10)

  if (is.null(m1)) return(NULL)

  coef1   <- m1$coefficients["gt.bin"]
  p.f     <- m1$prob["gt.bin"]
  p.firth <- ifelse(coef1 > 0, p.f / 2, 1 - p.f / 2)

  pre <- list(
    geno    = geno,
    cov.sub = covariates,
    aff.sub = phenotype,
    p.firth = p.firth
  )

  # ---- 2. Clustering ----
  clus <- cluster_score(score, size.threshold = size.threshold)

  if (is.null(clus)) {
    warning("Clustering failed. Returning NULL.")
    return(NULL)
  }

  # ---- 3. Statistical Analysis ----
  stats.out <- analyze_set(pre, clus,
                           test        = test,
                           acat.weight = acat.weight)

  # ---- 4. Tier Summary ----
  K <- max(clus$group.assignments)
  score.mat <- if (is.vector(score)) matrix(score, ncol = 1) else as.matrix(score)
  phenotype.int <- if (is.data.frame(phenotype)) as.integer(phenotype[[1]]) else as.integer(phenotype)
  n.samples <- nchar(geno[1])
  case.idx  <- which(phenotype.int == 1L)
  ctrl.idx  <- which(phenotype.int == 0L)

  tier.summary.rows <- vector("list", K)
  for (k in seq_len(K)) {
    idx <- which(clus$group.assignments == k)
    s   <- score.mat[idx, 1]

    tier.geno <- geno[idx]
    sample.any <- rep(0L, n.samples)
    for (g in tier.geno) {
      ac.vec <- utf8ToInt(g) - utf8ToInt("0")
      sample.any <- pmax(sample.any, as.integer(ac.vec > 0L))
    }

    ac.per.var <- vapply(tier.geno, function(g) sum(utf8ToInt(g) - utf8ToInt("0")), integer(1))
    maf.per.var <- ac.per.var / (2L * n.samples)

    tier.summary.rows[[k]] <- data.frame(
      tier              = k,
      n.variants        = length(idx),
      score.min         = min(s),
      score.max         = max(s),
      score.mean        = round(mean(s), 2),
      maf.mean          = signif(mean(maf.per.var), 3),
      maf.median        = signif(median(maf.per.var), 3),
      n.carriers.total  = sum(sample.any),
      n.carriers.case   = sum(sample.any[case.idx]),
      n.carriers.ctrl   = sum(sample.any[ctrl.idx]),
      carrier.freq.case = round(sum(sample.any[case.idx]) / max(length(case.idx), 1L), 4),
      carrier.freq.ctrl = round(sum(sample.any[ctrl.idx]) / max(length(ctrl.idx), 1L), 4),
      stringsAsFactors  = FALSE
    )
  }
  tier.summary.df <- do.call(rbind, tier.summary.rows)

  list(
    results           = stats.out,
    group.assignments = clus$group.assignments,
    cluster.sizes     = clus$cluster.sizes,
    score.centers     = clus$score.centers,
    model             = clus$prediction.model,
    tier.summary      = tier.summary.df,
    info              = list(
      test   = test,
      weight = acat.weight,
      n_vars = length(geno)
    )
  )
}