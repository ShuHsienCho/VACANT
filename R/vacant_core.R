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
#' @param size.threshold Integer. Minimum cluster size for K-means (default: 10).
#' @param transform.method Character. Score transformation applied before
#'   clustering: \code{"none"} (default), \code{"raw_squared"},
#'   \code{"phred_to_chisq"}, \code{"log"}, or \code{"sigmoid"}.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{\code{results}}{Named numeric vector of p-values and effect estimates.}
#'     \item{\code{model}}{Pareto staircase prediction model (for use with
#'       \code{predict_vacant_cluster()}).}
#'     \item{\code{info}}{List of metadata (test type, weights, n variants, transform).}
#'   }
#'   Returns NULL if no carriers are found or the model fails.
#'
#' @seealso \code{\link{vacant}}, \code{\link{predict_vacant_cluster}}
#' @export
vacant_core <- function(geno,
                        score,
                        phenotype,
                        covariates       = NULL,
                        test             = c("multi", "uni"),
                        acat.weight      = c("score", "equal"),
                        size.threshold   = 10,
                        transform.method = c("none", "raw_squared",
                                             "phred_to_chisq", "log", "sigmoid")) {

  test             <- match.arg(test)
  acat.weight      <- match.arg(acat.weight)
  transform.method <- match.arg(transform.method)

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
  clus <- cluster_score(score, geno,
                        size.threshold   = size.threshold,
                        transform.method = transform.method)

  if (is.null(clus)) {
    warning("Clustering failed. Returning NULL.")
    return(NULL)
  }

  # ---- 3. Statistical Analysis ----
  stats.out <- analyze_set(pre, clus,
                           test        = test,
                           acat.weight = acat.weight)

  list(
    results = stats.out,
    model   = clus$prediction.model,
    info    = list(
      test      = test,
      weight    = acat.weight,
      n_vars    = length(geno),
      transform = transform.method
    )
  )
}
