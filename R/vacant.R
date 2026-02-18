#' Run the VACANT Analysis (Returning Model Object)
#'
#' @param geno Character vector. Genotype strings (one string per variant, containing genotypes for all samples).
#' @param score Numeric matrix or vector. Annotation scores.
#' @param phenotype Numeric vector. Phenotypes.
#' @param covariates Numeric matrix or data.frame. Covariates (default NULL).
#' @param test Character. "uni" or "multi" (default).
#' @param acat.weight Character. "score" or "equal" (default).
#' @param size.threshold Integer. Minimum cluster size.
#' @param transform.method Character. "none", "raw_squared", "phred_to_chisq", or "log".
#' @return A list containing analysis results ($results) and the prediction model ($model).
#' @import logistf stringi mclust stats
#' @export
vacant <- function(geno,
                   score,
                   phenotype,
                   covariates = NULL,
                   test = c("uni", "multi"),
                   acat.weight = c("score", "equal"),
                   size.threshold = 10,
                   transform.method = c("none", "raw_squared", "phred_to_chisq", "log", "sigmoid")) {

  test <- match.arg(test)
  acat.weight <- match.arg(acat.weight)
  transform.method <- match.arg(transform.method)

  # ---- 2. Initial Firth P-value (Aggregated) ----
  # Calculate Burden (GT sum) for initial association check
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

  # Safety check: ensure safe_logistf is available (helper function)
  if (exists("safe_logistf")) {
    m1 <- safe_logistf(phenotype ~ ., data = df.firth,
                       family = "binomial",
                       control = logistf.control(maxit = 100),
                       max_retries = 10)
  } else {
    # Fallback to standard logistf if helper missing
    m1 <- tryCatch({
      logistf::logistf(phenotype ~ ., data = df.firth,
                       control = logistf.control(maxit = 100))
    }, error = function(e) NULL)
  }

  if (is.null(m1)) return(NULL) # Safety check

  coef1  <- m1$coefficients["gt.bin"]
  p.f    <- m1$prob["gt.bin"]
  p.firth <- ifelse(coef1 > 0, p.f/2, 1 - p.f/2)

  pre <- list(
    geno         = geno,
    cov.sub      = covariates,
    aff.sub      = phenotype,
    p.firth      = p.firth
  )

  # ---- 3. Clustering ----
  # Pass transform.method to cluster_score
  # Ensure cluster_score uses the updated signature
  clus <- cluster_score(score, geno, size.threshold = size.threshold, transform.method = transform.method)

  if (is.null(clus)) {
    warning("Clustering failed. Returning NULL.")
    return(NULL)
  }

  # ---- 4. Analysis ----
  # Ensure analyze_set function exists in the package/environment
  if (!exists("analyze_set")) {
    stop("Function 'analyze_set' not found. Please ensure package is loaded correctly.")
  }
  stats.out <- analyze_set(pre, clus, test = test, acat.weight = acat.weight)

  # [New Return Structure]
  # We return a List containing both Statistics and the Prediction Model
  return(list(
    results = stats.out,            # The numeric vector (P-values, Betas)
    model   = clus$prediction.model,# The Pareto anchors for future use
    info    = list(                 # Metadata
      test = test,
      weight = acat.weight,
      n_vars = length(geno),
      transform = transform.method
    )
  ))
}
