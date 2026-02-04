#' Run the VACANT Analysis (Returning Model Object)
#'
#' @return A list containing analysis results ($results) and the prediction model ($model).
#' @import logistf stringi mclust stats
#' @export
vacant <- function(geno,
                   score,
                   phenotype,
                   covariates = NULL,
                   test = c("uni", "multi"),
                   acat.weight = c("score", "equal"),
                   size.threshold = 10) {

  test <- match.arg(test)
  acat.weight <- match.arg(acat.weight)

  # ... (Argument Validation & Firth Step 2 logic remains the same) ...
  # [Copy Step 1 & 2 from your previous code]

  # ---- 2. Initial Firth P-value (Aggregated) ----
  gt  <- if (length(geno) == 1) {
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
                     family = "binomial",
                     control = logistf.control(maxit = 100),
                     max_retries = 10)

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
  clus <- cluster_score(score, geno, size.threshold)

  if (is.null(clus)) {
    warning("Clustering failed. Returning NULL.")
    return(NULL)
  }

  # ---- 4. Analysis ----
  stats.out <- analyze_set(pre, clus, test = test, acat.weight = acat.weight)

  # [New Return Structure]
  # We return a List containing both Statistics and the Prediction Model
  return(list(
    results = stats.out,            # The numeric vector (P-values, Betas)
    model   = clus$prediction.model,# The Pareto anchors for future use
    info    = list(                 # Metadata
      test = test,
      weight = acat.weight,
      n_vars = length(geno)
    )
  ))
}
