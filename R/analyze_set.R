#' Analyze score-based clusters with Firth regression
#'
#' Performs either sequential univariate tests or a joint multivariate test.
#'
#' @param pre List. Pre-processing outputs (phenotype, genotypes, null model info).
#' @param clus List. Output from \code{cluster_score()}, containing assignments
#'   and calculated scalar weights (\code{ac.weights}).
#' @param test Character. "uni" or "multi".
#' @param acat.weight Character. "score" or "equal".
#'
#' @return Named numeric vector with results.
#' @importFrom logistf logistf logistf.control
#' @export
analyze_set <- function(pre, clus,
                        test = c("uni", "multi"),
                        acat.weight = c("score", "equal")) {

  test <- match.arg(test)
  acat.weight <- match.arg(acat.weight)

  # ---- 1. Setup ----
  k <- length(clus$cluster.sizes)
  aff <- pre$aff.sub
  covariates <- pre$cov.sub

  # Weighted average score (L2 magnitude of centers) for summary
  center.mags <- sqrt(rowSums(clus$score.centers^2))
  tx.score <- sum(center.mags * clus$cluster.sizes) / sum(clus$cluster.sizes)

  results.out <- c(
    p.firth    = pre$p.firth,
    clusters   = k,
    tx.score   = tx.score
  )

  # ---- 2. Univariate Version ----
  if (test == "uni") {
    clus.order <- seq_len(k)
    p.score.uni <- numeric(k)

    for (i in seq_len(k)) {
      keep   <- clus.order[i:k]
      sel    <- clus$group.assignments %in% keep
      ac.sel <- pre$geno[sel]

      gt.i <- if (length(ac.sel) == 1) {
        utf8ToInt(ac.sel) - utf8ToInt("0")
      } else {
        Reduce(`+`, lapply(ac.sel, function(s) utf8ToInt(s) - utf8ToInt("0")))
      }
      gt.bin.i <- ifelse(gt.i > 0, 1, 0)

      if (is.null(covariates)) {
        df.i <- cbind(phenotype = aff, gt.bin.i)
      } else {
        df.i <- cbind(phenotype = aff, gt.bin.i, covariates)
      }

      if (sum(gt.bin.i) > 0) {
        m1.i <- safe_logistf(phenotype ~ ., data = df.i,
                             family = "binomial",
                             control = logistf.control(maxit = 100),
                             max_retries = 10)

        coef1.i   <- m1.i$coefficients["gt.bin.i"]
        p0        <- m1.i$prob["gt.bin.i"]
        p.firth.i <- ifelse(coef1.i > 0, p0 / 2, 1 - p0 / 2)
      } else {
        p.firth.i <- .Machine$double.neg.eps
      }

      p.score.uni[i] <- max(p.firth.i, .Machine$double.neg.eps)
    }

    weights <- if (acat.weight == "score") clus$ac.weights else rep(1 / k, k)
    acat.uni <- acat_p(acat_t(p.score.uni, weights), weights)

    results.out <- c(
      results.out,
      setNames(p.score.uni, paste0("p.score.uni", seq_len(k))),
      acat.uni = acat.uni
    )
  }

  # ---- 3. Multivariate Version ----
  if (test == "multi") {
    geno.mat <- do.call(cbind, lapply(seq_len(k), function(i) {
      sel    <- clus$group.assignments == i
      ac.sel <- pre$geno[sel]
      geno   <- if (length(ac.sel) == 1) {
        utf8ToInt(ac.sel) - utf8ToInt("0")
      } else {
        Reduce(`+`, lapply(ac.sel, function(s) utf8ToInt(s) - utf8ToInt("0")))
      }
      ifelse(geno > 0, 1, 0)
    }))

    if (is.null(covariates)) {
      df2 <- cbind(phenotype = aff, geno.mat)
    } else {
      df2 <- cbind(phenotype = aff, geno.mat, covariates)
    }

    gt_cols <- paste0("gt", seq_len(k))
    colnames(df2)[2:(k + 1)] <- gt_cols

    # Null model (covariates only)
    m0 <- safe_logistf(phenotype ~ ., data = df2[, -c(2:(k + 1)), drop = FALSE],
                       family = "binomial", control = logistf.control(maxit = 100))

    # Full multivariate model
    m1m <- safe_logistf(phenotype ~ ., data = df2,
                        family = "binomial",
                        control = logistf.control(maxit = 100),
                        max_retries = 10)

    # Per-tier statistics: beta, p-value, OR, 95% CI
    cluster_stats <- numeric()
    if (!is.null(m1m)) {
      for (i in seq_len(k)) {
        term <- paste0("gt", i)
        if (term %in% names(m1m$coefficients)) {
          beta      <- m1m$coefficients[term]
          p_val_raw <- m1m$prob[term]
          p_one_sided <- ifelse(beta > 0, p_val_raw / 2, 1 - p_val_raw / 2)
          if (p_one_sided == 0) p_one_sided <- .Machine$double.neg.eps

          # OR and 95% CI from logistf profile likelihood CI
          idx.term  <- which(names(m1m$coefficients) == term)
          ci.lower  <- m1m$ci.lower[idx.term]
          ci.upper  <- m1m$ci.upper[idx.term]
          or.est    <- exp(beta)
          or.lower  <- exp(ci.lower)
          or.upper  <- exp(ci.upper)

          cluster_stats <- c(cluster_stats,
                             setNames(beta,        paste0("beta.", i)),
                             setNames(p_one_sided, paste0("p.", i)),
                             setNames(or.est,      paste0("OR.", i)),
                             setNames(or.lower,    paste0("OR.lower.", i)),
                             setNames(or.upper,    paste0("OR.upper.", i)))
        } else {
          cluster_stats <- c(cluster_stats,
                             setNames(0,  paste0("beta.", i)),
                             setNames(1,  paste0("p.", i)),
                             setNames(1,  paste0("OR.", i)),
                             setNames(NA, paste0("OR.lower.", i)),
                             setNames(NA, paste0("OR.upper.", i)))
        }
      }

      # Gene-level p-value: penalized LRT only
      p.lrt.multi <- anova(m0, m1m, method = "PLR")$pval

    } else {
      p.lrt.multi <- 1
      cluster_stats <- rep(NA, k * 5)
    }

    results.out <- c(
      results.out,
      cluster_stats,
      p.lrt.multi = p.lrt.multi
    )
  }

  return(results.out)
}
