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
#' @import logistf
#' @export
analyze_set <- function(pre, clus, 
                        test = c("uni", "multi"), 
                        acat.weight = c("score", "equal")) {
  
  test <- match.arg(test)
  acat.weight <- match.arg(acat.weight)
  
  # ---- 1. Setup ----
  k <- length(clus$cluster.sizes)
  aff = pre$aff.sub
  covariates = pre$cov.sub
  
  # Calculate weighted average score (using L2 magnitude) for summary
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
    p.score.uni   <- numeric(k)
    
    for (i in seq_len(k)) {
      keep   <- clus.order[ i:k ]
      sel    <- clus$group.assignments %in% keep
      ac.sel <- pre$geno.maf[sel]
      
      gt.i   <- if (length(ac.sel) == 1) {
        utf8ToInt(ac.sel) - utf8ToInt("0")
      } else {
        Reduce(`+`, lapply(ac.sel,function(s) utf8ToInt(s) - utf8ToInt("0")))
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
        
        # [Output] Beta and P for Univariate
        # We can store them if needed, but standard VACANT focused on ACAT.
        # Let's just keep p.firth logic for ACAT.
        coef1.i   <- m1.i$coefficients["gt.bin.i"]
        p0        <- m1.i$prob["gt.bin.i"]
        p.firth.i <- ifelse(coef1.i > 0, p0/2, 1 - p0/2)
      } else {
        p.firth.i <- .Machine$double.neg.eps
      }
      
      p.score.uni[i] <- max(p.firth.i, .Machine$double.neg.eps)
    }
    
    weights <- if (acat.weight == "score") clus$ac.weights else rep(1/k, k)
    acat.uni <- acat_p(acat_t(p.score.uni, weights), weights)
    
    results.out <- c(
      results.out,
      setNames(p.score.uni, paste0("p.score.uni", seq_len(k))),
      acat.uni = acat.uni
    )
  }
  
  # ---- 3. Multivariate Version (Modified) ----
  if (test == "multi") {
    geno.mat <- do.call(cbind, lapply(seq_len(k), function(i) {
      sel <- clus$group.assignments == i
      ac.sel <- pre$geno.maf[sel]
      geno   <- if (length(ac.sel) == 1) {
        utf8ToInt(ac.sel) - utf8ToInt("0")
      } else {
        Reduce(`+`, lapply(ac.sel,function(s) utf8ToInt(s) - utf8ToInt("0")))
      } 
      ifelse(geno > 0, 1, 0)
    }))
    
    if (is.null(covariates)) {
      df2 <- cbind(phenotype = aff, geno.mat)
    } else {
      df2 <- cbind(phenotype = aff, geno.mat, covariates)
    }
    
    # Name genotype columns gt1, gt2, ..., gtk
    gt_cols <- paste0("gt", seq(k))
    colnames(df2)[2:(k+1)] = gt_cols
    
    # Null model
    m0 <- safe_logistf(phenotype ~ ., data = df2[, -c(2:(k+1)), drop=FALSE], 
                       family = "binomial", control = logistf.control(maxit=100))
    
    # Full multivariate model
    m1m <- safe_logistf(phenotype ~ ., data = df2,
                        family = "binomial",
                        control = logistf.control(maxit = 100),
                        max_retries = 10)
    
    # [New] Extract Per-Cluster Statistics
    cluster_stats <- numeric()
    if (!is.null(m1m)) {
      for (i in 1:k) {
        term <- paste0("gt", i)
        if (term %in% names(m1m$coefficients)) {
          beta <- m1m$coefficients[term]
          p_val_raw <- m1m$prob[term]
          # One-sided P-value (assuming risk increasing)
          p_one_sided <- ifelse(beta > 0, p_val_raw/2, 1 - p_val_raw/2)
          if(p_one_sided == 0) p_one_sided <- .Machine$double.neg.eps
          
          cluster_stats <- c(cluster_stats, 
                             setNames(beta, paste0("beta.", i)), 
                             setNames(p_one_sided, paste0("p.", i)))
        } else {
          # Should not happen unless variable dropped due to constancy
          cluster_stats <- c(cluster_stats, 
                             setNames(0, paste0("beta.", i)), 
                             setNames(1, paste0("p.", i)))
        }
      }
      
      # Prepare ACAT inputs
      p.multi <- m1m$prob[grepl("gt", names(m1m$prob))]
      coef.multi <- m1m$coefficients[grepl("gt", names(m1m$coefficients))]
      p.multi <- ifelse(coef.multi > 0, p.multi/2, 1 - p.multi/2)
      p.multi[p.multi == 0] <- .Machine$double.neg.eps
      
      w.equal <- rep(1/k, k) 
      acat.multi <- acat_p(acat_t(p.multi, w.equal), w.equal)
      p.lrt.multi <- anova(m0, m1m, method = "PLR")$pval
      
    } else {
      # Fallback if model failed completely
      acat.multi <- 1
      p.lrt.multi <- 1
      cluster_stats <- rep(NA, k*2) # Fill placeholders
    }

    results.out <- c(
      results.out,
      cluster_stats, # [Added] Individual betas and p-values
      acat.multi = acat.multi,
      p.lrt.multi = p.lrt.multi
    )
  }
  
  return(results.out)
}