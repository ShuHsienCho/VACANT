#' Calculate ACAT test statistic
#'
#' Given a vector of p-values and corresponding weights, computes the aggregated
#' Cauchy test statistic.
#'
#' @param p Numeric vector of p-values (each in (0,1)).
#' @param w Numeric vector of positive weights of the same length as \code{p}.
#' @return Numeric scalar: the ACAT test statistic.
#' @export
acat_t <- function(p, w) {
  sum(w * tan((0.5 - p) * pi))
}

#' Convert ACAT statistic back to a p-value
#'
#' Given an ACAT test statistic and its weights, returns the combined p-value.
#'
#' @param t Numeric scalar: ACAT test statistic.
#' @param w Numeric vector of weights used to compute \code{t}.
#' @return Numeric scalar: the combined p-value.
#' @export
acat_p <- function(t, w) {
  0.5 - atan(t / sum(w)) / pi
}

#' Perform k-means clustering with retry on error
#'
#' Attempts k-means clustering on \code{data}. If an error occurs,
#' prints the error message, waits one second, and retries recursively.
#'
#' @param data Numeric matrix or data frame.
#' @param mod A model-like object containing a \code{classification} vector.
#' @return An object of class \code{\link[stats]{kmeans}}.
#' @importFrom stats kmeans
#' @export
perform_kmeans <- function(data, mod) {
  tryCatch({
    k <- length(unique(mod$classification))
    initial <- stats::kmeans(data, centers = k)
    centers_sorted <- sort(initial$centers)
    stats::kmeans(data, centers = centers_sorted)
  }, error = function(e) {
    message("Error in perform_kmeans: ", e$message, "; retrying in 1 second...")
    Sys.sleep(1)
    perform_kmeans(data, mod)
  })
}

#' Extract a single character from a genotype string
#'
#' Given a genotype string vector and an index, returns the allele at that position.
#'
#' @param x Character vector of genotype strings.
#' @param idx Integer. Position of allele to extract.
#' @return Character scalar of the extracted allele.
#' @import stringi
#' @export
#' 
extract_sub <- function(x, idx) {
  stringi::stri_c(stringi::stri_sub(x, from = idx, to = idx), collapse = "")
}


#' Safe Firth Logistic Regression with Retry on Singular Fisher Information
#'
#' Wraps \code{logistf::logistf} to perform Firth-corrected logistic regression,
#' automatically retrying if the Fisher information matrix becomes singular.
#'
#' @param formula A \code{formula} specifying the logistic regression model.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param family Character or function. The model family (default "binomial").
#' @param control A control object created by \code{logistf.control()}.
#' @param max_retries Integer. Maximum number of retry attempts (default 5).
#' @return An object of class \code{logistf}.
#' @importFrom logistf logistf
#' @export
safe_logistf <- function(formula, data, family="binomial", control, max_retries=5) {
  attempt <- 1
  repeat {
    warn_msg <- NULL
    fit <- withCallingHandlers(
      try(logistf::logistf(formula, data=data, family=family, control=control), silent=TRUE),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("Determinant of Fisher information matrix was numerically 0", msg)) {
          warn_msg <<- msg
          invokeRestart("muffleWarning")
        }
      }
    )
    if (is.null(warn_msg) && !inherits(fit, "try-error")) {
      return(fit)
    }
    message(sprintf("Warning in logistf (attempt %d/%d): %s. Retrying...",
                    attempt, max_retries, warn_msg))
    attempt <- attempt + 1
    if (attempt > max_retries) {
      stop("logistf failed after ", max_retries, " attempts due to singular Fisher information.")
    }
  }
}

#' Efficiently find Pareto Anchors (Lower-Left Frontier)
#'
#' This helper function implements a "Sort-and-Scan" algorithm to identify
#' non-dominated points (Pareto frontier) in O(N log N) or O(N * k) time.
#'
#' @param points Numeric matrix of points (e.g., scores of risk cluster variants).
#' @return A matrix of anchor points defining the lower-left boundary, 
#'   subset from the input \code{points}.
#' @export
find_pareto_anchors_optimized <- function(points) {
  if (is.null(points) || nrow(points) == 0) return(NULL)
  if (nrow(points) == 1) return(points)
  
  # Deduplicate to reduce computation
  pts <- unique(points)
  
  # 1. Sort: Sort by the first dimension (X), then subsequent dimensions
  ord <- do.call(order, as.data.frame(pts))
  pts.sorted <- pts[ord, , drop = FALSE]
  
  # 2. Fast Scan
  # Pre-allocate memory (worst case: all points are anchors)
  anchors <- matrix(NA_real_, nrow = nrow(pts.sorted), ncol = ncol(pts.sorted))
  n.anchors <- 0
  
  # Algorithm:
  # Since points are sorted by X, we only need to check if the current point
  # is dominated by any *previously found* anchors.
  
  # Optimization for 2D case (O(N))
  if (ncol(pts) == 2) {
    min.y.so.far <- Inf
    for (i in seq_len(nrow(pts.sorted))) {
      curr.y <- pts.sorted[i, 2]
      # Logic: If Y is smaller than any Y seen so far, it is a new anchor
      # because its X is guaranteed to be larger than previous anchors.
      if (curr.y < min.y.so.far) {
        n.anchors <- n.anchors + 1
        anchors[n.anchors, ] <- pts.sorted[i, ]
        min.y.so.far <- curr.y
      }
    }
  } else {
    # General N-D case (O(N * k))
    for (i in seq_len(nrow(pts.sorted))) {
      pt <- pts.sorted[i, ]
      is.dominated <- FALSE
      
      if (n.anchors > 0) {
        # Check against existing anchors
        # Dominated if: Anchor <= pt (all dimensions)
        curr.anchors <- anchors[1:n.anchors, , drop = FALSE]
        diffs <- sweep(curr.anchors, 2, pt, "-")
        
        # If any anchor is <= current point in all dimensions, current point is dominated
        if (any(rowSums(diffs <= 1e-9) == ncol(pts))) {
          is.dominated <- TRUE
        }
      }
      
      if (!is.dominated) {
        n.anchors <- n.anchors + 1
        anchors[n.anchors, ] <- pt
      }
    }
  }
  
  # Return valid rows
  return(anchors[1:n.anchors, , drop = FALSE])
}