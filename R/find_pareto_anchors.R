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