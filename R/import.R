# ==============================================================================
# imports.R
# Centralized import declarations to avoid namespace conflicts.
# All @importFrom statements live here; individual functions do NOT
# use broad @import to prevent masking between data.table / dplyr / mclust.
# ==============================================================================

#' @importFrom data.table fread fwrite rbindlist
#' @importFrom dplyr left_join
#' @importFrom parallel mclapply
#' @importFrom stats kmeans na.omit
#' @importFrom stringi stri_c stri_sub stri_count_fixed stri_extract_all_regex
#' @importFrom logistf logistf logistf.control
#' @importFrom mclust Mclust mclustBIC mclustBICupdate hcRandomPairs
NULL
