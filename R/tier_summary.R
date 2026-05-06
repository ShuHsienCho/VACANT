# ==============================================================================
# tier_summary.R
# Post-analysis tier characterization for VACANT results.
#
# Three functions:
#   tier_summary()    - base summary (always available, no external data needed)
#   enrich_vep()      - add VEP functional consequence breakdown (optional)
#   enrich_clinvar()  - add ClinVar classification counts (optional)
#
# Design:
#   tier_summary() produces a list with $tier (per-tier) and $variant
#   (per-variant) data frames. enrich_*() functions take this output and
#   append columns, returning the same structure. This allows chaining:
#
#     result <- vacant(...)
#     ts <- tier_summary(result, score, geno, phenotype)
#     ts <- enrich_vep(ts, vep.table)
#     ts <- enrich_clinvar(ts, clinvar.table)
# ==============================================================================


#' Summarize tier composition from a VACANT result
#'
#' Produces per-tier and per-variant summary tables from a completed VACANT
#' analysis.  The output supports downstream visualization (box plots, scatter
#' plots) and addresses reviewer concerns about post-hoc tier interpretability
#' and clinical actionability.
#'
#' @param vacant.result List. Output from \code{vacant()} or
#'   \code{vacant_core()}.  Must contain \code{group.assignments},
#'   \code{cluster.sizes}, and \code{results}.
#' @param score Numeric matrix or vector.  Raw annotation scores, same order
#'   as variants passed to \code{vacant()}.  Only variants that passed MAF
#'   filtering should be included (same length as
#'   \code{vacant.result$group.assignments}).
#' @param geno Character vector.  Genotype strings, one per variant.  Each
#'   character = one sample's allele count (0/1/2).  Same order as \code{score}.
#' @param phenotype Integer vector.  Binary phenotype (0 = control, 1 = case),
#'   one per sample.  Order must match genotype string positions.
#' @param variant.info Data frame or NULL.  Optional variant metadata
#'   (e.g., chr, pos, ref, alt).  If supplied, must have
#'   \code{nrow(variant.info) == length(geno)} and columns are prepended
#'   to the variant-level table.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{tier}}{Data frame, one row per tier.  Columns: tier,
#'       n.variants, n.carriers.total, n.carriers.case, n.carriers.control,
#'       carrier.freq.case, carrier.freq.control, score.mean, score.median,
#'       score.sd, score.q25, score.q75, score.min, score.max, and
#'       (if available from results) beta, OR, OR.lower, OR.upper, p.}
#'     \item{\code{variant}}{Data frame, one row per variant.  Columns:
#'       variant.idx, tier, score (one column per score dimension),
#'       carrier.count, carrier.count.case, carrier.count.control, plus
#'       any columns from \code{variant.info}.}
#'   }
#'
#' @export
tier_summary <- function(vacant.result,
                         score,
                         geno,
                         phenotype,
                         variant.info = NULL) {

  # ---- Input validation ------------------------------------------------------
  assignments <- vacant.result$group.assignments
  if (is.null(assignments)) {
    stop("vacant.result must contain group.assignments. ",
         "Use the updated vacant_core() that returns cluster assignments.")
  }

  if (is.vector(score) && !is.list(score)) {
    score.mat <- matrix(score, ncol = 1)
    colnames(score.mat) <- "score"
  } else {
    score.mat <- as.matrix(score)
  }
  if (is.null(colnames(score.mat))) {
    colnames(score.mat) <- paste0("score.", seq_len(ncol(score.mat)))
  }

  n.var <- length(assignments)
  if (nrow(score.mat) != n.var) {
    stop("score has ", nrow(score.mat), " rows but group.assignments has ",
         n.var, " elements.")
  }
  if (length(geno) != n.var) {
    stop("geno has ", length(geno), " elements but group.assignments has ",
         n.var, " elements.")
  }
  if (!is.null(variant.info) && nrow(variant.info) != n.var) {
    stop("variant.info has ", nrow(variant.info),
         " rows but expected ", n.var, ".")
  }

  phenotype <- as.integer(phenotype)
  n.samples <- nchar(geno[1])
  if (length(phenotype) != n.samples) {
    stop("phenotype has ", length(phenotype),
         " elements but genotype strings have ", n.samples, " characters.")
  }

  n.cases    <- sum(phenotype == 1L)
  n.controls <- sum(phenotype == 0L)
  case.idx   <- which(phenotype == 1L)
  ctrl.idx   <- which(phenotype == 0L)

  K <- max(assignments)

  # ---- Per-variant carrier counts and allele counts -------------------------
  carrier.count      <- integer(n.var)
  carrier.count.case <- integer(n.var)
  carrier.count.ctrl <- integer(n.var)
  allele.count       <- integer(n.var)

  n.chroms <- 2L * n.samples

  for (v in seq_len(n.var)) {
    ac.vec <- utf8ToInt(geno[v]) - utf8ToInt("0")
    carriers <- which(ac.vec > 0L)
    carrier.count[v]      <- length(carriers)
    carrier.count.case[v] <- sum(carriers %in% case.idx)
    carrier.count.ctrl[v] <- sum(carriers %in% ctrl.idx)
    allele.count[v]       <- sum(ac.vec)
  }

  maf <- allele.count / n.chroms

  # ---- Variant-level table --------------------------------------------------
  variant.df <- data.frame(
    variant.idx         = seq_len(n.var),
    tier                = assignments,
    score.mat,
    allele.count        = allele.count,
    maf                 = maf,
    carrier.count       = carrier.count,
    carrier.count.case  = carrier.count.case,
    carrier.count.ctrl  = carrier.count.ctrl,
    stringsAsFactors    = FALSE
  )

  if (!is.null(variant.info)) {
    variant.df <- cbind(as.data.frame(variant.info), variant.df)
  }

  # ---- Per-tier summary -----------------------------------------------------
  # Extract effect estimates from results if available
  res <- vacant.result$results

  tier.rows <- vector("list", K)
  for (k in seq_len(K)) {
    idx <- which(assignments == k)
    s   <- score.mat[idx, 1]  # primary score dimension for summary stats

    # Unique carriers per tier (sample carries >= 1 variant in this tier)
    tier.geno   <- geno[idx]
    sample.any  <- rep(0L, n.samples)
    for (g in tier.geno) {
      ac.vec <- utf8ToInt(g) - utf8ToInt("0")
      sample.any <- pmax(sample.any, as.integer(ac.vec > 0L))
    }
    n.car.total <- sum(sample.any)
    n.car.case  <- sum(sample.any[case.idx])
    n.car.ctrl  <- sum(sample.any[ctrl.idx])

    row.k <- data.frame(
      tier                 = k,
      n.variants           = length(idx),
      n.carriers.total     = n.car.total,
      n.carriers.case      = n.car.case,
      n.carriers.control   = n.car.ctrl,
      carrier.freq.case    = n.car.case / max(n.cases, 1L),
      carrier.freq.control = n.car.ctrl / max(n.controls, 1L),
      maf.mean             = mean(maf[idx]),
      maf.median           = median(maf[idx]),
      score.mean           = mean(s),
      score.median         = median(s),
      score.sd             = sd(s),
      score.q25            = quantile(s, 0.25, names = FALSE),
      score.q75            = quantile(s, 0.75, names = FALSE),
      score.min            = min(s),
      score.max            = max(s),
      stringsAsFactors     = FALSE
    )

    # Effect estimates from analyze_set (multi test)
    beta.name     <- paste0("beta.", k)
    or.name       <- paste0("OR.", k)
    or.lower.name <- paste0("OR.lower.", k)
    or.upper.name <- paste0("OR.upper.", k)
    p.name        <- paste0("p.", k)

    if (!is.null(res) && beta.name %in% names(res)) {
      row.k$beta     <- res[beta.name]
      row.k$OR       <- res[or.name]
      row.k$OR.lower <- res[or.lower.name]
      row.k$OR.upper <- res[or.upper.name]
      row.k$p        <- res[p.name]
    }

    tier.rows[[k]] <- row.k
  }

  tier.df <- do.call(rbind, tier.rows)
  rownames(tier.df) <- NULL

  structure(
    list(tier = tier.df, variant = variant.df),
    class = "vacant_tier_summary"
  )
}


#' Print method for VACANT tier summary
#' @param x A \code{vacant_tier_summary} object.
#' @param ... Ignored.
#' @export
print.vacant_tier_summary <- function(x, ...) {
  cat("VACANT Tier Summary\n")
  cat(sprintf("  %d tiers, %d variants\n\n", nrow(x$tier), nrow(x$variant)))
  print(x$tier, row.names = FALSE)
  invisible(x)
}


# ==============================================================================
# Enrichment: VEP functional consequence
# ==============================================================================

#' Add VEP functional consequence to a tier summary
#'
#' Joins a VEP annotation table to the variant-level output of
#' \code{tier_summary()}, then computes per-tier consequence breakdown.
#'
#' @param ts A \code{vacant_tier_summary} object (output of
#'   \code{tier_summary()}).
#' @param vep.table Data frame.  Must contain a column named
#'   \code{consequence} (character; e.g., "missense_variant",
#'   "stop_gained") and columns matching \code{join.cols} for variant
#'   identification.
#' @param join.cols Character vector.  Column names present in both
#'   \code{ts$variant} (via \code{variant.info}) and \code{vep.table} to
#'   use as join keys.  Default: \code{c("chr", "pos", "ref", "alt")}.
#'
#' @return The input \code{ts} with:
#'   \itemize{
#'     \item \code{ts$variant$consequence} added.
#'     \item \code{ts$tier} gains one column per unique consequence
#'       (count), plus \code{n.consequence.annotated} (variants matched).
#'   }
#'
#' @importFrom dplyr left_join
#' @export
enrich_vep <- function(ts,
                       vep.table,
                       join.cols = c("chr", "pos", "ref", "alt")) {

  if (!inherits(ts, "vacant_tier_summary")) {
    stop("ts must be output of tier_summary().")
  }

  missing.var <- setdiff(join.cols, colnames(ts$variant))
  if (length(missing.var) > 0L) {
    stop("ts$variant is missing join columns: ",
         paste(missing.var, collapse = ", "),
         ". Supply variant.info with these columns to tier_summary().")
  }
  missing.vep <- setdiff(join.cols, colnames(vep.table))
  if (length(missing.vep) > 0L) {
    stop("vep.table is missing join columns: ",
         paste(missing.vep, collapse = ", "))
  }
  if (!"consequence" %in% colnames(vep.table)) {
    stop("vep.table must contain a 'consequence' column.")
  }

  # Join
  vep.sub <- vep.table[, c(join.cols, "consequence"), drop = FALSE]
  ts$variant <- suppressMessages(
    dplyr::left_join(ts$variant, vep.sub, by = join.cols)
  )

  # Per-tier consequence breakdown
  K <- max(ts$variant$tier)
  csq.all <- sort(unique(na.omit(ts$variant$consequence)))

  for (k in seq_len(K)) {
    tier.csq <- ts$variant$consequence[ts$variant$tier == k]
    ts$tier$n.consequence.annotated[k] <- sum(!is.na(tier.csq))
    for (csq in csq.all) {
      col.name <- paste0("n.", gsub("[^a-zA-Z0-9]", ".", csq))
      ts$tier[[col.name]][k] <- sum(tier.csq == csq, na.rm = TRUE)
    }
  }

  ts
}


# ==============================================================================
# Enrichment: ClinVar classification
# ==============================================================================

#' Add ClinVar clinical significance to a tier summary
#'
#' Joins a ClinVar annotation table to the variant-level output of
#' \code{tier_summary()}, then computes per-tier classification counts.
#'
#' @param ts A \code{vacant_tier_summary} object.
#' @param clinvar.table Data frame.  Must contain a column named
#'   \code{clinvar.class} (character; e.g., "Pathogenic",
#'   "Likely_pathogenic", "VUS", "Likely_benign", "Benign") and
#'   columns matching \code{join.cols}.
#' @param join.cols Character vector.  Join keys.  Default:
#'   \code{c("chr", "pos", "ref", "alt")}.
#'
#' @return The input \code{ts} with:
#'   \itemize{
#'     \item \code{ts$variant$clinvar.class} added.
#'     \item \code{ts$tier} gains columns: n.pathogenic, n.likely.path,
#'       n.vus, n.likely.benign, n.benign, n.clinvar.annotated.
#'   }
#'
#' @importFrom dplyr left_join
#' @export
enrich_clinvar <- function(ts,
                           clinvar.table,
                           join.cols = c("chr", "pos", "ref", "alt")) {

  if (!inherits(ts, "vacant_tier_summary")) {
    stop("ts must be output of tier_summary().")
  }

  missing.var <- setdiff(join.cols, colnames(ts$variant))
  if (length(missing.var) > 0L) {
    stop("ts$variant is missing join columns: ",
         paste(missing.var, collapse = ", "),
         ". Supply variant.info with these columns to tier_summary().")
  }
  missing.cv <- setdiff(join.cols, colnames(clinvar.table))
  if (length(missing.cv) > 0L) {
    stop("clinvar.table is missing join columns: ",
         paste(missing.cv, collapse = ", "))
  }
  if (!"clinvar.class" %in% colnames(clinvar.table)) {
    stop("clinvar.table must contain a 'clinvar.class' column.")
  }

  # Standardize labels
  label.map <- c(
    "Pathogenic"        = "P",
    "Likely_pathogenic" = "LP",
    "Likely pathogenic" = "LP",
    "Uncertain_significance" = "VUS",
    "Uncertain significance" = "VUS",
    "VUS"               = "VUS",
    "Likely_benign"     = "LB",
    "Likely benign"     = "LB",
    "Benign"            = "B"
  )

  cv.sub <- clinvar.table[, c(join.cols, "clinvar.class"), drop = FALSE]
  cv.sub$clinvar.std <- ifelse(
    cv.sub$clinvar.class %in% names(label.map),
    label.map[cv.sub$clinvar.class],
    "Other"
  )

  ts$variant <- suppressMessages(
    dplyr::left_join(ts$variant,
                     cv.sub[, c(join.cols, "clinvar.class", "clinvar.std")],
                     by = join.cols)
  )

  # Per-tier ClinVar counts
  K <- max(ts$variant$tier)
  std.levels <- c("P", "LP", "VUS", "LB", "B", "Other")
  std.col.names <- c("n.pathogenic", "n.likely.path", "n.vus",
                     "n.likely.benign", "n.benign", "n.other.clinvar")

  for (k in seq_len(K)) {
    tier.cv <- ts$variant$clinvar.std[ts$variant$tier == k]
    ts$tier$n.clinvar.annotated[k] <- sum(!is.na(tier.cv))
    for (i in seq_along(std.levels)) {
      ts$tier[[std.col.names[i]]][k] <- sum(tier.cv == std.levels[i],
                                            na.rm = TRUE)
    }
  }

  ts
}
