# =============================================================================
# COMPACT — Comparative Pattern Counts
# compact_run.R
#
# Core analysis functions for identifying and comparing gene expression
# patterns across biological conditions.
#
# Author: Lakshmi Kuttippurathu, PhD
# License: MIT
# Repository: https://github.com/lakshmikc/COMPACT
# =============================================================================


# ---- Dependencies -----------------------------------------------------------
# install.packages(c("gtools", "ggplot2", "reshape2", "pheatmap", "RColorBrewer"))


#' Generate Synthetic Example Data for COMPACT
#'
#' Creates a minimal synthetic gene expression dataset for demonstrating
#' COMPACT analysis. Simulates a time-series experiment across multiple
#' dietary conditions (chow, high-fat, ethanol) with 3 time points.
#'
#' @param n_genes   Integer. Number of genes to simulate (default: 200)
#' @param n_samples Integer. Number of replicates per group (default: 3)
#' @param seed      Integer. Random seed for reproducibility (default: 42)
#'
#' @return A list with:
#'   \item{expr}{Numeric matrix of log2 fold-change values (genes x samples)}
#'   \item{sample_info}{Data frame describing each sample's group membership}
#'
#' @examples
#' example_data <- compact_simulate()
#' str(example_data)
#'
#' @export
compact_simulate <- function(n_genes = 200, n_samples = 3, seed = 42) {
  set.seed(seed)

  diets      <- c("CHO", "HF", "EtOH")
  timepoints <- c("1h", "6h", "24h")

  groups    <- expand.grid(Diet = diets, Time = timepoints, stringsAsFactors = FALSE)
  group_ids <- paste(groups$Diet, groups$Time, sep = "_")

  total_samples <- nrow(groups) * n_samples
  sample_names  <- paste0(rep(group_ids, each = n_samples), "_rep", seq_len(n_samples))

  # Background noise
  expr <- matrix(
    rnorm(n_genes * total_samples, mean = 0, sd = 0.5),
    nrow     = n_genes,
    ncol     = total_samples,
    dimnames = list(paste0("Gene_", seq_len(n_genes)), sample_names)
  )

  # Add structured signal to a subset of genes
  cho_24h_cols <- grep("CHO_24h", colnames(expr))
  hf_late_cols <- grep("HF_6h|HF_24h", colnames(expr))
  etoh_cols    <- grep("EtOH", colnames(expr))

  expr[1:20,  cho_24h_cols] <- rnorm(20 * length(cho_24h_cols), mean =  1.2, sd = 0.3)
  expr[21:40, hf_late_cols] <- rnorm(20 * length(hf_late_cols), mean =  0.9, sd = 0.3)
  expr[41:60, etoh_cols]    <- rnorm(20 * length(etoh_cols),    mean = -0.9, sd = 0.3)

  sample_info <- data.frame(
    sample_id = sample_names,
    diet      = rep(rep(diets, each = length(timepoints)), each = n_samples),
    timepoint = rep(rep(timepoints, times = length(diets)), each = n_samples),
    stringsAsFactors = FALSE
  )

  message("Simulated: ", n_genes, " genes x ", total_samples, " samples")
  return(list(expr = expr, sample_info = sample_info))
}


# --- Internal helpers ---------------------------------------------------------

#' Calculate per-group mean expression
#' @keywords internal
.calc_group_averages <- function(expr, sample_info) {
  groups    <- unique(sample_info[, c("diet", "timepoint")])
  group_ids <- paste(groups$diet, groups$timepoint, sep = "_")

  avg_matrix <- matrix(
    NA, nrow = nrow(expr), ncol = nrow(groups),
    dimnames = list(rownames(expr), group_ids)
  )

  for (i in seq_len(nrow(groups))) {
    idx      <- which(sample_info$diet      == groups$diet[i] &
                      sample_info$timepoint == groups$timepoint[i])
    col_idx  <- which(colnames(expr) %in% sample_info$sample_id[idx])
    avg_matrix[, i] <- rowMeans(expr[, col_idx, drop = FALSE])
  }

  return(list(avg_expr = avg_matrix, group_labels = groups))
}


#' Discretize expression into -1 / 0 / +1 ternary patterns
#' @keywords internal
.discretize_patterns <- function(avg_expr, group_labels, fc_threshold) {
  discrete <- apply(avg_expr, 2, function(col) {
    ifelse(col >  fc_threshold,  1L,
    ifelse(col < -fc_threshold, -1L, 0L))
  })
  rownames(discrete) <- rownames(avg_expr)
  return(list(discrete = discrete, group_labels = group_labels))
}


#' Build per-diet pattern strings (e.g. "1_0_-1") for each gene
#' @keywords internal
.build_pattern_strings <- function(pattern_data, diets) {
  pattern_strings <- list()
  discrete        <- pattern_data$discrete
  group_labels    <- pattern_data$group_labels

  for (d in diets) {
    d_cols               <- which(group_labels$diet == d)
    d_matrix             <- discrete[, d_cols, drop = FALSE]
    pattern_strings[[d]] <- apply(d_matrix, 1, paste, collapse = "_")
  }
  return(pattern_strings)
}


#' Build pairwise COMPACT count matrix between two diet conditions
#' @keywords internal
.build_compact_matrix <- function(pattern_strings, diet1, diet2, symmetric) {
  p1 <- pattern_strings[[diet1]]
  p2 <- pattern_strings[[diet2]]

  if (symmetric) {
    all_p <- sort(unique(c(p1, p2)))
    row_p <- col_p <- all_p
  } else {
    row_p <- sort(unique(p1))
    col_p <- sort(unique(p2))
  }

  mat <- matrix(0L,
    nrow = length(row_p), ncol = length(col_p),
    dimnames = list(paste0(diet1, ":", row_p), paste0(diet2, ":", col_p))
  )

  for (i in seq_along(row_p)) {
    for (j in seq_along(col_p)) {
      mat[i, j] <- sum(p1 == row_p[i] & p2 == col_p[j])
    }
  }
  return(mat)
}


# --- Main function ------------------------------------------------------------

#' Run COMPACT Analysis
#'
#' Main entry point for COMPACT (Comparative Pattern Counts) analysis.
#' Given a log2 fold-change expression matrix and sample metadata, this function:
#'
#' \enumerate{
#'   \item Averages expression values across replicates within each group
#'   \item Discretizes values into ternary expression patterns (-1, 0, +1)
#'   \item Builds a compact pattern string for each gene in each condition
#'   \item Computes pairwise COMPACT matrices comparing all condition pairs
#' }
#'
#' @param expr         Numeric matrix of log2 fold-change values (genes x samples).
#'                     Rows are genes, columns are samples.
#' @param sample_info  Data frame with at minimum these columns:
#'                     \code{sample_id}, \code{diet}, \code{timepoint}.
#'                     Each row describes one sample in \code{expr}.
#' @param fc_threshold Log2 fold-change cutoff for calling a gene up or down.
#'                     Default: \code{log2(1.5)} ≈ 0.585.
#' @param symmetric    Logical. If \code{TRUE}, COMPACT matrices use the union
#'                     of patterns from both conditions as row and column labels,
#'                     yielding a square matrix. Default: \code{FALSE}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{compact_matrices}}{Named list of pairwise count matrices, one
#'     per diet pair (e.g., \code{"CHO:HF"}, \code{"CHO:EtOH"}, \code{"HF:EtOH"})}
#'   \item{\code{pattern_strings}}{Named list of per-gene pattern strings per diet}
#'   \item{\code{pattern_counts}}{Named list of pattern frequency tables per diet}
#'   \item{\code{avg_expr}}{Group-averaged expression matrix}
#'   \item{\code{discrete}}{Discretized (-1/0/+1) expression matrix}
#' }
#'
#' @examples
#' # 1. Simulate example data
#' example_data <- compact_simulate()
#'
#' # 2. Run analysis
#' results <- compact_run(
#'   expr        = example_data$expr,
#'   sample_info = example_data$sample_info
#' )
#'
#' # 3. View pattern counts for CHO
#' results$pattern_counts[["CHO"]]
#'
#' # 4. View CHO vs HF COMPACT matrix
#' results$compact_matrices[["CHO:HF"]]
#'
#' # 5. Visualize results
#' compact_plot(results)
#'
#' @seealso \code{\link{compact_simulate}}, \code{\link{compact_plot}}
#' @export
compact_run <- function(expr,
                        sample_info,
                        fc_threshold = log2(1.5),
                        symmetric    = FALSE) {

  # Input validation
  req_cols <- c("sample_id", "diet", "timepoint")
  missing  <- setdiff(req_cols, colnames(sample_info))
  if (length(missing) > 0) {
    stop("sample_info is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("'expr' must be a numeric matrix or data frame.")
  }
  missing_samples <- setdiff(sample_info$sample_id, colnames(expr))
  if (length(missing_samples) > 0) {
    stop("These sample_ids in sample_info were not found in expr: ",
         paste(head(missing_samples, 5), collapse = ", "))
  }

  diets <- unique(sample_info$diet)

  message("=== COMPACT Analysis ===")
  message("  Genes:      ", nrow(expr))
  message("  Samples:    ", ncol(expr))
  message("  Conditions: ", paste(diets, collapse = ", "))
  message("  Timepoints: ", paste(unique(sample_info$timepoint), collapse = ", "))
  message("  FC threshold (log2): ", round(fc_threshold, 3))

  message("\n[1/4] Calculating group averages...")
  avg_result <- .calc_group_averages(expr, sample_info)

  message("[2/4] Discretizing expression patterns...")
  pattern_data <- .discretize_patterns(
    avg_result$avg_expr, avg_result$group_labels, fc_threshold
  )

  message("[3/4] Building pattern strings per condition...")
  pattern_strings <- .build_pattern_strings(pattern_data, diets)
  pattern_counts  <- lapply(pattern_strings, table)

  message("[4/4] Computing pairwise COMPACT matrices...")
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' is required. Install with: install.packages('gtools')")
  }
  diet_pairs       <- gtools::combinations(length(diets), 2, diets)
  compact_matrices <- list()

  for (i in seq_len(nrow(diet_pairs))) {
    d1  <- diet_pairs[i, 1]
    d2  <- diet_pairs[i, 2]
    key <- paste0(d1, ":", d2)
    compact_matrices[[key]] <- .build_compact_matrix(pattern_strings, d1, d2, symmetric)
    message("  ", key, " matrix: ",
            nrow(compact_matrices[[key]]), " x ", ncol(compact_matrices[[key]]))
  }

  message("\nDone. Use compact_plot(results) to visualize.")

  return(list(
    compact_matrices = compact_matrices,
    pattern_strings  = pattern_strings,
    pattern_counts   = pattern_counts,
    avg_expr         = avg_result$avg_expr,
    discrete         = pattern_data$discrete
  ))
}
