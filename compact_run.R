#' Run COMPACT Analysis
#'
#' Main function to run Comparative Pattern Counts analysis on
#' a gene expression matrix across biological conditions.
#'
#' @param expr A numeric matrix or data frame (genes x samples)
#' @param groups A character vector of group labels for each sample
#' @param threshold Adjusted p-value cutoff (default: 0.05)
#' @param pattern Pattern type to analyze: "all", "shared", or "unique" (default: "all")
#' @param normalization Normalization method: "DESeq2", "TMM", or "none" (default: "DESeq2")
#'
#' @return A list containing:
#'   \item{patterns}{Data frame of identified expression patterns}
#'   \item{counts}{Pattern count summary per condition}
#'   \item{stats}{Statistical results for each gene}
#'
#' @examples
#' data(compact_example)
#' results <- compact_run(
#'   expr   = compact_example$expr,
#'   groups = compact_example$groups
#' )
#' head(results$patterns)
#'
#' @export
compact_run <- function(expr,
                        groups,
                        threshold     = 0.05,
                        pattern       = "all",
                        normalization = "DESeq2") {

  # --- Input validation ---
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("'expr' must be a matrix or data frame.")
  }
  if (ncol(expr) != length(groups)) {
    stop("Number of columns in 'expr' must match length of 'groups'.")
  }

  message("Running COMPACT analysis...")
  message("  Samples: ", ncol(expr))
  message("  Genes:   ", nrow(expr))
  message("  Groups:  ", paste(unique(groups), collapse = ", "))

  # --- Normalization ---
  expr_norm <- .normalize(expr, method = normalization)

  # --- Pattern detection ---
  patterns <- .detect_patterns(expr_norm, groups, threshold, pattern)

  # --- Pattern counting ---
  counts <- .count_patterns(patterns, groups)

  message("COMPACT analysis complete. ", nrow(patterns), " patterns identified.")

  return(list(
    patterns = patterns,
    counts   = counts,
    stats    = attr(patterns, "stats")
  ))
}


#' Plot COMPACT Results
#'
#' Generate summary visualizations of COMPACT pattern analysis results.
#'
#' @param results Output from \code{compact_run()}
#' @param type Plot type: "heatmap", "barplot", or "both" (default: "both")
#' @param top_n Number of top patterns to display (default: 50)
#'
#' @return A ggplot2 object or list of plots
#'
#' @examples
#' data(compact_example)
#' results <- compact_run(compact_example$expr, compact_example$groups)
#' compact_plot(results)
#'
#' @export
compact_plot <- function(results, type = "both", top_n = 50) {

  if (!all(c("patterns", "counts") %in% names(results))) {
    stop("Input must be output from compact_run().")
  }

  plots <- list()

  if (type %in% c("heatmap", "both")) {
    plots$heatmap <- .plot_heatmap(results$patterns, top_n)
  }

  if (type %in% c("barplot", "both")) {
    plots$barplot <- .plot_counts(results$counts)
  }

  if (length(plots) == 1) return(plots[[1]])
  return(plots)
}


# ---- Internal helper functions (not exported) ----

.normalize <- function(expr, method) {
  # Placeholder — replace with your normalization logic
  if (method == "none") return(as.matrix(expr))
  message("  Normalizing with method: ", method)
  # Add DESeq2 / TMM normalization here
  return(as.matrix(expr))
}

.detect_patterns <- function(expr_norm, groups, threshold, pattern) {
  # Placeholder — replace with your pattern detection logic
  message("  Detecting patterns (type: ", pattern, ")...")
  data.frame(gene = rownames(expr_norm), pattern_id = NA, adj_pval = NA)
}

.count_patterns <- function(patterns, groups) {
  # Placeholder — replace with your pattern counting logic
  table(patterns$pattern_id)
}

.plot_heatmap <- function(patterns, top_n) {
  message("Generating heatmap for top ", top_n, " patterns...")
  # Add pheatmap / ggplot2 heatmap logic here
}

.plot_counts <- function(counts) {
  message("Generating pattern count barplot...")
  # Add ggplot2 barplot logic here
}
