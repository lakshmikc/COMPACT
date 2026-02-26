# =============================================================================
# COMPACT — Comparative Pattern Counts
# compact_plot.R
#
# Visualization functions for COMPACT analysis results.
#
# Author: Lakshmi Kuttippurathu, PhD
# License: MIT
# Repository: https://github.com/lakshmikc/COMPACT
# =============================================================================


#' Plot COMPACT Analysis Results
#'
#' Generates visualizations of COMPACT output. Supports three plot types:
#'
#' \describe{
#'   \item{\code{"heatmap"}}{A heatmap of the COMPACT count matrix for a chosen
#'     condition pair, showing how many genes share each cross-condition pattern.
#'     Cells with zero counts are shown in white; higher counts in red.}
#'   \item{\code{"barplot"}}{A bar chart of the most frequent expression patterns
#'     within a single condition, useful for identifying dominant gene behaviors.}
#'   \item{\code{"both"}}{Renders both plots side by side (default).}
#' }
#'
#' @param results    Output list from \code{\link{compact_run}}.
#' @param type       Character. One of \code{"heatmap"}, \code{"barplot"}, or
#'                   \code{"both"} (default: \code{"both"}).
#' @param pair       Character. Which condition pair to plot for the heatmap,
#'                   e.g. \code{"CHO:HF"}. Defaults to the first available pair.
#' @param condition  Character. Which condition to plot for the barplot,
#'                   e.g. \code{"CHO"}. Defaults to the first available condition.
#' @param top_n      Integer. Number of top patterns to show in the barplot
#'                   (default: 15).
#' @param save_pdf   Logical. If \code{TRUE}, saves plots to a PDF file
#'                   (default: \code{FALSE}).
#' @param pdf_file   Character. File path for PDF output if \code{save_pdf = TRUE}
#'                   (default: \code{"compact_plots.pdf"}).
#'
#' @return Invisibly returns a list of ggplot2 objects. Plots are also printed
#'   to the active graphics device.
#'
#' @examples
#' example_data <- compact_simulate()
#' results      <- compact_run(example_data$expr, example_data$sample_info)
#'
#' # Plot both (default)
#' compact_plot(results)
#'
#' # Heatmap only for a specific pair
#' compact_plot(results, type = "heatmap", pair = "CHO:HF")
#'
#' # Top pattern barplot for EtOH condition
#' compact_plot(results, type = "barplot", condition = "EtOH", top_n = 10)
#'
#' # Save to PDF
#' compact_plot(results, save_pdf = TRUE, pdf_file = "my_compact_plots.pdf")
#'
#' @seealso \code{\link{compact_run}}, \code{\link{compact_simulate}}
#' @export
compact_plot <- function(results,
                         type      = "both",
                         pair      = NULL,
                         condition = NULL,
                         top_n     = 15,
                         save_pdf  = FALSE,
                         pdf_file  = "compact_plots.pdf") {

  # --- Validate input ---------------------------------------------------------
  required <- c("compact_matrices", "pattern_counts")
  missing  <- setdiff(required, names(results))
  if (length(missing) > 0) {
    stop("Input does not look like compact_run() output. Missing: ",
         paste(missing, collapse = ", "))
  }
  if (!type %in% c("heatmap", "barplot", "both")) {
    stop("'type' must be one of: 'heatmap', 'barplot', 'both'")
  }
  if (!requireNamespace("ggplot2",  quietly = TRUE)) stop("Install ggplot2:  install.packages('ggplot2')")
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Install reshape2: install.packages('reshape2')")

  # Default pair and condition
  if (is.null(pair))      pair      <- names(results$compact_matrices)[1]
  if (is.null(condition)) condition <- names(results$pattern_counts)[1]

  plots <- list()

  # --- Heatmap ----------------------------------------------------------------
  if (type %in% c("heatmap", "both")) {

    if (!pair %in% names(results$compact_matrices)) {
      stop("Pair '", pair, "' not found. Available: ",
           paste(names(results$compact_matrices), collapse = ", "))
    }

    mat      <- results$compact_matrices[[pair]]
    mat_long <- reshape2::melt(mat, varnames = c("Pattern1", "Pattern2"), value.name = "Count")

    p_heatmap <- ggplot2::ggplot(mat_long,
        ggplot2::aes(x = Pattern2, y = Pattern1, fill = Count)) +
      ggplot2::geom_tile(color = "grey90", linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(Count > 0, Count, "")),
        size = 3, color = "black"
      ) +
      ggplot2::scale_fill_gradient(
        low      = "white",
        high     = "#CC0000",
        na.value = "white",
        name     = "Gene\nCount"
      ) +
      ggplot2::labs(
        title    = paste0("COMPACT Matrix: ", pair),
        subtitle = "Each cell = number of genes sharing that cross-condition pattern",
        x        = gsub(".*:", "", pair),
        y        = gsub(":.*", "", pair)
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y      = ggplot2::element_text(size = 8),
        plot.title       = ggplot2::element_text(face = "bold"),
        panel.grid       = ggplot2::element_blank()
      )

    plots$heatmap <- p_heatmap
    print(p_heatmap)
  }

  # --- Barplot ----------------------------------------------------------------
  if (type %in% c("barplot", "both")) {

    if (!condition %in% names(results$pattern_counts)) {
      stop("Condition '", condition, "' not found. Available: ",
           paste(names(results$pattern_counts), collapse = ", "))
    }

    cnt_table <- results$pattern_counts[[condition]]
    cnt_df    <- as.data.frame(cnt_table, stringsAsFactors = FALSE)
    colnames(cnt_df) <- c("Pattern", "Count")
    cnt_df    <- cnt_df[order(cnt_df$Count, decreasing = TRUE), ]
    cnt_df    <- head(cnt_df, top_n)
    cnt_df$Pattern <- factor(cnt_df$Pattern, levels = rev(cnt_df$Pattern))

    p_barplot <- ggplot2::ggplot(cnt_df,
        ggplot2::aes(x = Pattern, y = Count, fill = Count)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::scale_fill_gradient(low = "#FDDBC7", high = "#B2182B") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title    = paste0("Top ", top_n, " Expression Patterns: ", condition),
        subtitle = "Pattern format: timepoint1_timepoint2_timepoint3 | -1=down, 0=unchanged, 1=up",
        x        = "Pattern",
        y        = "Number of Genes"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank()
      )

    plots$barplot <- p_barplot
    print(p_barplot)
  }

  # --- Save to PDF ------------------------------------------------------------
  if (save_pdf) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for PDF saving.")
    pdf(pdf_file, width = 10, height = 7)
    for (p in plots) print(p)
    dev.off()
    message("Plots saved to: ", pdf_file)
  }

  invisible(plots)
}
