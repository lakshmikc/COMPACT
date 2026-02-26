# =============================================================================
# COMPACT Demo Script
# examples/demo.R
#
# End-to-end example of COMPACT analysis using simulated data.
# Run this script after installing COMPACT to verify everything works.
#
# Author: Lakshmi Kuttippurathu, PhD
# =============================================================================

# --- 0. Load functions (if not installed as a package) -----------------------
source("R/compact_run.R")
source("R/compact_plot.R")

# Install dependencies if needed:
# install.packages(c("gtools", "ggplot2", "reshape2"))


# --- 1. Simulate example data ------------------------------------------------
# Generates a synthetic dataset:
#   - 200 genes
#   - 3 dietary conditions: CHO (chow), HF (high-fat), EtOH (ethanol)
#   - 3 timepoints: 1h, 6h, 24h post treatment
#   - 3 biological replicates per group

cat("\n--- Step 1: Simulate data ---\n")
example_data <- compact_simulate(n_genes = 200, n_samples = 3, seed = 42)

# Inspect the data structure
cat("Expression matrix dimensions:", dim(example_data$expr), "\n")
cat("Sample info preview:\n")
print(head(example_data$sample_info, 9))


# --- 2. Run COMPACT analysis -------------------------------------------------
# This will:
#   - Average replicates within each group
#   - Discretize into -1 / 0 / +1 patterns using log2(1.5) threshold
#   - Build pattern strings per dietary condition
#   - Compute all pairwise COMPACT matrices (CHO:HF, CHO:EtOH, HF:EtOH)

cat("\n--- Step 2: Run COMPACT ---\n")
results <- compact_run(
  expr         = example_data$expr,
  sample_info  = example_data$sample_info,
  fc_threshold = log2(1.5),
  symmetric    = FALSE
)


# --- 3. Explore results ------------------------------------------------------

cat("\n--- Step 3: Explore results ---\n")

# View available condition pairs
cat("Condition pairs analyzed:\n")
print(names(results$compact_matrices))

# View top patterns in CHO condition
cat("\nTop patterns in CHO condition:\n")
print(sort(results$pattern_counts[["CHO"]], decreasing = TRUE)[1:10])

# View CHO vs HF COMPACT matrix
cat("\nCHO vs HF COMPACT matrix:\n")
print(results$compact_matrices[["CHO:HF"]])

# Find genes with a specific pattern in CHO (e.g., up at all timepoints: "1_1_1")
target_pattern <- "1_1_1"
genes_with_pattern <- names(which(results$pattern_strings[["CHO"]] == target_pattern))
cat("\nGenes with pattern '", target_pattern, "' in CHO (", length(genes_with_pattern), " genes):\n", sep = "")
print(head(genes_with_pattern, 10))


# --- 4. Visualize results ----------------------------------------------------

cat("\n--- Step 4: Visualize ---\n")

# Heatmap of CHO vs HF COMPACT matrix
compact_plot(results, type = "heatmap", pair = "CHO:HF")

# Top 15 patterns in the EtOH condition
compact_plot(results, type = "barplot", condition = "EtOH", top_n = 15)

# Both plots for CHO vs EtOH
compact_plot(results, type = "both", pair = "CHO:EtOH", condition = "CHO")

# Save all plots to PDF
compact_plot(results, type = "both", save_pdf = TRUE, pdf_file = "compact_demo_plots.pdf")

cat("\nDemo complete! Check compact_demo_plots.pdf for saved figures.\n")

