# Installation Guide — COMPACT

## System Requirements

| Requirement | Minimum Version |
|-------------|----------------|
| R           | 4.0.0          |
| RStudio     | 1.4+ (optional but recommended) |
| RAM         | 8 GB (16 GB recommended for large datasets) |
| OS          | Linux, macOS, Windows |

---

## Step 1 — Install R

Download and install R from [https://cran.r-project.org](https://cran.r-project.org).

---

## Step 2 — Install Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install()
```

---

## Step 3 — Install Dependencies

```r
# CRAN packages
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "reshape2",
  "pheatmap",
  "RColorBrewer",
  "devtools"
))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "SummarizedExperiment"
))
```

---

## Step 4 — Install COMPACT

```r
devtools::install_github("lakshmikc/COMPACT")
```

---

## Step 5 — Verify Installation

```r
library(COMPACT)
packageVersion("COMPACT")

# Run a quick test with built-in example data
data(compact_example)
results <- compact_run(expr = compact_example$expr, groups = compact_example$groups)
print(results)
```

---

## Troubleshooting

**Issue:** `devtools::install_github` fails
**Fix:** Make sure you have a working internet connection and that the `devtools` package is up to date:
```r
update.packages("devtools")
```

**Issue:** Bioconductor package installation errors
**Fix:** Make sure your Bioconductor version matches your R version:
```r
BiocManager::valid()
```

**Issue:** Memory errors on large datasets
**Fix:** Increase R's memory limit or subset your data before running:
```r
options(java.parameters = "-Xmx8g")  # Increase to 8GB
```

---

## Getting Help

- 📖 See the [README](README.md) for usage examples
- 🐛 Report bugs via [GitHub Issues](https://github.com/lakshmikc/COMPACT/issues)
- 📬 Contact: lakshmikc@gmail.com
