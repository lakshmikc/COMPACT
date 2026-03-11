# COMPACT
COMPACT (Comparative Pattern Counts) is an open-source R-based pipeline for analyzing high-throughput omics data.
WHY is it different? Because it can identify both signifincat AND subtle patterns (which we othersie miss!!) across biological conditions.

# COMPACT — Comparative Pattern Counts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r&logoColor=white)]()
[![DOI](https://img.shields.io/badge/Published-Peer--Reviewed-brightgreen)]()


> ⚠️ **Work in Progress** — This repository is actively under development and not yet ready for production use. Code, documentation, and examples are being added incrementally. Watch or star the repo to follow along!



> **COMPACT** (Comparative Pattern Counts) is an open-source R-based pipeline for analyzing high-throughput omics data.
> WHY is it different? Because it can identify both signifincat AND subtle patterns (which we othersie miss!!) across biological conditions.

---

## 🧬 What Does COMPACT Do?

Biological datasets are complex, noisy, and often high-dimensional. COMPACT was designed to cut through that complexity by:

- Identifying **both sgnificant and subtle expression patterns** across experimental conditions
- Performing **comparative pattern counting** to reveal condition-specific and shared signatures
- Enabling **systems-level interpretation** of transcriptomic data
- Providing **publication-ready visualizations** of genomic patterns

Whether you're working with RNA-seq data from disease vs. control comparisons, time-series experiments, or multi-condition studies, COMPACT gives you a structured, reproducible framework to make sense of the patterns in your data.

---

## 📦 Installation

COMPACT is an R-based pipeline. To get started, ensure you have R (≥ 4.0) installed.

### Install from GitHub

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install COMPACT
devtools::install_github("lakshmikc/COMPACT")
```

### Dependencies

COMPACT requires the following R packages:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "reshape2", "pheatmap"))

# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR", "limma"))
```

---

## 🚀 Quick Start

```r
library(COMPACT)

# Load your expression matrix (genes x samples)
# Rows = genes, Columns = samples
expr_matrix <- read.csv("your_expression_data.csv", row.names = 1)

# Define your sample groups
groups <- c("control", "control", "treated", "treated")

# Run COMPACT analysis
results <- compact_run(
  expr      = expr_matrix,
  groups    = groups,
  threshold = 0.05,     # adjusted p-value cutoff
  pattern   = "all"     # analyze all pattern types
)

# View top patterns
head(results$patterns)

# Generate summary plot
compact_plot(results)
```

---

## 📁 Repository Structure

```
COMPACT/
├── R/                  # Core R functions
│   ├── compact_run.R   # Main analysis function
│   ├── compact_plot.R  # Visualization functions
│   └── utils.R         # Helper utilities
├── data/               # Example datasets
├── examples/           # Demo scripts
├── docs/               # Extended documentation
├── LICENSE
└── README.md
```

---

## 📖 Citation

If you use COMPACT in your research, please cite:

> **Kuttippurathu L**, et al. *A novel comparative pattern analysis approach identifies chronic alcohol mediated dysregulation of transcriptomic dynamics during liver regeneration.* [JBMC Genomics, 2016]. DOI: [[link to paper](https://pubmed.ncbi.nlm.nih.gov/27012785/)]

---

## 🤝 Contributing

Contributions are welcome! If you find a bug, have a feature request, or want to contribute code:

1. Fork the repository
2. Create a new branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m 'Add your feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request

---

## 📬 Contact

**Lakshmi Kuttippurathu, PhD**
- 🔗 [LinkedIn](https://www.linkedin.com/in/lakshmikc/)
- 📚 [Publications (ORCID)](https://orcid.org/0000-0001-6612-9040)
- 📧 lakshmikc@encurebio.com

---

*Built with ❤️ for the open-science community. If this tool helps your research, please consider starring the repository ⭐*
