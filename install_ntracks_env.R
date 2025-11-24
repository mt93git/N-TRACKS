# N-TRACKS Environment Installer
# Run this once to install all dependencies for the Morphometric Profiler App.

# 1. Define CRAN Packages
cran_pkgs <- c(
  "shiny", "shinythemes", "shinyjs", "dplyr", "ggplot2", 
  "readr", "DT", "pheatmap", "grid", "RColorBrewer", 
  "zip", "cowplot", "ggpubr", "plotly"
)

# 2. Define Bioconductor/Complex Packages
bioc_pkgs <- c("Seurat") 

# 3. Installation Loop
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

# Install CRAN
invisible(lapply(cran_pkgs, install_if_missing))

# Install Bioconductor / Seurat
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Seurat", quietly = TRUE)) BiocManager::install("Seurat")

message("? N-TRACKS Environment Ready. Launch the app in src/r_integration/")
