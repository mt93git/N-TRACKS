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

# 4. Download Bio-Formats Java Bridge (Open Microscopy Environment - GNU GPL v2+)
lib_dir <- file.path(getwd(), "lib")
if (!dir.exists(lib_dir)) dir.create(lib_dir, recursive = TRUE)
jar_path <- file.path(lib_dir, "bioformats_package.jar")

if (!file.exists(jar_path)) {
  message("--> Downloading Bio-Formats Java Bridge from Open Microscopy Environment (GPL v2+)...")
  tryCatch({
    download.file(
      url = "https://downloads.openmicroscopy.org/bio-formats/7.0.0/artifacts/bioformats_package.jar",
      destfile = jar_path,
      mode = "wb"
    )
    message("--> Bio-Formats successfully installed to lib/bioformats_package.jar")
  }, error = function(e) {
    warning("Automated download failed. Please manually download bioformats_package.jar from https://www.openmicroscopy.org/bio-formats/ and place it in the lib/ folder.")
  })
} else {
  message("--> Bio-Formats Java Bridge already present in lib/.")
}

message("==> N-TRACKS Environment Ready. Launch the profiler with shiny::runApp('src/r_integration/ntracks_morphometric_profiler.R')")
