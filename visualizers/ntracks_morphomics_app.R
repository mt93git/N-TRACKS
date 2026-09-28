# =========================================================================
# N-TRACKS: Morphomics Dashboard Launcher
# =========================================================================
# Convenience launcher redirecting to src/r_integration/ntracks_morphomics_app.R

script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
target_file <- file.path(dirname(script_dir), "src", "r_integration", "ntracks_morphomics_app.R")

if (!file.exists(target_file)) {
  target_file <- "src/r_integration/ntracks_morphomics_app.R"
}

source(target_file)
