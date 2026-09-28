# N-TRACKS (v2.0)

[![MATLAB Version](https://img.shields.io/badge/MATLAB-R2021b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Bio-Formats](https://img.shields.io/badge/Bio--Formats-Integrated-green.svg)](https://www.openmicroscopy.org/bio-formats/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Neutrophil Tracking & Recognition Analysis of Cell Kinetics & Shape**

> A hybrid computational framework engineered to bridge multi-terabyte 3D/4D confocal live-cell imaging and single-cell phenotypic profiling.  
> Developed by **Maxence Tricaud**.

---

## 1. Overview & Biological Mission

**N-TRACKS** provides an end-to-end computational pipeline for high-throughput tracking and morphometric phenotyping of motile and resident leukocytes (neutrophils, macrophages) in noisy 3D confocal time-series (`.nd2`, `.tif`).

Beyond standard centroid tracking, **N-TRACKS** quantifies granular 3D morphological dynamics (sphericity variance, prolate/oblate aspect ratios, volume stability, directional tortuosity) to link spatial cellular behavior to downstream phenotypic programs.

---

## 2. Pipeline Architecture

```
N-TRACKS/
├── src/
│   ├── matlab/
│   │   ├── core/
│   │   │   ├── ntracks_core_engine.m     # 3D watershed segmentation & Hungarian LAP assignment
│   │   │   └── ntracks_stitcher.m        # Multi-tile and time-series track stitcher
│   │   ├── gui/
│   │   │   └── ntracks_qc_gui.m          # Human-in-the-loop quality control interface
│   │   └── workflows/
│   │       ├── run_batch_pipeline.m      # Automated headless batch processor
│   │       └── run_qc_pipeline.m         # Automated QC verification harness
│   └── r_integration/
│       ├── ntracks_morphometric_profiler.R # 2200-line analytical Shiny dashboard
│       └── ntracks_morphomics_app.R        # Modular phenotypic profiler
├── data/
│   └── output_examples/
│       └── MASTER_Stitched_Harmonized_Sepsis_ACME.csv # Curated benchmark track dataset
├── lib/
│   └── (bioformats_package.jar)         # Dynamic Java bridge (auto-downloaded during setup)
├── install_ntracks_env.R                # Automated R dependency installer & Bio-Formats fetcher
├── LICENSE                              # MIT Open-Source License
└── README.md
```

---

## 3. Algorithmic Principles

### A. 3D Boundary Extraction (Adaptive Watershed)
Located in `src/matlab/core/ntracks_core_engine.m`:
* Ingests multi-slice 3D Z-stacks via a dynamically loaded local **Bio-Formats** JVM bridge.
* Applies multi-scale intensity thresholding with morphological gradient filtering and 3D distance-transform watershed to resolve overlapping cell boundaries in low signal-to-noise acquisitions.
* Filters out imaging debris using physical volume constraints (`min_volume`, `max_volume`) and 3D sphericity cutoffs (`min_sphericity`).

### B. Spatiotemporal Kinetic Tracking (Hungarian LAP)
* Solves the global **Linear Assignment Problem (LAP)** across contiguous timeframes via the Kuhn-Munkres algorithm.
* Handles cell entry/exit, transient occlusion, and motility shifts using a parameterized cost matrix (`costOfNonAssign = 20`).

### C. Downstream Phenotypic Profiling (R / Shiny)
Located in `src/r_integration/ntracks_morphometric_profiler.R`:
* Integrates curated kinetic trajectories (`MASTER_Stitched_Harmonized_Sepsis_ACME.csv`) into a high-dimensional exploratory dashboard.
* Extracts multi-condition velocity distributions, migration tortuosity, and elongation dynamics across experimental cohorts.

---

## 4. Quick Start Guide

### Step 1: Install R Dependencies
```R
Rscript install_ntracks_env.R
```

### Step 2: Run Headless MATLAB Batch Processing
Open MATLAB and execute:
```matlab
addpath(genpath('src/matlab'));
run_batch_pipeline();
```

### Step 3: Launch Interactive Profiler
In R / RStudio:
```R
shiny::runApp("src/r_integration/ntracks_morphometric_profiler.R")
```

---

## 5. Credits & Attribution

* **v2.0 Integration & Multi-Omics Architecture:** Maxence Tricaud.
* **Foundational Tracking Engine:** Derived and evolved from the ACME core framework ([Palomino et al.](https://github.com/miguel55/ACME)), substantially refactored for automated HPC batch execution, dynamic Bio-Formats self-loading, interactive QC GUI, and R phenotypic harmonization.

---

## 6. License & Open-Source Attributions

* **N-TRACKS Core Framework:** Distributed under the **MIT License**. Copyright (c) 2025-2026 Maxence Tricaud. See `LICENSE` for full details.
* **Third-Party Dependency (Bio-Formats):** Uses the Bio-Formats Java library, developed by the Open Microscopy Environment (OME) and distributed under the **GNU General Public License (GPL) v2 or later**. Bio-Formats is not tracked or bundled in this repository; it is dynamically downloaded to `lib/` during environment setup.
* **Foundational Tracking Algorithm:** Tracking core derived from the ACME framework ([Palomino et al.](https://github.com/miguel55/ACME)).
