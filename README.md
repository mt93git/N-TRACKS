# N-TRACKS (v2.2)

[![Python Version](https://img.shields.io/badge/Python-%3E%3D3.10-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-%3E%3D2.0-EE4C2C.svg)](https://pytorch.org/)
[![MATLAB Version](https://img.shields.io/badge/MATLAB-R2021b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Bio-Formats](https://img.shields.io/badge/Bio--Formats-Dynamic--Bridge-green.svg)](https://www.openmicroscopy.org/bio-formats/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Neutrophil Tracking & Recognition Analysis of Cell Kinetics & Shape**

> A unified multi-language computational framework engineered to bridge multi-terabyte 3D/4D confocal live-cell imaging, deep learning kinetic forecasting (CNN-LSTM), Hungarian LAP assignment, and Markovian behavioral state transition dynamics.  
> Developed by **Maxence Tricaud**.

---

## 1. Overview & Biological Mission

**N-TRACKS** provides an end-to-end computational pipeline for high-throughput tracking and morphometric phenotyping of motile and resident leukocytes (neutrophils, macrophages) in noisy 3D confocal time-series (`.nd2`, `.tif`).

The platform solves key challenges in high-resolution intravital and confocal microscopy:
1. **Volumetric Feature Learning:** Deep 3D Convolutional Neural Networks (3D-CNN) capture volumetric shape gradients, volume fluctuations, and surface curvature.
2. **Temporal Dynamics & Motion Forecasting:** Bidirectional Long Short-Term Memory (Bi-LSTM) networks forecast next-step displacements to resolve dense cellular collisions and ambiguous crossings.
3. **Robust Combinatorial Linking:** Global Linear Assignment Problem (LAP) solved via Kuhn-Munkres optimization ensures globally optimal spatio-temporal cell identification.
4. **Stochastic State Transitions:** Discrete-time Markov chain modeling quantifies stochastic transitions between discrete phenotypic behavioral regimes (*Stationary*, *Directed-Chemotactic*, *Meandering*, *Swarming*, *Extravasating*).

---

## 2. Pipeline Architecture

```
N-TRACKS/
├── models/
│   └── cnn_lstm/                         # Deep Learning Spatiotemporal Subsystem
│       ├── __init__.py
│       ├── architecture.py               # 3D-CNN + Bi-LSTM neural network & multi-task loss
│       ├── infer_tracks.py               # Kinetic velocity forecasting & state prediction CLI
│       └── train.py                      # Training harness with multi-task supervision
├── r_markov/                             # Markovian Behavioral Transition Subsystem
│   ├── markov_transition_analysis.R      # Markov transition matrix computation & Chapman-Kolmogorov
│   └── markov_transition_cohort_analysis.R # Cross-cohort transition networks (igraph / ggplot2)
├── matlab_core/ (and src/matlab/)        # 3D Image Segmentation & Hungarian Tracking Engine
│   ├── ntracks_core_engine.m             # 3D distance watershed & Kuhn-Munkres LAP solver
│   ├── ntracks_stitcher.m                # Multi-tile and time-series stitcher
│   └── ntracks_qc_gui.m                  # Quality control curation interface
├── visualizers/                          # Exploratory Dashboards & Morphomics
│   ├── 01_RUN_APP.R                      # Unified Shiny launcher
│   ├── ntracks_morphomics_app.R          # High-dimensional morphokinetic profiler
│   └── ntracks_seurat_explorer.R         # Single-cell phenotypic manifold explorer
├── configs/                              # Calibration configs & batch manifests
├── schemas/                              # Pydantic validation schemas for metadata and tracks
├── data/output_examples/                 # Curated benchmark tracking datasets
│   └── MASTER_Stitched_Harmonized_Sepsis_ACME.csv
├── reference_outputs/transition_matrices/# Published Markov transition matrices & heatmaps
├── install_ntracks_env.R                 # Automated R dependencies & Bio-Formats JVM bridge
├── pyproject.toml                        # Python dependencies (PyTorch, Pandas, SciPy)
├── LICENSE                               # MIT Open-Source License
└── README.md
```

---

## 3. The Three Algorithmic Pillars

### Pillar 1: Deep Learning Spatiotemporal Kinetics (CNN-LSTM)
Located in `models/cnn_lstm/`:
* **3D-CNN Encoder:** Processes 3D voxel patches $(C=1, D, H, W)$ to extract dense shape, texture, and boundary representations.
* **Bidirectional LSTM:** Sequential recurrence over temporal window $t=1\dots T$ models cell momentum, directional persistence, and acceleration.
* **Multi-Task Heads:** Jointly optimizes displacement vector regression $(\Delta x, \Delta y, \Delta z)$, behavioral classification, and Markov transition logits $P_{ij}$.

### Pillar 2: 3D Hungarian Linear Assignment Problem (LAP)
Located in `matlab_core/` and `src/matlab/`:
* Ingests multi-channel 3D Z-stacks via dynamically linked Bio-Formats JVM bridge (auto-downloaded, no heavy binaries tracked in git).
* Computes adaptive 3D Euclidean distance-transform watershed segmentation.
* Global bipartite matching across time frames using the Kuhn-Munkres algorithm with parameterized non-assignment penalties to handle cell entry, exit, and temporary occlusion.

### Pillar 3: Markovian State Transition Dynamics
Located in `r_markov/`:
* Converts continuous trajectory metrics into discrete morphokinetic states.
* Computes empirical and regularized transition probability matrices $P = [p_{ij}]$, where $p_{ij} = \mathbb{P}(S_{t+1} = j \mid S_t = i)$.
* Validates Markovian property via Chapman-Kolmogorov tests and extracts stationary distribution vectors $\pi P = \pi$ across experimental cohorts.

---

## 4. Quick Start & Verification

### A. Verify CNN-LSTM Deep Learning Inference
```bash
# Run forward inference with the 3D-CNN + Bi-LSTM model:
python3 -m models.cnn_lstm.infer_tracks --seq-len 10
```

### B. Verify Markov Transition Modeling
```bash
# Execute Markov behavioral transition analysis:
Rscript r_markov/markov_transition_cohort_analysis.R
```

### C. Launch Interactive Morphomics Dashboard
```R
# In R or RStudio:
source("install_ntracks_env.R")
shiny::runApp("src/r_integration/ntracks_morphomics_app.R")
# (or via the visualizers convenience launcher):
shiny::runApp("visualizers/ntracks_morphomics_app.R")
```

---

## 5. Software Architecture & Reproducibility Statement

### Institutional Environment & Proprietary Tooling Rationale
The upstream 3D confocal volumetric ingestion, adaptive distance-transform watershed segmentation, and initial Kuhn-Munkres bipartite matching were originally developed and calibrated within an institutional **MATLAB environment (R2021b+)**. This architecture was selected due to the laboratory's specialized microscope hardware integration and GPU-accelerated 3D Image Processing Toolboxes required for processing multi-gigabyte `.nd2` live-cell confocal acquisitions.

### Open-Source Modularity & Downstream Accessibility
Recognizing the commercial licensing barrier and reproducibility friction inherent to MATLAB, the **N-TRACKS** architecture was explicitly engineered as a **decoupled, modular ecosystem**:
1. **Open Standard Data Exchange:** All segmented coordinate trajectories, frame-by-frame kinematics, and morphometric attributes are serialized into open, transparent tabular standards (`.csv`), validated against formal Pydantic schemas ([`schemas/`](schemas/)).
2. **Autonomous Deep Learning Engine (100% Open-Source Python):** The spatiotemporal **3D-CNN + Bi-LSTM** displacement and state prediction models ([`models/cnn_lstm/`](models/cnn_lstm/)) operate entirely on standard **PyTorch**, executing independently without any MATLAB dependency.
3. **Markovian Transition Suite (100% Open-Source R):** The Chapman-Kolmogorov validation, transition matrix estimation, and directed network graphs ([`r_markov/`](r_markov/)) run autonomously in standard **R**.
4. **Interactive Morphometrics Dashboard (100% Open-Source R/Shiny):** Downstream cohort visualization and feature distribution inspection operate independently in R/Shiny.

*Rather than superficially substituting an incomplete open-source wrapper over complex proprietary hardware drivers, N-TRACKS embraces transparency: commercial drivers are sequestered upstream, while all computational modeling, deep learning kinetics, and statistical transition analyses are 100% open, reproducible, and verifiable by the community.*

---

## 6. Credits & Open-Source Attributions

* **Architecture, Deep Learning & Multi-Omics Integration:** Maxence Tricaud.
* **Third-Party Dependency (Bio-Formats):** Uses the Bio-Formats Java library, developed by the Open Microscopy Environment (OME) and distributed under the **GNU General Public License (GPL) v2 or later**. Bio-Formats is not tracked or bundled in this repository; it is dynamically downloaded to `lib/` during environment setup.
* **Tracking Algorithm Heritage:** Tracking core derived from the ACME framework ([Palomino et al.](https://github.com/miguel55/ACME)), substantially refactored for automated HPC execution, dynamic Bio-Formats self-loading, interactive QC GUI, and R phenotypic harmonization.

---

## 7. License

Distributed under the **MIT License**. Copyright (c) 2025-2026 Maxence Tricaud. See `LICENSE` for details.
