# N-TRACKS (v2.0)

**N**eutrophil **T**racking & **R**ecognition **A**nalysis of **C**ell **K**inetics & **S**hape

## 1. Overview
**N-TRACKS** is a hybrid computational pipeline designed to bridge the gap between high-content 4D imaging and single-cell transcriptomics. Built to process >10TB datasets, it integrates a robust MATLAB segmentation engine with an R-based multi-omics harmonization layer.

**The Mission:**
To move beyond simple cell counting and establish a direct link between *morphological behavior* (Motility, Sphericity) and *gene expression programs* (scRNA-seq clusters).

## 2. Architecture

### ?? The Engine (MATLAB)
Located in src/matlab/:
* **
tracks_core_engine.m**: The primary algorithm. Uses modified watershed segmentation to detect low-contrast neutrophil boundaries in noisy biological environments.
* **un_batch_pipeline.m**: An HPC-optimized wrapper allowing parallel processing of multiple experimental conditions.
* **
tracks_qc_gui.m**: A "Human-in-the-Loop" Quality Control interface for verifying segmentation masks before downstream analysis.

### ?? The Bridge (R)
Located in src/r_integration/:
* **
tracks_transcriptomic_integrator.R**: The handshake protocol. This script ingests morpho-kinetic data and maps it onto Seurat objects, enabling the discovery of "Fast-Mover" vs. "Arrested" transcriptomic states.

## 3. Usage Strategy

### Step A: Extraction (MATLAB)
Run un_batch_pipeline.m to process raw .nd2 or .tif files. The system utilizes Bio-Formats to handle massive files without memory overflow.

### Step B: Validation (GUI)
Launch 
tracks_qc_gui.m to visualize random samples of the output. Flag and remove artifacts (debris, doublets) to ensure Ground Truth quality.

### Step C: Integration (R)
Load your Seurat object and the N-TRACKS output CSV into 
tracks_transcriptomic_integrator.R. The script will normalize the kinetic features and append them as metadata columns.

## 4. Credits & Attribution
**N-TRACKS** is an evolutionary fork of the **ACME** framework.

* **Original Algorithm:** Developed by [Miguel Palomino et al.](https://github.com/miguel55/ACME).
* **v2.0 Integration & Architecture:** Developed by **Maxence Tricaud (Yale University)**.
    * *Contribution:* Refactoring for HPC batch processing, development of the QC GUI, and creation of the R-based multi-omics integration layer.

## 5. Requirements
* **MATLAB:** R2021b+ (Image Processing Toolbox)
* **R:** v4.0+ (Seurat, Tidyverse)
* **Dependencies:** Bio-Formats JAR (included in /lib)
