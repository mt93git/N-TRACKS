# N-TRACKS Feature Calibration & Unit Standardization Specification

## 1. Executive Summary & Calibration Scope
This specification defines the exact physical units, mathematical derivations, and empirical calibration scalars (K_VOL, K_VEL, K_SPH) applied across both Pipeline 1 (Windows) and Pipeline 2 (Mac Intel / Headless Linux) outputs. 

To prevent accidental double-scaling or erroneous metric comparisons during downstream execution on Yale's HPC cluster (ACME framework), this document details the calibration state of each individual column across all generated reference matrices.

---

## 2. Benchmark Calibration Scalars (Phase 8 Standardization)
All scalars were determined through empirical calibration runs against standardized fluorescent polystyrene microspheres (10 µm diameter) and validated multi-position acquisitions:

| Scalar | Value | Purpose | Application Point |
| :--- | :--- | :--- | :--- |
| **K_VOL** | **1.107220** | 3D Voxel Volumetric Expansion Correction | Applied inside MATLAB core segmentation engine prior to shape filtering. |
| **K_VEL** | **1.288988** | Kinetic Displacement-to-Velocity Scalar | Applied during trajectory displacement summation and velocity calculations. |
| **K_SPH** | **1.004726** | Normalized Sphericity Surface-Volume Ratio Alignment | Applied during final statistical aggregation and rolling window feature extraction. |

---

## 3. Detailed Matrix Calibration Ledger

### A. `MASTER_RESULTS_FRAME_BY_FRAME.csv` (Pipeline 2 Instantaneous Tracking Matrix)
> [!CAUTION]
> **HIGH RISK OF DOUBLE-SCALING ON HPC**: Notice the asymmetry between `Volume` and `Centroids` in this file!

| Column Name | Physical Dimension | Calibration State | HPC Processing Instruction |
| :--- | :--- | :--- | :--- |
| `Centroid_X` | µm / Voxel | **RAW (Unscaled)** | Compute step Euclidean displacements, then multiply by K_VEL (1.288988). |
| `Centroid_Y` | µm / Voxel | **RAW (Unscaled)** | Compute step Euclidean displacements, then multiply by K_VEL (1.288988). |
| `Centroid_Z` | µm / Slice | **RAW (Unscaled)** | Compute step Euclidean displacements, then multiply by K_VEL (1.288988). |
| `Volume` | µm³ | **ALREADY CALIBRATED** | **DO NOT SCALE.** Pre-multiplied by K_VOL (1.107220) in MATLAB `regionprops3`. Scaling again yields 1.226x error. |
| `SurfaceArea` | µm² | **RAW (Unscaled)** | Direct mesh surface area from isosurface marching cubes. |
| `Sphericity` | Dimensionless [0, 1] | **INTERMEDIATE** | Derived from [pi^(1/3) * (6 * Volume)^(2/3)] / SurfaceArea using calibrated volume. To match Phase 8, multiply by K_SPH (1.004726). |
| `AxisLength_1..3` | µm | **RAW (Unscaled)** | Principal axis lengths derived from 3D inertia tensor. |
| `MeanIntensity` | A.U. [0, 65535] | **RAW (Unscaled)** | Mean 16-bit fluorophore intensity within segmented 3D mask. |

---

### B. `MASTER_RESULTS_RAW.csv` & `MASTER_Stitched_Harmonized_Sepsis_ACME.csv` (Trajectory-Level Summary Matrices)
In these final aggregated datasets, all summary features are **fully calibrated and finalized**:

- **Velocity & Distance Metrics**:
  - `TraveledDistance` = Sum(Step_Displacement) * K_VEL (µm) [CALIBRATED]
  - `MeanVelocity`, `StdVelocity`, `MinVelocity`, `MaxVelocity` = Mean_v * K_VEL (µm/min) [CALIBRATED]
  - `LinearDistance` = Euclidean_Distance(Start, End) * K_VEL (µm) [CALIBRATED]
  - `Tortuosity` = TraveledDistance / LinearDistance (Dimensionless ≥ 1.0)
- **Morphological Metrics**:
  - `MeanVolume`, `StdVolume` = Pre-scaled by K_VOL (µm³) [CALIBRATED]
  - `MeanSphericity`, `StdSphericity` = Pre-scaled by K_SPH [CALIBRATED]
  - `MeanSurfaceArea`, `StdSurfaceArea` = Measured in µm²
- **Spatial & Directional Metrics**:
  - `DistanceToSurface_Mean`: µm distance from cell centroid to nearest boundary of z-stack envelope.
  - `AngleOfDirection_Sum`: Cumulative directional turning angle in degrees (°).

---

## 4. Verification Check
To verify that your HPC downstream pipeline is consuming calibrated data without double-scaling, execute:
```bash
python3 02_Pipeline_V2_Production_4D/runners/validate_all_schemas.py
```
This script asserts calibration integrity against `02_Pipeline_V2_Production_4D/configs/tracking_calibration.json`.
