# Neutrophil Morphokinetic Reference Guide: Colors, Parameters & Behavioral Signatures

This reference guide details the locked visual palettes, the 34 extracted morphokinetic parameters, and the characteristic feature signatures defining each behavioral state (B1–B5).

---

## 1. Behavioral Classification & Color Palette

Each neutrophil behavior is locked to a definitive biological designation and a fixed hex color code.

```
[B1: Yellow]         [B2: Green]          [B3: Purple]         [B4: Pink]           [B5: Blue]
Sessile / Spherical  Confined / Erratic   Directional / Migratory  Flattened / Dynamic  Adherent / Oscillatory
#E8E654              #5EC468              #C780CF              #DB5390              #5DACF5
```

| ID | Definitive Label | Hex Code | Color | Physical Definition |
| :--- | :--- | :--- | :--- | :--- |
| **B1** | `B1: Sessile/Spherical` | `#E8E654` | Yellow | Resting, unpolarized, spherical cells with minimal displacement. |
| **B2** | `B2: Confined/Erratic` | `#5EC468` | Green | Localized tumbling and searching within confined 3D microenvironments. |
| **B3** | `B3: Directional/Migratory` | `#C780CF` | Purple | Hyper-kinetic, elongated front-to-back polarity with rapid directional translocation. |
| **B4** | `B4: Flattened/Dynamic` | `#DB5390` | Pink | Dynamic deformation, oblate flattening, and crawling along surfaces. |
| **B5** | `B5: Adherent/Oscillatory` | `#5DACF5` | Blue | Volumetric expansion (swelling), firm adhesion, high contact area, and localized membrane oscillation. |

---

## 2. Behavioral Morphokinetic Signatures

### B1: Sessile / Spherical
* **Color:** `#E8E654` (Yellow)
* **Enriched Metrics:** `MeanSphericity`, `Tortuosity`, `AngleOfDirection_Std`
* **Depleted Metrics:** `MinVelocity`, `MeanProlateEllipt`, `TraveledDistance`
* **Physical Characteristics:** Near-perfect spherical geometry (Sphericity Ψ ≥ 0.85); minimal net translocation; high tortuosity and angular variance driven by localized Brownian motion rather than active migration.

---

### B2: Confined / Erratic
* **Color:** `#5EC468` (Green)
* **Enriched Metrics:** `StdAR3` (Aspect Ratio 3 variance), `StdProlateEllipt`
* **Depleted Metrics:** `StdAR1`, `MeanMidAxis`
* **Physical Characteristics:** Frequent directional shifts and tumbling; moderate displacement; dynamic fluctuations in longitudinal aspect ratios as cells navigate confined matrices.

---

### B3: Directional / Migratory
* **Color:** `#C780CF` (Purple)
* **Enriched Metrics:** `LinearDistance`, `TraveledDistance`, `MaxVelocity`, `MeanVelocity`
* **Depleted Metrics:** `StdOblateEllipt`, `StdProlateEllipt`
* **Physical Characteristics:** Maximized cumulative path length and linear displacement; highest peak and average velocities; persistent front-to-back polarized morphology with low shape variance.

---

### B4: Flattened / Dynamic
* **Color:** `#DB5390` (Pink)
* **Enriched Metrics:** `MeanAR3`, `MeanOblateEllipt`, `StdAR2`, `MeanAR2`, `StdMajorAxis`
* **Depleted Metrics:** `MeanMinAxis`, `CVvelocity`
* **Physical Characteristics:** Flattened oblate morphology ("pancake" profile); high aspect-ratio deformability; significant elongation and retraction along the major axis during active substrate crawling.

---

### B5: Adherent / Oscillatory *(Pathogenic)*
* **Color:** `#5DACF5` (Blue)
* **Enriched Metrics:** `MeanVolume`, `StdVolume`, `MeanSurfaceArea`, `StdSurfaceArea`, `StdMajorAxis`
* **Depleted Metrics:** `MeanSphericity` (Oblate distortion), `StdAR3`
* **Physical Characteristics:** Maximized 3D volume and surface area; marked reduction in sphericity; near-zero net migration (sessile) accompanied by localized, high-frequency membrane oscillations.

---

## 3. Morphokinetic Parameter Dictionary (34 Features)

The quantitative feature space consists of 34 continuous metrics divided into three domains:

### 3.1. Volumetric & Surface Geometry (10 Features)

| Feature Name | Description |
| :--- | :--- |
| `MeanVolume` | Average 3D cell volume across all frames of the track (in voxels / calibrated µm³). |
| `StdVolume` | Standard deviation of volume across time (indicates swelling or compression). |
| `MeanSurfaceArea` | Average external surface area computed from 3D boundary voxels. |
| `StdSurfaceArea` | Standard deviation of external surface area across time. |
| `MeanSphericity` | Average ratio of the surface area of an equivalent sphere to the actual cell surface area (range 0 to 1). |
| `StdSphericity` | Standard deviation of sphericity across time. |
| `MeanMajorAxis` | Average length of the primary axis of the fitted equivalent bounding ellipsoid. |
| `StdMajorAxis` | Standard deviation of the primary axis length (measures dynamic stretching). |
| `MeanMidAxis` | Average length of the intermediate axis of the bounding ellipsoid. |
| `StdMidAxis` | Standard deviation of the intermediate axis length. |

---

### 3.2. Ellipticity & Aspect Ratios (8 Features)

| Feature Name | Description |
| :--- | :--- |
| `MeanOblateEllipt` | Degree of oblate flattening (pancake-like flattening index). |
| `StdOblateEllipt` | Temporal standard deviation of oblate ellipticity. |
| `MeanProlateEllipt` | Degree of prolate elongation (cigar-like elongation index). |
| `StdProlateEllipt` | Temporal standard deviation of prolate ellipticity. |
| `MeanAR1` | Primary aspect ratio: Major Axis / Mid Axis. |
| `StdAR1` | Temporal standard deviation of Aspect Ratio 1. |
| `MeanAR2` | Secondary aspect ratio: Major Axis / Minor Axis. |
| `StdAR2` | Temporal standard deviation of Aspect Ratio 2. |
| `MeanAR3` | Tertiary aspect ratio: Mid Axis / Minor Axis. |
| `StdAR3` | Temporal standard deviation of Aspect Ratio 3. |

---

### 3.3. Kinematics & Spatial Trajectory (16 Features)

| Feature Name | Description |
| :--- | :--- |
| `MeanVelocity` | Average instantaneous centroid speed per unit time (µm/min). |
| `StdVelocity` | Standard deviation of instantaneous speed across the track. |
| `MaxVelocity` | Maximum instantaneous speed observed along the trajectory. |
| `MinVelocity` | Minimum instantaneous speed observed along the trajectory. |
| `CVvelocity` | Coefficient of variation of velocity (StdVelocity / MeanVelocity). |
| `LinearDistance` | Euclidean straight-line distance from track start to track end (µm). |
| `TraveledDistance` | Total cumulative path length integrated across all consecutive timepoints (µm). |
| `Tortuosity` | Ratio of total path length to net linear distance (TraveledDistance / LinearDistance). |
| `MeanDisplacement` | Mean frame-to-frame centroid displacement. |
| `StdDisplacement` | Standard deviation of frame-to-frame centroid displacement. |
| `AngleOfDirection_Mean` | Mean direction of travel relative to the reference coordinate plane (radians). |
| `AngleOfDirection_Std` | Standard deviation of angular direction (measures directional persistence). |
| `MeanDistanceToSurface` | Mean vertical distance between the cell centroid and the coverslip/substrate plane (µm). |
| `StdDistanceToSurface` | Standard deviation of vertical distance to the substrate. |

---

## 4. Feature Summary Across Behaviors

| Feature Group | B1 (Yellow) | B2 (Green) | B3 (Purple) | B4 (Pink) | B5 (Blue) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Volume (`MeanVolume`)** | Low | Moderate | Moderate | Moderate | **Highest** |
| **Surface Area (`MeanSurfaceArea`)** | Low | Moderate | Moderate | Moderate | **Highest** |
| **Sphericity (`MeanSphericity`)** | **Highest** | Moderate | Moderate | Low | **Lowest** |
| **Velocity (`MeanVelocity`)** | **Lowest** | Moderate | **Highest** | Moderate | Low |
| **Linear Distance (`LinearDistance`)** | Low | Low | **Highest** | Moderate | Low |
| **Oblate Ellipticity (`MeanOblateEllipt`)** | Low | Low | Low | **Highest** | High |
| **Prolate Ellipticity (`StdProlateEllipt`)** | Low | **High** | Low | Moderate | Low |
| **Major Axis Variance (`StdMajorAxis`)** | Low | Moderate | Moderate | **High** | **High** |
| **Directional Persistence (`AngleOfDirection_Std`)** | Low (High variance) | Moderate | **Highest** (Low variance) | Moderate | Low (High variance) |
