# Germ-Granule-Aging

Code accompanying **“Germ Granule Remodeling Reveals Germline Age in Caenorhabditis elegans”** by Purbasha Dasgupta, Ryan J. Lee, Erica Dinneen, and Carlos G. Silva-García.

This project investigates age-associated changes in PGL-1-marked P granules and MUT-16-marked Mutator foci in pachytene germline nuclei of *C. elegans*. It includes study-specific image-analysis and aging-clock scripts, together with GGProfiler, an interactive Shiny application for granule profiling.

## Overview

- **Germ-granule remodeling:** extraction and statistical analysis of age-associated changes in granule architecture across two experimental batches.
- **GG-Clock:** an elastic-net model that predicts adult age from twelve animal-level granule features, with leave-one-animal-out cross-validation and projection of intermittent-fasting animals onto the ad libitum-trained model.
- **GGProfiler:** an interactive application for image cropping, granule segmentation, quality control, feature extraction, condition comparisons, and figure export.

GGProfiler is a general profiling application. The study-specific scripts perform the manuscript's aging-clock analyses; launching GGProfiler does not train or run GG-Clock.

## Code organization

The paths below assume the analysis scripts are placed in `GG-Clock/` and the application in `GGProfiler/`.

| File | Purpose |
| --- | --- |
| `GG-Clock/01_Germ_granule _remodeling.R` | Re-extracts germ granule features for Batches 1 and 2 animals and fits age-associated mixed-effects models, including batch-adjusted analyses, sensitivity analyses, and batch-by-age interaction tests. |
| `GG-Clock/02_Primary_Clock_Batch 1.R` | Processes Batch 1 images, aggregates features per animal, trains and cross-validates the primary clock with Batch 2 clock, projects IF animals, estimates feature contributions, and exports figures and statistics. |
| `GGProfiler/06_GGProfiler.R` | Defines the GGProfiler Shiny application. |


## Requirements

R and an environment capable of running R scripts are required. RStudio is convenient for opening the GGProfiler project and using the application locally.

| Component | Packages |
| --- | --- |
| GGProfiler | `shiny`, `ggplot2`, `DT`, and Bioconductor package `EBImage`; `zip` is optional for ZIP export. |
| Primary GG-Clock | `EBImage`, `glmnet`, `ggplot2`, `ggbeeswarm`, `patchwork`, `ragg`, and `openxlsx`. |
| Remodeling analysis | `EBImage`, `dplyr`, `tidyr`, `readr`, `lme4`, `lmerTest`, `ggplot2`, `scales`, and `ragg`. |
The analysis scripts attempt to install missing dependencies. Internet access is needed for installation. Figure-export settings use Times New Roman and a Cairo-capable PDF device; adjust the font setting if that font is unavailable.

If a project includes `renv.lock` and its activation files, open that project and restore its recorded package environment before running:

```r
install.packages("renv")  # if not already installed
renv::restore()
```

A lockfile in `GGProfiler/` describes that project; do not assume it also covers every dependency of the separate clock analyses. Record `sessionInfo()` with analysis outputs.

## Run GGProfiler

1. Download the repository and extract it locally.
2. Open `GGProfiler.Rproj`, if supplied, or set R's working directory to `GGProfiler/`.
3. Restore the project environment as above, or install the required packages:

```r
install.packages(c("shiny", "ggplot2", "DT", "zip", "BiocManager"))
BiocManager::install("EBImage", update = FALSE, ask = FALSE)
```

4. Launch the application from the `GGProfiler/` directory:

```r
app <- source("06_GGProfiler.R")$value
shiny::runApp(app)
```

Keep the `www/` folder alongside the script so the application can find its display assets.

### Image preparation and workflow

Use matched, aligned, single-plane images of the same nucleus. Separate grayscale TIFFs are recommended. Images for different markers from the same nucleus must have matching dimensions and crops. The application does not register images or automatically identify nuclear boundaries.

1. **Setup:** choose one to four channels and configure marker names, image sources, and segmentation settings.
2. **Crop nuclei:** optionally create matching nuclear crops from larger single-plane images.
3. **Conditions:** define experimental groups and upload each marker's images.
4. **QC:** inspect pairing and segmentation before extraction.
5. **Extract:** calculate and download measurements.
6. **Compare / Visualize:** compare conditions and export plots and overlays.

Files are paired by their final numeric ID before the extension. For example, `PGL1_1.tif` pairs with `MUT16_1.tif`. IDs must be unique within a condition, including across different animals; do not restart numbering for each animal within the same condition.

Downloads include per-image features, condition summaries, object measurements, processing logs, pairing manifests, statistics, figures, and an RDS analysis snapshot. The snapshot does not embed the original uploaded images.

**Statistical unit:** GGProfiler's built-in comparisons use image-level observations and do not model repeated nuclei within an animal. For experiments with several nuclei per animal, retain animal identifiers and use appropriate animal-level aggregation or hierarchical models outside the app.

## Interpretation and reproducibility

Spatial measurements are in pixels, and object areas are in square pixels. The twelve clock features describe counts, areas, eccentricity, nearest-neighbor spacing variation, PGL-1 radial dispersion, directional cross-channel association, and the PGL-1/MUT-16 area ratio. Historical variable names such as `red_size` and `green_size` refer to area.

GG-Clock was developed in the imaging and reproductive context described in the manuscript. Applying it to other imaging settings, genotypes, or reproductive contexts requires validation. Inspect QC overlays and retain the input data, analysis settings, software environment, and exact code version used for each analysis.

## Data availability

Raw microscopy images are available from the corresponding author upon request, as stated in the manuscript. Feature definitions and detailed methods are provided in Supplementary File 1; supplementary statistics are provided in Supplementary File 2. Full reproduction requires the corresponding input images or derived tables in addition to the scripts.

## Citation and contact

Please cite the accompanying manuscript:

Dasgupta, P., Lee, R. J., Dinneen, E., and Silva-García, C. G. *Germ Granule Remodeling Reveals Germline Age in Caenorhabditis elegans.* Manuscript.

Corresponding author: Carlos G. Silva-García, **gio@brown.edu**, Brown University.
