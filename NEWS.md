## SUMO 0.2.0

### 🔧 Major Enhancements
- Introduced `demo_multiomics_analysis()` for end-to-end demonstration of MOFA-based integrative analysis using both SUMO-generated and real-world datasets (e.g., CLL).
- Added PowerPoint reporting with automatic visualization export using `officer` and `rvg`.

### 🆕 New Functions
- `compute_means_vars()`: Computes sample- and feature-level statistics for benchmarking and noise calibration across multiple omics datasets.
- `plot_weights()`: Refactored with modular plotting helpers and new options for:
  - `show.legend` (toggle legend)
  - `"integrated"` mode (combined view of both omics)
  - Signal region annotation from `sim_object$signal_annotation`
- `plot_factor()`: Now supports both scatter and histogram visualization of factor scores. Includes support for plotting all factors and controlling legend display.

### 🎨 Visualization Improvements
- Modularized plotting into helper functions (`build_scatter_plot()`, `build_histogram_plot()`).
- Added consistent theming, Viridis color scales, and better axis labeling across plots.
- Improved annotation of signal-vs-noise regions for clarity in benchmarking.

### 🧪 Simulation Engine Upgrades
- Updated `OmixCraftHD()` to support:
  - Per-factor specification of means and variances for samples and features (e.g., `signal.samples = c(3, 0.5)`)
  - Automatic derivation of signal masks used for evaluation and annotation
- Better error handling and defaults for edge cases (e.g., NULL factors)

### 🧼 Internal Cleanup
- Removed deprecated argument `signal_vector`
- Added `globalVariables()` suppressions for clean CRAN checks
- Improved code readability and consistency across plotting functions

### 📦 Infrastructure
- Package now includes robust examples for every major exported function
- Improved documentation headers for `Roxygen2` and pkgdown compatibility
