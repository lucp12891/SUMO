#' SUMO: Simulation Utilities for Multi-Omics Data
#'
#' It provides tools for simulating complex multi-omics datasets, enabling researchers to generate data that mirrors the biological intricacies observed in real-world omics studies. This package addresses a critical gap in current bioinformatics by offering flexible and customizable methods for synthetic multi-omics data generation, supporting method development, validation, and benchmarking.
#'
#' **Key Features:**
#' - **Multi-Omics Simulation**: Generate multi-layered datasets with shared and modality-specific structures.
#' - **Flexible Generation Engine**: Fine control over samples, noise levels, signal distributions, and latent factor structures.
#' - **Pipeline-Friendly Design**: Seamlessly integrates with existing multi-omics analysis workflows and packages (e.g., `SummarizedExperiment`, `MultiAssayExperiment`).
#' - **Visualization Tools**: Built-in plotting functions for exploring synthetic signals, factor structures, and noise.
#'
#' **Main Functions:**
#' - `OmixCraftHD()`: Simulates synthetic high-dimensional multi-omics datasets.
#' - `plot_simData()`: Visualizes generated data at different levels.
#' - `plot_factor()`: Displays factor scores across samples for signal inspection.
#' - `plot_weights()`: Visualizes feature loadings to assess signal versus noise.
#' - `demo_multiomics_analysis()`: Full demo function for applying MOFA on SUMO-generated or real-world data.
#' - `compute_means_vars()`: Estimate parameters from the real experimental dataset.
#'
#' @docType package
#' @name SUMO
#' @keywords benchmarking multi-omics models
"_PACKAGE"
