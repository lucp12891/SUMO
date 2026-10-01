# SUMO

**Simulation Utilities for Multi-Omics Data**

SUMO is an R package for generating synthetic multi-omics datasets with known latent factors, sample signal regions, and feature loadings. It supports method development, benchmarking of integrative clustering and factor models, and teaching.

Each omics layer has its own feature space, while samples are shared across layers. Signals can affect all layers, individual layers, or subsets of layers.

## Installation

Requires **R 4.2 or later**. Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("lucp12891/SUMO", dependencies = NA)
library(SUMO)
```

This installs the required dependencies. Python is not needed for data simulation or the basic plotting functions.

## Quick start

Generate three omics layers for 100 samples, with three latent factors and different feature counts:

```r
library(SUMO)

sim <- simulateMultiOmics(
  vector_features = c(3000, 2500, 2000),
  n_samples = 100,
  n_factors = 3,
  snr = 3,
  signal.samples = c(5, 1),
  signal.features = list(
    c(3, 0.05),
    c(2.5, 0.05),
    c(2, 0.05)
  ),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)

# Dimensions of each layer: samples x features
lapply(sim$omics, dim)

# Inspect which layers contain each factor
sim$factor_map

# Display all layers together
plot_simData(sim, data = "merged")

# Display one layer with reproducible sample and feature permutations
plot_simData(sim, data = "omic2", permute = TRUE, permute_seed = 123)

# Inspect the simulated sample-level factor scores
plot_factor(sim, factor_num = 1, type = "scatter")
```

Matrices in `sim$omics` contain **samples in rows and features in columns**. Check the orientation required by your downstream method before using them.

## Configuring the simulation

| Argument | Purpose |
| --- | --- |
| `vector_features` | Feature count for each omics layer; provide at least two layers. |
| `n_samples` | Number of shared samples. |
| `n_factors` | Number of latent factors to generate. |
| `snr` | Signal-to-noise variance ratio used to scale signal or noise. |
| `signal.samples` | Mean and standard deviation of sample scores within signal regions. |
| `signal.features` | One mean/standard deviation vector per layer, used to generate signal loadings. |
| `factor_structure` | Distribution of factors across layers, as described below. |
| `num.factor` | `"multiple"` or `"single"`; single-factor mode uses one factor. |
| `seed` | Random seed for reproducible simulation. |
| `real_stats` | Use supplied background means and variances when `TRUE`. |
| `real_means_vars` | One named `c(mean = ..., var = ...)` vector per layer when `real_stats = TRUE`. |

Supply `signal.features` explicitly as a list with the same length as `vector_features`, even though its formal default is `NULL`.

### Factor structures

| Setting | Factor allocation |
| --- | --- |
| `"shared"` | Each factor affects every layer. |
| `"unique"` | Each factor affects one randomly selected layer. |
| `"mixed"` | Each factor affects a randomly selected, nonempty subset of layers. |
| `"partial"` | Each factor affects at least two layers and fewer than all layers when there are more than two. With two layers, both are affected. |

Single-factor mode supports `"shared"`, `"unique"`, and `"partial"`. The `"custom"` option is not yet implemented.

The simulator produces continuous data through latent factor contributions and Gaussian background noise. Known signal regions and factor allocations provide reference information for evaluating recovery; they should not be treated automatically as clinical labels.

## Simulation output

| Component | Contents |
| --- | --- |
| `omics` | Named list of matrices: `omic1`, `omic2`, and so on. |
| `concatenated_datasets[[1]]` | Matrix formed by joining the feature columns of all layers. |
| `list_alphas` | Sample-level scores for each latent factor. |
| `list_betas` | Feature loadings by layer and factor. |
| `signal_annotation$samples` | Sample indices assigned to each factor's signal region. |
| `signal_annotation$features` | Signal feature indices by layer and factor. |
| `factor_map` | Layers affected by each factor. |
| `factor_structure` | Requested factor structure. |

## Other functions

| Function | Purpose |
| --- | --- |
| `simulate_twoOmicsData()` | Generate two-layer simulations using the alternative two-omics interface. |
| `as_multiomics()` | Convert legacy two-omics outputs to the standard multi-omics structure. |
| `plot_simData()` | Plot merged or individual-layer heatmaps, with optional permutations. |
| `plot_factor()` | Plot sample-level factor scores as scatter plots or histograms. |
| `plot_weights()` | Plot feature loadings as scatter plots or histograms. |
| `compute_means_vars()` | Summarize overall, row-wise, and column-wise means and standard deviations. |
| `demo_multiomics_analysis()` | Run a MOFA2 demonstration on simulated or CLL data, with optional PowerPoint export. |
| `sumo_setup_mofa()` | Configure a reticulate Python environment for MOFA training. |

For legacy two-omics simulations, standardize the output before using functions that expect `omics`:

```r
legacy <- simulate_twoOmicsData(
  vector_features = c(4000, 3000),
  n_samples = 100,
  n_factors = 2,
  snr = 2.5,
  num.factor = "multiple",
  advanced_dist = "mixed"
)

standard <- as_multiomics(legacy)
plot_simData(standard, data = "merged")
```

## Optional MOFA2 workflow

The analysis demonstration requires additional packages. Install the Bioconductor dependencies separately:

```r
install.packages("BiocManager")
BiocManager::install(c("MOFA2", "MOFAdata", "basilisk"))

# Train a model; requires a working Python backend for MOFA2
demo_multiomics_analysis(
  data_type = "SUMO",
  export_pptx = FALSE,
  use_pretrained = "never"
)
```

Training prefers MOFA2's basilisk backend. An interactive reticulate setup is also available through `sumo_setup_mofa()`. Set `export_pptx = TRUE` to write a PowerPoint report.

Pretrained model loading is supported only when the corresponding artifacts are available. Check availability with `sumo_pretrained_mofa_available()` before requesting pretrained mode.

## Documentation

- [HTML walkthrough](SUMO_vignette_v2.0.html): download and open in a browser.
- [Function reference files](man/).
- [Changes in version 1.2.2](NEWS.md).

In R, open function help with:

```r
?SUMO
?simulateMultiOmics
?simulate_twoOmicsData
?demo_multiomics_analysis
```

## Authors and contact

**Author and maintainer:** Bernard Isekah Osang'ir

**Contributors:** Ziv Shkedy, Surya Gupta, and Jürgen Claesen

- Email: [Bernard.Osangir@sckcen.be](mailto:Bernard.Osangir@sckcen.be)
- ORCID: [0000-0002-5557-3602](https://orcid.org/0000-0002-5557-3602)
- Bug reports and feature requests: [GitHub Issues](https://github.com/lucp12891/SUMO/issues)

## License

SUMO is licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE.md).
