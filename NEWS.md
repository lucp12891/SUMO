## SUMO 1.2.0

### 🧪 Major Simulation Engine Upgrade

* ✅ **Support for simulating >2 omics datasets:** A core feature upgrade enabling generation of multi-omics datasets with **three or more layers**, such as transcriptomics, proteomics, and metabolomics, in a single simulation run.
* 🔄 **Function renamed for clarity:**

  * Previous: `OmixCraftHD` → Now: `simulate_twoOmicsData` (for 2-omics simulation).
  * New: `OmixCraftMultiHD` (alias `simulateMultOmics`) introduced to simulate an **arbitrary number of omics** using the same generative framework.
* 💡 Non-overlapping latent factors with flexible per-omic and per-factor signal regions are now seamlessly handled for multiple datasets.

### 🧬 Factor Model Design Enhancements

* Fully **modular generative factor model** allows:

  * Signal assignment via latent factors per omic.
  * `num.factor = "single"` or `"multiple"` control.
  * `factor_structure` for each factor: `shared`, `unique`, `mixed`, or `partial`.
* Sample blocks and feature blocks are now simulated with **sequential, non-overlapping indices**, respecting biological plausibility and avoiding signal bleed.

### 🛠️ Refinements and Validation

* Robust input checks ensure valid parameterization for complex simulation designs.
* All alpha (sample score) vectors and beta (feature weight) vectors now support sparse signal encoding and customized signal-to-noise ratios (`snr`).
* Each omics layer can be assigned its own mean and standard deviation for signal blocks, increasing realism.

### 📈 Documentation and Usability

* Improved documentation and usage examples for:

  * `simulateMultiOmics()`
  * Factor assignment logic and feature-sample block allocation
* Example usage includes:

  ```r
  sim_object <- simulateMultiOmics(
    vector_features = c(3000, 2500, 2000),
    n_samples = 100,
    n_factors = 3,
    snr = 3,
    signal.samples = c(5, 1),
    signal.features = list(c(3, 0.3), c(2.5, 0.25), c(2, 0.2)),
    factor_structure = "mixed",
    num.factor = "multiple",
    seed = 123
  )
  ```

### 📦 Infrastructure

* Full compatibility with `Roxygen2`, `pkgdown`, and `devtools::check()` (0 errors, 0 warnings, 0 notes).
* Added CRAN-safe global variable declarations.

