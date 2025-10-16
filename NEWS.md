## SUMO 1.2.2

### New

* **Pretrained artifacts for fast demos/tests**
  Added cached, pre-trained demo objects (simulated outputs and factor structures) to speed up examples, vignettes, and CI cold starts.
* **Optional backend support: basilisk & reticulate**
  Introduced backend hooks for isolated Python environments via **basilisk** and interoperability via **reticulate** (aka “reticulare” in earlier notes). Includes auto-detection and clearer setup guidance.
* **Demo mode in `simulate_twoOmicsData()`**
  `multiomics_demo = "SHOW"` provides a ready-to-run two-factor toy example for quick multi-omics walkthroughs.

### Enhancements

* **Legacy → standard conversion pipeline**
  Polished `as_multiomics()` to reliably upgrade legacy `simulate_twoOmicsData()` outputs into the current schema (`omics`, per-omic `list_betas`, `signal_annotation$features`, `factor_map`, etc.), including:

  * normalization of factor-score names (`alpha{n}`),
  * mapping `delta{n}` → `beta{n}` for omic 2,
  * inference of sample/feature annotations when explicit indices are missing,
  * robust split of `concatenated_datasets` back into `omic1`/`omic2` using counts or `omic1_`/`omic2_` prefixes.
* **Internal helpers made non-exported (with docs scaffolding)**
  `%||%` for defaults and `.is_multiomics_like()` for schema checks are now internal, with roxygen headers (no code changes).
* **Stronger standardization & resilience**
  Automatic row/column naming for omics matrices; safer handling of empty/`NA` index vectors; improved messages for sample-partition constraints; sturdier concatenated-split heuristics.

### Documentation

* **Vignette updates**
  Expanded “Getting Started” and migration guidance: using demo mode, converting legacy objects with `as_multiomics()`, and configuring basilisk/reticulate backends.
* **Function docs**
  Roxygen scaffolding added for the internal utilities and converter while preserving the original function bodies.

### Fixes

* Consistent `delta`→`beta` mapping across omics.
* Deterministic naming and factor-label normalization.
* More informative errors when sample segments are too small or factors are unmapped.

### Breaking Changes

* **None.** Legacy outputs remain supported; call `as_multiomics()` to standardize them for downstream workflows.
