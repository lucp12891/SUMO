## SUMO 1.2.1

### Minor Enhancements
- **Noise modeling using real data statistics**:
  Added `real_stats = TRUE` and `real_means_vars` arguments to `simulateMultiOmics()`.
  This allows users to simulate omics datasets where the background noise is modeled
  to match real-world omic-specific mean and variance profiles.
- **Signal-free omics still receive noise**:
  Omics datasets that are not mapped to any latent factor now still receive
  noise values, ensuring more realistic simulation structure across all omics.
- **Improved error handling**:
  Enhanced error messages for better clarity when issues arise during simulation.
  This includes more informative messages when omics are not mapped to latent factors.
- **Documentation updates**:
  Updated documentation to reflect the new features and improvements in the package.
  This includes examples and explanations for the new `real_stats` functionality.
- **Code optimization**:
  Optimized internal code for better performance and efficiency in simulations.
  This includes improvements in how noise is generated and applied to omics datasets.
