# Load necessary libraries
library(testthat)
library(SUMO)  # Load your package
#output_sim <- simulate_twoOmicsData(vector_features = c(2000,3000), sigmas_vector=c(4,5), n_samples=100, n_factors=5, num.factor='multiple', advanced_dist='mixed')
output_obj <- simulate_twoOmicsData(
   vector_features = c(4000,3000),
   n_samples = 100,
   n_factors = 2,
   snr = 2.5,
   num.factor = 'multiple',
   advanced_dist = 'mixed')
test_that("plot_factor runs without error", {
  # Test that the plot function runs without errors
  expect_silent(plot_factor(output_obj, factor_num = 1))
})
#expect_silent(plot_factor(output_obj, factor_num = 1))
