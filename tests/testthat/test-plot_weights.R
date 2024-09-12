# Load necessary libraries
library(testthat)
library(SUMO)  # Load your package
output_sim <- OmixCraftHD(vector_features = c(2000,3000), sigmas_vector=c(4,5), n_samples=100, n_factors=5, num.factor='multiple', advanced_dist='mixed')

test_that("plot_weights runs without error", {
  # Test that the plot function runs without errors
  expect_no_error(plot_weights(output_sim, factor_num = 1, data = 'omic.one', type = 'scatter'))
})
#expect_silent(plot_weights(output_sim, factor_num = 1, data = 'omic.one', type = 'scatter'))
