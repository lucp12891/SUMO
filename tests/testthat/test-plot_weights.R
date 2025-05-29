library(testthat)
library(SUMO)  # Load your package
test_that("plot_weights runs without error", {
  output_obj <- simulate_twoOmicsData(
     vector_features = c(4000, 3000),
     n_samples = 100,
     n_factors = 2,
     signal.samples = NULL,
     signal.features.one = NULL,
     signal.features.two = NULL,
     snr = 2.5,
     num.factor = 'multiple',
     advanced_dist = 'mixed'
     )

  # Extract a valid factor number dynamically
  available_factors <- as.numeric(gsub("[^0-9]", "", names(output_obj$list_betas)))
  expect_true(length(available_factors) > 0)

  test_factor <- available_factors[1]

  expect_silent(
    plot_weights(
      sim_object = output_obj,
      factor_num = test_factor,
      data = "omic.one",
      type = "scatter",
      show.legend = FALSE
    )
  )
})
