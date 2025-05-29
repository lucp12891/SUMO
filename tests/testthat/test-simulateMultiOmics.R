test_that("generate data", {
  result <- simulateMultiOmics(
               vector_features = c(3000, 2500, 2000),
               n_samples = 100,
               n_factors = 3,
               snr = 3,
               signal.samples = c(5, 1),
               signal.features = list(
                 c(3, 0.3),   # omic1 signal mean/sd
                 c(2.5, 0.25),# omic2 signal mean/sd
                 c(2, 0.2)    # omic3 signal mean/sd
               ),
               factor_structure = "mixed",
               num.factor = "multiple",
               seed = 123
               )
  expect_type(result, "list")
  expect_true(length(result) > 0)
  })
