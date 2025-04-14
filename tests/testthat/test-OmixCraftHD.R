test_that("generate data", {
  #result <- OmixCraftHD(vector_features = c(2000,2000), sigmas_vector=c(4,5), n_samples=100, n_factors=3, num.factor='multiple', advanced_dist='mixed')
  result <- OmixCraftHD(
     vector_features = c(3000, 2500),
     n_samples = 100,
     n_factors = 5,
     signal.samples = c(3, 0.5),          # mean = 3 samples, variance = 0.5
     signal.features.one = c(5, 2.0),     # mean = 5  in omic1, variance = 2.0
     signal.features.two = c(4.5, 0.05),  # mean = 4.5 in omic2, variance = 0.05
     snr = 2,
     num.factor = 'multiple',
     advanced_dist = 'mixed'
   )
  expect_type(result, "list")
  expect_true(length(result) > 0)
  })
