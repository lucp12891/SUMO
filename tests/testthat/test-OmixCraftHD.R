test_that("generate data", {
  result <- OmixCraftHD(vector_features = c(2000,2000), sigmas_vector=c(4,5), n_samples=100, n_factors=3, num.factor='multiple', advanced_dist='mixed')
  expect_type(result, "list")
  expect_true(length(result) > 0)
  })
