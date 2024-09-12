test_that("dividing features for omic two", {
  result<-feature_selection_two(n_features_two = 2500, num.factor = "multiple", no_factor = 4)
  expect_type(result, "list")
  expect_true(length(result) > 0)
})
