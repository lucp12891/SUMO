test_that("dividing features for omic one", {
  result<-feature_selection_one(n_features_one = 2000, num.factor = "multiple", no_factor = 3)
  expect_type(result, "list")
  expect_true(length(result) > 0)
})
