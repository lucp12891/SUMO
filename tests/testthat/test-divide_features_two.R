test_that("divide features for omic two", {
  result<-divide_features_two(n_features_two = 2000, num.factor = 'single')
  expect_type(result, "list")
  expect_true(length(result) > 0)
})
