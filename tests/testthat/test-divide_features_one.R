test_that("divide features for omic one", {
  result<-divide_features_one(n_features_one = 2000, num.factor = 'single')
  expect_type(result, "list")
  expect_true(length(result) > 0)
})

