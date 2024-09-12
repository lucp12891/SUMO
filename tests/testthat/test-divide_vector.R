test_that("divide vectors", {
  result <- divide_samples(100,3,2)
  expect_type(result, "list")
  expect_true(length(result) > 0)
})
