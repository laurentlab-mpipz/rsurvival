context("Is At Least NA ratio")

test_that("NAs are correctly counted", {
  expect_equal(IsAtLeastNARatio(c(1, 2, 3, NA), 0.24), TRUE)
  expect_equal(IsAtLeastNARatio(c(1, 2, 3, NA), 0.26), FALSE)
  expect_equal(IsAtLeastNARatio(c(NA), 1), TRUE)
  expect_equal(IsAtLeastNARatio(c(1), 0), TRUE)
})

test_that("NA is returned when vect is empty", {
  expect_equal(IsAtLeastNARatio(c(), 0.5), NA)
})

test_that("warning is thrown when min.ratio is not relevant", {
  expect_warning(IsAtLeastNARatio(c(1,2,3), -1), "min.ratio")
  expect_warning(IsAtLeastNARatio(c(1,2,3), 42), "min.ratio")
})

test_that("error is thrown when min.ratio is not a number", {
  expect_error(IsAtLeastNARatio(c(1,2,3), "a"), "min.ratio")
  expect_error(IsAtLeastNARatio(c(1,2,3), c(1,2,3)), "min.ratio")
  expect_error(IsAtLeastNARatio(c(1,2,3), TRUE), "min.ratio")
  expect_error(IsAtLeastNARatio(c(1,2,3), NULL), "min.ratio")
  expect_error(IsAtLeastNARatio(c(1,2,3), NA), "min.ratio")
})