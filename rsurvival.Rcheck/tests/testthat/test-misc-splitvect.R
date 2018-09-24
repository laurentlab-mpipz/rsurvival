context("Split Vector")

vect <- c(0,1,2,3,4,5,6,7,8,9)
surv <- c(1,0,1,0,1,1,0,1,1,1)
surv <- as.logical(surv)

split <- SplitVect(vect, surv)

test_that("SplitVect returns a list of vector", {
  expect_is(split, "list")
  expect_is(split$alive, "numeric")
  expect_is(split$dead, "numeric")
})

test_that("returned vector have the right length", {
  expect_equal(length(split$alive), 7)
  expect_equal(length(split$dead), 3)
})

test_that("warning is thrown when parameter survival is too short", {
  expect_warning(SplitVect(vect, c(TRUE, FALSE)), "survival")
  expect_warning(SplitVect(vect, NA), "survival")
})

test_that("error is thrown when parameter survival is too long", {
  expect_error(SplitVect(vect, as.logical(c(0:10))),  "survival")
})