context("Split Genotype")

gt <- readRDS('res-splitgt/gt.rds')
surv <- readRDS('res-splitgt/survivals.rds')
splitted <- SplitGt(gt, surv)

test_that("SplitGt returns a list of data frame", {
  expect_is(gt, "list")
  expect_is(gt$alive, "data.frame")
  expect_is(gt$dead, "data.frame")
})

test_that("returned data frames contain strings", {
  expect_is(gt$alive[1,1], "character")
  expect_is(gt$dead[1,1], "character")
})

test_that("dimensions of returned data frames are correct", {
  expect_equal(dim(gt$alive), c(nrow(gt), sum(as.numeric(surv))))
  expect_equal(dim(gt$dead), c(nrow(gt), sum(as.numeric(!surv)))))
})

test_that("warning is thrown when survival is too short", {
  expect_warning(SplitGt(gt, surv[1:42]), "survival")
})

test_that("error is thrown when survival is too long", {
  expect_error(SplitGt(gt, c(surv, TRUE)), "survival")
})

test_that("error is thrown when survival is not a logical vector", {
  expect_error(SplitGt(gt, "ribs"), "survival")
  expect_error(SplitGt(gt, 42), "survival")
})

test_that("error is thrown when gt is not a matrix or data.frame", {
  expect_error(SplitGt("tenders", surv), "gt")
  expect_error(c(1:42), surv), "gt")
})