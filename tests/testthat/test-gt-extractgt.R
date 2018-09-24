context("Extract Genotype")

vcf  <- readRDS('res-extractgt/vcf.rds')
surv <- readRDS('res-extractgt/survivals.rds')
gt   <- ExtractGt(vcf)

test_that("ExtractGt returns a data frame", {
  expect_is(gt, "data.frame")
})

test_that("returned data frame contains strings", {
  expect_is(gt[1,1], "character")
  expect_is(gt[dim(gt)[1], dim(gt)[2]], "character")
})

test_that("dimensions of returned data frame are correct", {
  expect_equal(dim(gt), c(3109, 109))
  expect_equal(sum(is.na(gt)), 23028)
})

test_that("the min.depth setting works", {
  censored.gt <- ExtractGt(vcf, min.depth = 5)
  expect_equal(dim(censored.gt), c(3109, 109))
  expect_equal(sum(is.na(censored.gt)), 162520)
})

test_that("the min.sample.qual setting works", {
  omitted.gt <- ExtractGt(vcf, min.sample.qual = 0.8)
  expect_equal(dim(omitted.gt), c(3109, 96))
  expect_equal(sum(is.na(omitted.gt)), 11467)
})

test_that("the min.variant.qual setting works", {
  omitted.gt <- ExtractGt(vcf, min.variant.qual = 0.9)
  expect_equal(dim(omitted.gt), c(2264, 109))
  expect_equal(sum(is.na(omitted.gt)), 10536)
})

test_that("the include.depth option include the depth matrix", {
  depth.included.gt <- ExtractGt(vcf, include.depth = TRUE)
  expect_is(depth.included.gt, "list")
  expect_equal(length(depth.included.gt), 2)
  expect_equal(dim(depth.included.gt$gt), dim(depth.included.gt$gt))
})

test_that("the depth matrix is affected by omitting / censoring settings", {
  depth.included.gt <- ExtractGt(vcf, include.depth = TRUE,
                                 min.variant.qual = 0.9, min.sample.qual = 0.8,
                                 min.depth = 5)
  expect_equal(dim(depth.included.gt$gt), dim(depth.included.gt$gt))
})

test_that("warning is thrown when min.depth is not a positive integer", {
  expect_warning(ExtractGt(vcf, min.depth = -1), "min.depth")
  expect_warning(ExtractGt(vcf, min.depth = 0.42), "min.depth")
})

test_that("error is thrown when min.depth is not a numeric or NULL", {
  expect_error(ExtractGt(vcf, min.depth = "tortilla"), "min.depth")
  expect_error(ExtractGt(vcf, min.depth = NA), "min.depth")
})

test_that("warning is thrown when min.variant.qual is not between 0 and 1", {
  expect_warning(ExtractGt(vcf, min.variant.qual = -1), "min.variant.qual")
  expect_warning(ExtractGt(vcf, min.variant.qual = 42), "min.variant.qual")
})

test_that("error is thrown when min.variant.qual is not a numeric or NULL", {
  expect_error(ExtractGt(vcf, min.variant.qual = "burritos"), "min.variant.qual")
  expect_error(ExtractGt(vcf, min.variant.qual = NA), "min.variant.qual")
})

test_that("warning is thrown when min.sample.qual is not between 0 and 1", {
    expect_warning(ExtractGt(vcf, min.variant.qual = -1), "min.sample.qual")
  expect_warning(ExtractGt(vcf, min.variant.qual = 777), "min.sample.qual")
})

test_that("error is thrown when min.sample.qual is not a numeric or NULL", {
  expect_error(ExtractGt(vcf, min.variant.qual = "nachos"), "min.sample.qual")
  expect_error(ExtractGt(vcf, min.variant.qual = NA), "min.sample.qual")
})