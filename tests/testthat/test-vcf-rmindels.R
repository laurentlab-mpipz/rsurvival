context("Remove Indels")

vcf <- readRDS('res-rmindels\\vcf.rds')
vcf.list <- list("alive" = vcf, "dead" = vcf)
cleared.vcf <- RmIndels(vcf, verbose = FALSE)
cleared.vcf.list <- RmIndels(vcf.list, verbose = FALSE)

test_that("RmIndels returns a vcfR object", {
  expect_is(cleared.vcf, "vcfR")
})

test_that("error is thrown when parameter vcf is not a vcfR Object", {
  expect_error(RmIndels(NA), "vcf")
  expect_error(RmIndels("pineapple"), "vcf")
  expect_error(RmIndels(c(1,2,3)), "vcf")
})

test_that("the returned vcfR object has the correct length", {
  expect_true(nrow(cleared.vcf) == 135)
})

test_that("the returned vcfR object contains no indels", {
  expect_equal(vcfR::extract.indels(vcf), cleared.vcf)
})

test_that("lists of vcfR objects are correctly handled", {
  expect_is(cleared.vcf.list$alive, "vcfR")
  expect_is(cleared.vcf.list$dead, "vcfR")
  expect_true(nrow(cleared.vcf.list$alive) == 135)
  expect_true(nrow(cleared.vcf.list$dead) == 135)
  expect_equal(vcfR::extract.indels(vcf.list$alive), cleared.vcf.list$alive)
  expect_equal(vcfR::extract.indels(vcf.list$dead), cleared.vcf.list$dead)
})