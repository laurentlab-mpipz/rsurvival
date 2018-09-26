context("Remove Non Biallelics")

vcf <- readRDS('res-rmnonbiallelics/vcf.rds')
vcf.list <- list("alive" = vcf, "dead" = vcf)
cleared.vcf <- RmNonBiallelics(vcf, verbose = FALSE)
cleared.vcf.list <- RmNonBiallelics(vcf.list, verbose = FALSE)

test_that("RmNonBiallelics returns a vcfR object", {
  expect_is(cleared.vcf, "vcfR")
})

test_that("error is thrown when parameter vcf is not a vcfR Object", {
  expect_error(RmNonBiallelics(NA), "vcf")
  expect_error(RmNonBiallelics("pineapple"), "vcf")
  expect_error(RmNonBiallelics(c(1,2,3)), "vcf")
})

test_that("the returned vcfR object has the correct length", {
  expect_true(nrow(cleared.vcf) == 141)
})

test_that("the returned vcfR object contains only biallelics", {
  expect_equal(sum(!vcfR::is.biallelic(cleared.vcf)), 0)
})

test_that("lists of vcfR objects are correctly handled", {
  expect_is(cleared.vcf.list$alive, "vcfR")
  expect_is(cleared.vcf.list$dead, "vcfR")
  expect_true(nrow(cleared.vcf.list$alive) == 141)
  expect_true(nrow(cleared.vcf.list$dead) == 141)
  expect_equal(sum(!vcfR::is.biallelic(cleared.vcf.list$alive)), 0)
  expect_equal(sum(!vcfR::is.biallelic(cleared.vcf.list$dead)), 0)
})