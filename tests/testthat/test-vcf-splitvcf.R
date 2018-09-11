context("Split VCF")

vcf <- readRDS('res-splitvcf\\vcf.rds')
survivals <- readRDS('res-splitvcf\\survivals.rds')
split <- SplitVcf(vcf, survivals)
nbcol.vcf <- dim(vcf)["gt_cols"]

test_that("SplitVcf returns vcfR objects in a list of vcfR objects", {
  expect_is(split, "list")
  expect_is(split$alive, "vcfR")
  expect_is(split$dead, "vcfR")
})

test_that("SplitVcf returns vcfR objects with a correct length", {
  expect_equal(nrow(split$alive), nrow(vcf))
  expect_equal(nrow(split$dead), nrow(vcf))
  expect_equal(nrow(split$alive), nrow(vcf))
  expect_equal(dim(split$alive)['gt_cols'] + dim(split$dead)['gt_cols'],
                nbcol.vcf + 1)
  expect_true(dim(split$alive)['gt_cols'] == 1 + sum(survivals))
  expect_true(dim(split$dead)['gt_cols'] == nbcol.vcf - sum(survivals))
})

test_that("error is thrown when parameter vcf is not a vcfR object", {
  expect_error(SplitVcf(42, survivals), "vcf")
  expect_error(SplitVcf("brownie", survivals), "vcf")
})

test_that("error is thrown when parameter survivals is not a logical vector", {
  expect_error(SplitVcf(vcf, 42), "survival")
  expect_error(SplitVcf(vcf, "brownie"), "survival")
})

test_that("warning is thrown if survivals vector is too short", {
  expect_warning(SplitVcf(vcf, c(TRUE)), "shorter")
})

test_that("error is thrown if survivals vector is too long", {
  expect_error(SplitVcf(vcf, c(survivals, TRUE)), "longer")
})