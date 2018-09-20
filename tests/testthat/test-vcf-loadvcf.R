context("Load VCF File")

test_that("LoadVcf returns a vcfR object", {
  vcf <- LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = FALSE,
                  only.snp = FALSE, verbose = FALSE)
  expect_equal(typeof(vcf), "S4")
  expect_is(vcf, "vcfR")
  expect_true(nrow(vcf) == 144)
})

test_that("LoadVcf filters non-biallelic variants", {
  vcf <- LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = TRUE,
                  only.snp = FALSE, verbose = FALSE)  
  expect_equal(typeof(vcf), "S4")
  expect_is(vcf, "vcfR")
  expect_true(nrow(vcf) == 141)
})

test_that("LoadVcf filters indel variants", {
  vcf <- LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = FALSE,
                  only.snp = , verbose = FALSE)  
  expect_equal(typeof(vcf), "S4")
  expect_is(vcf, "vcfR")
  expect_true(nrow(vcf) == 135)
})

test_that("LoadVcf can apply both filters", {
  vcf <- LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = TRUE,
                  only.snp = TRUE, verbose = FALSE)  
  expect_equal(typeof(vcf), "S4")
  expect_is(vcf, "vcfR")
  expect_true(nrow(vcf) == 132)
})

test_that("LoadVcf throws errors when logical inputs are not valid", {
  expect_error(LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = "sweets",
                only.snp = FALSE, verbose = FALSE),
                "only.biallelic")
  expect_error(LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = 42,
            only.snp = FALSE, verbose = FALSE),
            "only.biallelic")
  expect_error(LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = FALSE,
                only.snp = "sweets", verbose = FALSE),
                "only.snp")
  expect_error(LoadVcf('res-loadvcf/lollipop.vcf', only.biallelic = FALSE,
                only.snp = 42, verbose = FALSE),
                "only.snp")
})

test_that("LoadVcf throws errors when file.path input is not valid", {
  expect_error(LoadVcf('res-loadvcf\\donotexist.vcf', verbose = FALSE),
                "file.path")
  expect_error(LoadVcf('res-loadvcf\\lollipop.cinnamon', verbose = FALSE),
              "file.path")
  expect_error(LoadVcf(NA),
            "file.path")
  expect_error(LoadVcf(42),
            "file.path")
})