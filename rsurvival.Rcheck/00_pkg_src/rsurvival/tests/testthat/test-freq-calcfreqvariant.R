context("Calculate Frequency of Variant")

tol <- 1e-4

gt <- readRDS('res-calcfreqvariant/gt.rds')

freq.standard <- CalcFreqVariant(gt[1, ])
freq.nothing  <- CalcFreqVariant(gt[1, ], genotypic = FALSE, allelic = FALSE)
freq.geno     <- CalcFreqVariant(gt[1, ], genotypic = TRUE, allelic = FALSE)
freq.allel    <- CalcFreqVariant(gt[1, ], genotypic = FALSE, allelic = TRUE)
freq.geno.al  <- CalcFreqVariant(gt[1, ], genotypic = TRUE, allelic = TRUE)
freq.nototals <- CalcFreqVariant(gt[1, ], genotypic = TRUE, allelic = TRUE,
                                 totals = FALSE)

freq.allel.rel  <- CalcFreqVariant(gt[1, ], genotypic = FALSE, allelic = TRUE,
                                    absolute = FALSE)
freq.allel.perc <- CalcFreqVariant(gt[1, ], genotypic = FALSE, allelic = TRUE,
                                    absolute = FALSE, percentage = TRUE)
freq.al.noextra <- CalcFreqVariant(gt[1, ], genotypic = FALSE, allelic = TRUE,
                                    absolute = FALSE, extrapolate = FALSE)

CountMatchName <- function(freq, regex){
  result <- sum(as.numeric(grepl(regex, names(freq))))   
}

test_that("CalcFreqVariant returns a vector", {
  expect_null(dim(freq.standard))
})

test_that("returned vector is numeric", {
  expect_true(is.numeric(freq.standard))
})

test_that("length of returned vector is correct", {

  expect_equal(length(freq.nothing), 0)
  expect_equal(length(freq.geno), 5)
  expect_equal(length(freq.allel), 4)
  expect_equal(length(freq.geno.al), 9)
  expect_equal(length(freq.nototals), 7)

})

test_that("column names are correct", {

  expect_equal(CountMatchName(freq.geno, "gt"), 5)
  expect_equal(CountMatchName(freq.geno, "count"), 5)
  expect_equal(CountMatchName(freq.geno, "al"), 0)
  expect_equal(CountMatchName(freq.geno, "freq"), 0)
  expect_equal(CountMatchName(freq.geno, "perc"), 0)

  expect_equal(CountMatchName(freq.allel.rel, "gt"), 0)
  expect_equal(CountMatchName(freq.allel.rel, "count"), 1)
  expect_equal(CountMatchName(freq.allel.rel, "al"), 4)
  expect_equal(CountMatchName(freq.allel.rel, "freq"), 3)
  expect_equal(CountMatchName(freq.allel.rel, "perc"), 0)

  expect_equal(CountMatchName(freq.allel.perc, "gt"), 0)
  expect_equal(CountMatchName(freq.allel.perc, "count"), 1)
  expect_equal(CountMatchName(freq.allel.perc, "al"), 4)
  expect_equal(CountMatchName(freq.allel.perc, "freq"), 0)
  expect_equal(CountMatchName(freq.allel.rel, "perc"), 0)

})

test_that("count values are correct", {

  expect_equal(as.numeric(freq.geno["count.gt.HOMOREF"]), 9)
  expect_equal(as.numeric(freq.geno["count.gt.HETERO"]), 28)
  expect_equal(as.numeric(freq.geno["count.gt.HOMOALT"]), 66)
  expect_equal(as.numeric(freq.geno["count.gt.MISSVAL"]), 6)
  expect_equal(as.numeric(freq.geno["count.gt.TOTAL"]), 109)

  expect_equal(as.numeric(freq.allel["count.al.REF"]), 46)
  expect_equal(as.numeric(freq.allel["count.al.ALT"]), 160)
  expect_equal(as.numeric(freq.allel["count.al.MISSVAL"]), 12)
  expect_equal(as.numeric(freq.allel["count.al.TOTAL"]), 218)

})

test_that("frequency values are correct", {

  expect_equal(as.numeric(freq.allel.rel["freq.al.REF"]), 0.2233,
               tolerance = tol)
  expect_equal(as.numeric(freq.allel.rel["freq.al.ALT"]), 0.7767,
               tolerance = tol)
  expect_equal(as.numeric(freq.allel.rel["freq.al.MISSVAL"]), 0.055,
               tolerance = tol)

})

test_that("percentage values are correct", {

  expect_equal(as.numeric(freq.allel.perc["perc.al.REF"]), 22.33,
               tolerance = tol)
  expect_equal(as.numeric(freq.allel.perc["perc.al.ALT"]), 77.67,
               tolerance = tol)
  expect_equal(as.numeric(freq.allel.perc["perc.al.MISSVAL"]), 5.505,
               tolerance = tol)

})

test_that("frequency values without extrapolation are correct", {

  expect_equal(as.numeric(freq.al.noextra["freq.al.REF"]), 0.2111,
               tolerance = tol)
  expect_equal(as.numeric(freq.al.noextra["freq.al.ALT"]), 0.734,
               tolerance = tol)
  expect_equal(as.numeric(freq.al.noextra["freq.al.MISSVAL"]), 0.055,
               tolerance = tol)

})

test_that("min.freq.gt option works", {

  full.of.nas <- CalcFreqVariant(gt[1, ], min.freq.gt = 0.09)
  full.of.values <- CalcFreqVariant(gt[1, ], min.freq.gt = 0.081)

  expect_equal(sum(is.na(full.of.nas)), length(full.of.nas))
  expect_equal(sum(is.na(full.of.values)), 0) 

})

test_that("min.freq.al option works", {

  full.of.nas <- CalcFreqVariant(gt[1, ], min.freq.al = 0.25)
  full.of.values <- CalcFreqVariant(gt[1, ], min.freq.al = 0.21)

  expect_equal(sum(is.na(full.of.nas)), length(full.of.nas))
  expect_equal(sum(is.na(full.of.values)), 0) 

})

test_that("warning is thrown when min.freq.al is out of [0:1]", {

  expect_warning(CalcFreqVariant(gt[1,], min.freq.al = 3), "min.freq.al")
  expect_warning(CalcFreqVariant(gt[1,], min.freq.al = -1), "min.freq.al")

})

test_that("warning is thrown when min.freq.gt is out of [0:1]", {

  expect_warning(CalcFreqVariant(gt[1,], min.freq.gt = 3), "min.freq.gt")
  expect_warning(CalcFreqVariant(gt[1,], min.freq.gt = 3), "min.freq.gt")

})

test_that("error is thrown when min.freq.al is not a number", {

  expect_error(CalcFreqVariant(gt[1,], min.freq.al = "pear"), "min.freq.al")
  expect_error(CalcFreqVariant(gt[1,], min.freq.al = NA), "min.freq.al")

})

test_that("error is thrown when min.freq.gt is out of [0:1]", {

  expect_error(CalcFreqVariant(gt[1,], min.freq.gt = "apple"), "min.freq.gt")
  expect_error(CalcFreqVariant(gt[1,], min.freq.gt = NA), "min.freq.gt")

})