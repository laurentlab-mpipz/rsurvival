context("Slice data frame rows")

df <- data.frame("banana" = c(1:3), "sugar"= c(4:6), "apple" = c(7:9))

slice.one <- SliceDfRows(df, 2)
slice.multi <- SliceDfRows(df, c(1, 3))

test_that("SliceDfRows returns a data frame", {
  expect_is(slice.one, "data.frame")
  expect_is(slice.multi, "data.frame")
})

test_that("returned data frame has the correct dimensions", {
  expect_equal(dim(slice.one), c(1, 3))
  expect_equal(dim(slice.multi), c(2, 3))
})

test_that("returned data frame kept the row names", {
  expect_equal(rownames(slice.one), "2")
  expect_equal(rownames(slice.multi), c("1", "3"))
})

test_that("returned data frame kept the column names", {
  expect_equal(colnames(slice.one), c("banana", "sugar", "apple"))
  expect_equal(colnames(slice.multi), c("banana", "sugar", "apple"))
})
