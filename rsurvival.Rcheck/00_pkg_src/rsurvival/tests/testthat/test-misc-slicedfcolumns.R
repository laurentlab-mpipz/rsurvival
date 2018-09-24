context("Slice data frame columns")

df <- data.frame("banana" = c(1:3), "sugar"= c(4:6), "apple" = c(7:9))

slice.one <- SliceDfColumns(df, 2)
slice.multi <- SliceDfColumns(df, c(1, 3))

test_that("SliceDfColumns returns a data frame", {
  expect_is(slice.one, "data.frame")
  expect_is(slice.multi, "data.frame")
})

test_that("returned data frame has the correct dimensions", {
  expect_equal(dim(slice.one), c(3, 1))
  expect_equal(dim(slice.multi), c(3, 2))
})

test_that("returned data frame kept the row names", {
  expect_equal(rownames(slice.one), c("1", "2", "3"))
  expect_equal(rownames(slice.multi), c("1", "2", "3"))
})

test_that("returned data frame kept the column names", {
  expect_equal(colnames(slice.one), "sugar")
  expect_equal(colnames(slice.multi), c("banana", "apple"))
})
