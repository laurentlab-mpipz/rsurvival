context("Convert File To Logicals")

test_that("files made of 1s and 0s are supported", {
  expect_equal(ConvertFileToLogic('res-convertfiletologic/numbers',
                                  type = 'number'), 
                c(TRUE, FALSE, TRUE))
})

test_that("files made of 1s and 0s are supported", {
  expect_equal(ConvertFileToLogic('res-convertfiletologic/characters',
                                  type = 'character'),
                c(TRUE, TRUE, FALSE))
})

test_that("warning is thrown when result has no length", {
  expect_warning(ConvertFileToLogic('res-convertfiletologic/empty',
                                    type = 'number'),
                "empty", ignore.case = TRUE)
})

test_that("warning is thrown when result contains NAs", {
  expect_warning(ConvertFileToLogic('res-convertfiletologic/wrong',
                                    type = 'character'),
                "NA", ignore.case = TRUE)
})

test_that("error is thrown when type is unknown", {
  expect_error(ConvertFileToLogic('res-convertfiletologic/empty',
                                    type = 'cupcake'),
                "type", ignore.case = TRUE)
  expect_error(ConvertFileToLogic('res-convertfiletologic/empty',
                                    type = ''),
                "type", ignore.case = TRUE)
})

test_that("error is thrown when file does not exist", {
  expect_error(ConvertFileToLogic('res-convertfiletologic/cookies',
                                    type = 'number'),
                "file.path", ignore.case = TRUE)
})