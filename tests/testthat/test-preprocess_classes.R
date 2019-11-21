context("Tests for microarray-preprocessing classes")

test_that("MicroarrayPreprocessConfig", {
  expect_error(
    MicroarrayPreprocessConfig(), info = "No arguments to constructor"
  )

  mock_db <- list()
  class(mock_db) <- "OrgDb"

  expect_is(
    MicroarrayPreprocessConfig(
      "GSE12345", mock_db
    ),
    "MicroarrayPreprocessConfig",
    info = "MicroarrayPreprocessConfig() constructor returns the correct class"
  )

  # TODO:
  # user can pass the name of an OrgDb instead of an actual OrgDb, and the
  # constructor will get pass the correct database to MicroarrayPreprocessConfig
})
