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

  expect_equal(
    MicroarrayPreprocessConfig(
      "GSE12345", mock_db, keep_sample_fn = "identity"
    )@keep_sample_fn,
    expected = identity,
    info = paste(
      "`MicroarrayPreprocessConfig()` can take a `keep_sample_fn` as a string"
    )
  )

  expect_equal(
    MicroarrayPreprocessConfig(
      "GSE12345", mock_db, keep_probe_fn = "identity"
    )@keep_probe_fn,
    expected = identity,
    info = paste(
      "`MicroarrayPreprocessConfig()` can take a `keep_probe_fn` as a string"
    )
  )

  expect_equal(
    MicroarrayPreprocessConfig(
      "GSE12345", "org.Hs.eg.db::org.Hs.eg.db"
    )@entrezgene_db,
    expected = org.Hs.eg.db::org.Hs.eg.db,
    info = paste(
      "`MicroarrayPreprocessConfig()` can take an `entrezgene_db` as a string"
    )
  )
})
