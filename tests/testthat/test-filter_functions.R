###############################################################################

context("Tests for ExpressionSet filtering functions")

###############################################################################

# --- Custom test functions

check_fails_with_non_eset <- function(f_name = "") {
  f <- get(f_name)

  expect_error(
    f("Not an ExpressionSet"),
    info = paste(f_name, ": input should be an ExpressionSet, not a string")
  )

  expect_error(
    f(matrix(1:12, nrow = 4)),
    info = paste(f_name, ": input should be an ExpressionSet, not a matrix")
  )
}

check_fails_with_empty_dimensions <- function(f_name = "") {
  f <- get(f_name)

  no_probes <- random_eset(n_probes = 0)
  expect_error(
    object = f(no_probes),
    info = paste(f_name, ": input should have some probes")
  )

  no_samples <- random_eset(n_samples = 0)
  expect_error(
    object = f(no_samples),
    info = paste(f_name, ": input should have some samples")
  )

  empty <- random_eset(n_probes = 0, n_samples = 0)
  expect_error(
    object = f(empty),
    info = paste(f_name, ": input should not be empty")
  )
}

###############################################################################

# --- Tests

test_that("keep_all_probes", {

  # Input should be an ExpressionSet
  check_fails_with_non_eset("keep_all_probes")

  # Input should have some some probes / samples
  check_fails_with_empty_dimensions("keep_all_probes")

  # The returned values should just be [1 : num_probes]
  n_probes <- sample.int(50, size = 1)
  eset <- random_eset(n_probes)
  expect_equal(
    object = keep_all_probes(eset),
    expected = seq_len(n_probes),
    info = "keep_all_probes: returns the index along the probes (rows)"
  )
})

test_that("keep_all_samples", {
  # Input should be an ExpressionSet
  check_fails_with_non_eset("keep_all_samples")

  # Input should have some probes / samples
  check_fails_with_empty_dimensions("keep_all_samples")

  # The returned values should just be [1 : num_samples]
  n_samples <- sample.int(50, 1)
  eset <- random_eset(n_samples = n_samples)
  expect_equal(
    object = keep_all_samples(eset),
    expected = seq_len(n_samples),
    info = "keep_all_samples: returns the index along the samples (columns)"
  )
})

test_that("keep_all_entrez_probes", {
  # Input should be an ExpressionSet
  check_fails_with_non_eset("keep_all_entrez_probes")

  # Input should have some probes / samples
  check_fails_with_empty_dimensions("keep_all_entrez_probes")

  # Input should have an entrez.id column in it's featureData
  eset_empty_fdata <- random_eset()
  expect_error(
    keep_all_entrez_probes(eset_empty_fdata),
    info = paste(
      "keep_all_entrez_probes: input should have some columns in it's fData"
    )
  )

  eset_no_entrez <- random_eset(5)
  Biobase::fData(eset_no_entrez)$not_an_entrez_column <- letters[1:5]
  expect_error(
    keep_all_entrez_probes(eset_no_entrez),
    info = paste(
      "keep_all_entrez_probes: input should have an 'entrez.id' column in",
      "it's fData"
    )
  )

  # Row indexes where a non-empty / non-NA entrez.id is observed should be
  # returned
  eset_entrez <- random_eset(4, 3)
  Biobase::fData(eset_entrez)$entrez.id <- c(NA, "", "12345", "12345///786")
  expect_equal(
    keep_all_entrez_probes(eset_entrez),
    expected = c(3L, 4L),
    info = paste(
      "keep_all_entrez_probes: keep rows with non-NA and non-empty 'entrez.id'"
    )
  )
})
