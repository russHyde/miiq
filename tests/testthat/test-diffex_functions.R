###############################################################################

context("Test functions for differential expression analysis")

###############################################################################

test_that("limma_workflow: invalid input", {
  n_samples <- 4
  eset <- random_eset(n_probes = 50, n_samples = n_samples)
  design_fn <- function(gset) {
    matrix(1, nrow = n_samples, ncol = 1)
  }
  contrast_fn <- function(design) {
    matrix(1, nrow = 1, ncol = 1)
  }

  # gset should not be NULL
  expect_error(
    limma_workflow(NULL, design_fn, contrast_fn),
    info = "limma_workflow: gset arg should not be NULL"
  )

  # gset should be an ExpressionSet
  expect_error(
    limma_workflow("Not an ExpressionSet", design_fn, contrast_fn),
    info = "limma_workflow: gset arg should be an ExpressionSet"
  )

  # design_fn should not be NULL
  expect_error(
    limma_workflow(eset, NULL, contrast_fn),
    info = "limma_workflow: design_fn arg should not be NULL"
  )

  # design_fn should be a function
  expect_error(
    limma_workflow(eset, "Not a function", contrast_fn),
    info = "limma_workflow: design_fn arg should be a function"
  )

  # contrast_fn should not be NULL
  expect_error(
    limma_workflow(eset, design_fn, NULL),
    info = "limma_workflow: contrast_fn arg should not be NULL"
  )

  # contrast_fn should be a function
  expect_error(
    limma_workflow(eset, design_fn, "Not a function"),
    info = "limma_workflow: contrast_fn arg should be a function"
  )
})

test_that("limma_workflow: valid input", {
  # contents should include
  # - a design matrix
  # - a contrasts matrix
  # - two limma fits
  n_samples <- 4
  eset <- random_eset(n_probes = 50, n_samples = n_samples)
  design_fn <- function(gset) {
    matrix(1, nrow = n_samples, ncol = 1)
  }
  contrast_fn <- function(design) {
    matrix(1, nrow = 1, ncol = 1)
  }

  results <- expect_silent(
    limma_workflow(eset, design_fn, contrast_fn)
  )
  expect_is(
    object = results$design, class = "matrix"
  )
  expect_is(
    object = results$contrast, class = "matrix"
  )
  expect_is(
    results$fits_init, "MArrayLM"
  )
  expect_is(
    results$fits, "MArrayLM"
  )
})

test_that("`design_builder` creates a function that acts on ExpressionSets", {
  eset <- random_eset(n_samples = 12)
  Biobase::pData(eset)$group <- factor(rep(c("group1", "group2"), each = 6))
  design_function <- design_builder(
    treatment_cols = "group",
    design_fn = function(x) model.matrix(~ group, data = x)
  )
  expect_silent(design_function(eset))
  expect_error(
    design_function("Not an ExpressionSet"),
    info = "design_builder creates a function that acts on ExpressionSets"
  )
})

test_that(
  "fn returned by `design_builder` create appropriately sized designs", {
  eset <- random_eset(n_samples = 12)
  Biobase::pData(eset)$group <- factor(rep(c("A", "B"), each = 6))
  design_function <- design_builder(
    treatment_cols = "group",
    design_fn = function(x) matrix(0, nrow = 11, ncol = 2)
  )
  expect_error(
    design_function(eset),
    info = paste(
      "the design-function returned by design_builder should make a design",
      "with a row for each sample in the ExpressionSet"
    )
  )
})
