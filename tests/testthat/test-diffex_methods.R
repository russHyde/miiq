###############################################################################

context(
  "Tests for methods on differential-expression-related classes (DiffexConfig)"
)

###############################################################################

test_that("run_diffex_workflow: invalid input", {
  n_samples <- 4
  eset <- random_eset(n_probes = 50, n_samples = n_samples)
  design_fn <- function(gset) {
    matrix(1, nrow = n_samples, ncol = 1)
  }
  contrast_fn <- function(design) {
    matrix(1, nrow = 1, ncol = 1)
  }

  diffex_config <- DiffexConfig(design_fn, contrast_fn)

  # eset shouldn't be NULL
  expect_error(
    run_diffex_workflow(NULL, diffex_config),
    info = "eset should not be NULL in `run_diffex_workflow`"
  )

  # eset should be an ExpressionSet
  expect_error(
    run_diffex_workflow("Not an ExpressionSet", diffex_config),
    info = "eset should be an ExpressionSet in `run_diffex_workflow`"
  )

  # config shouldn't be NULL
  expect_error(
    run_diffex_workflow(eset, NULL),
    info = "config should not be NULL in `run_diffex_workflow`"
  )

  # config should be a DiffexConfig
  expect_error(
    run_diffex_workflow(eset, "Not a DiffexConfig"),
    info = "config should be a DiffexConfig in `run_diffex_workflow`"
  )

  # design should have a row for each sample (column) in the ESet
  wrong_sized_design_for_the_eset <- function(gset) {
    matrix(1, nrow = n_samples + 1, ncol = 1)
  }
  config_with_wrong_sized_design <- DiffexConfig(
    wrong_sized_design_for_the_eset, contrast_fn
  )
  expect_error(
    run_diffex_workflow(eset, config_with_wrong_sized_design),
    info = "design made by design_fn should have a row for each column in eset"
  )

  # contrast matrix should have a row for each design-coefficient (column of
  # the design matrix)
  wrong_sized_contrast_for_the_design <- function(design) {
    matrix(1, nrow = 2, ncol = 1)
  }
  config_with_wrong_sized_contrast <- DiffexConfig(
    design_fn, wrong_sized_contrast_for_the_design
  )
  expect_error(
    run_diffex_workflow(eset, config_with_wrong_sized_contrast),
    info = paste(
      "contrast made by contrast_fn should have a row for each column of the",
      "design"
    )
  )
})

test_that("run_diffex_workflow: valid input", {
  # contents of the results should include
  # - a design matrix
  # - a contrasts matrix
  # - two limma fits
  n_samples <- 4

  eset <- random_eset(n_probes = 50, n_samples = n_samples)
  design <- matrix(1, nrow = n_samples, ncol = 1)
  contrast <- matrix(1, nrow = 1, ncol = 1)
  diffex_config <- DiffexConfig(design, contrast)

  results <- expect_silent(
    run_diffex_workflow(eset, diffex_config)
  )
  expect_equal(
    object = results$design, design
  )
  expect_equal(
    object = results$contrast, contrast
  )
  expect_is(
    results$fits_init, "MArrayLM"
  )
  expect_is(
    results$fits, "MArrayLM"
  )
})
