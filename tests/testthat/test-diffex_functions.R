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

test_that(
  "fn returned by `contrast_builder` creates correct-sized contrasts", {

  # All that contrast_builder does is wrap a user-defined function in some
  # checks. For example, it checks that the dimensions of an input design
  # match those of the output contrasts matrix.

  design_2col <- data.frame(
    intercept = rep(1, 4), groupB = rep(c(0, 1), each = 2)
  )
  design_3col <- data.frame(
    intercept = rep(1, 4),
    groupB = rep(c(0, 1), each = 2),
    treatment = c(0, 1, 0, 1)
  )

  contrast_lambda <- function(design) {
    # the returned contrast matrix  has 3 rows
    # .. so it corresponds with a design with 3 columns
    matrix(1, nrow = 3, ncol = 2)
  }
  contrast_function <- contrast_builder(contrast_lambda)

  expect_error(
    contrast_function("Not a design matrix / df"),
    info = paste(
      "a built contrast function should be applied to a design df / matrix"
    )
  )

  expect_error(
    contrast_function(
      # note that the dimensions of this 'design' match those of the returned
      # contrast matrix
      data.frame(
        a = letters[1:4],
        b = 1:4,
        c = factor(letters[1:4]),
        stringsAsFactors = FALSE
      )
    ),
    info = "design data.frame should have all numeric columns"
  )

  expect_error(
    contrast_function(
      matrix(letters[1:9], nrow = 3, ncol = 3)
    ),
    info = "design matrix should be numeric"
  )

  expect_error(
    contrast_function(design),
    info = paste(
      "applying a built contrast function to a design should give a matrix",
      "with # rows == # design-columns"
    )
  )

  expect_equal(
    contrast_function(design_3col),
    contrast_lambda(design_3col),
    info = paste(
      "where the dims match up and the design is numeric, application of the",
      "built function should return the same matrix as the bare contrast",
      "lambda-function"
    )
  )
})

test_that(
  "a list of string-defined contrasts can be passed to contrast_builder", {
    design_3col <- data.frame(
      groupA = rep(c(1, 0), each = 2),
      groupB = rep(c(0, 1), each = 2),
      treatment = c(0, 1, 0, 1)
    )

    contrast_list <- list(
      delta_group = "groupB - groupA",
      treatment_effect = "treatment"
    )

    contrast_function <- contrast_builder(contrast_list)

    contrast_matrix <- matrix(
      c(c(-1, 1, 0), c(0, 0, 1)), nrow = 3, ncol = 2
    )
    colnames(contrast_matrix) <- c("delta_group", "treatment_effect")

    expect_equivalent(
      contrast_function(design_3col),
      contrast_matrix,
      info = "contrast builder can take a list of strings to define contrasts"
    )
  }
)
