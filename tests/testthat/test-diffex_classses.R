###############################################################################

context("DiffexConfig tests")

# DiffexConfig() should return
# - an S4 object of class DiffexConfig
# - with a `design_fn` entry
# - with a `contrast_fn` entry
# - where if the user submits a function as the `design` argument, then the
#   `design_fn` of the DiffexConfig should be the same function
# - and similarly, if the user submits a function as `contrast`, then the
#   `contrast_fn` of the returned object should be that function
#
# The arguments to DiffexConfig() should be
# - `design` - either a design matrix (or data-frame) or a function that
#   constructs a design matrix when applied to an ExpressionSet
# - `contrast` - either a contrast matrix, or a function that constructs a
#   contrast matrix when applied to a design matrix

###############################################################################

test_that("Construction of a DiffexConfig using design / contrast functions", {

  my_design_fn <- function(gset) {
    matrix(
      1, nrow = ncol(gset), ncol = 1,
      dimnames = list(colnames(gset), "intercept")
    )
  }
  my_contrast_fn <- function(design) {
    matrix(
      1, nrow = 1, ncol = 1,
      dimnames = list("intercept", "baseline_expression")
    )
  }

  expect_is(
    DiffexConfig(design = my_design_fn, contrast = my_contrast_fn),
    "DiffexConfig"
  )
  expect_equal(
    DiffexConfig(design = my_design_fn, contrast = my_contrast_fn)@design_fn,
    my_design_fn
  )
  expect_equal(
    DiffexConfig(design = my_design_fn, contrast = my_contrast_fn)@contrast_fn,
    my_contrast_fn
  )

  # Validity of input to DiffexConfig()
  expect_error(
    DiffexConfig(design = "Not a fn, matrix or df", contrast = my_contrast_fn),
    info = paste(
      "`design` argument of `DiffexConfig()` should be a function, or matrix",
      "/ data-frame."
    )
  )
  expect_error(
    DiffexConfig(design = my_design_fn, contrast = "Not a fn, matrix or df"),
    info = paste(
      "`contrast` argument of `DiffexConfig()` should be a function or matrix."
    )
  )
})

test_that("Construction of a DiffexConfig using design / contrast matrices", {
  eset <- random_eset(n_samples = 8)
  some_factor <- factor(rep(letters[1:2], each = 4))

  design_matrix <- model.matrix(~ -1 + some_factor)
  colnames(design_matrix) <- levels(some_factor)

  contrast_matrix <- limma::makeContrasts(
    b_vs_a = "b - a", levels = design_matrix
  )

  expect_is(
    DiffexConfig(design = design_matrix, contrast = contrast_matrix),
    "DiffexConfig"
  )

  dc <- DiffexConfig(design = design_matrix, contrast = contrast_matrix)

  expect_equal(
    dc@design_fn(eset),
    design_matrix,
    info = paste(
      "Calling the `design_fn()` of a `DiffexConfig` should return the same",
      "design-matrix as was used to construct that `DiffexConfig` (if a",
      "design-matrix was used)."
    )
  )
  expect_equal(
    dc@contrast_fn(design_matrix),
    contrast_matrix,
    info = paste(
      "If a contrast-matrix was used in `DiffexConfig()`, that matrix should",
      "be returned when calling the `contrast_fn()` of the `DiffexConfig`",
      "object."
    )
  )
})
