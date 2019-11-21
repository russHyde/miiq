###############################################################################

context("Test eset construction and normalisation functions")

###############################################################################

test_that(
  paste(
    "build_eset_normaliser should build a function that can construct an",
    "ExpressionSet"
  ), {
  my_normaliser <- build_eset_normaliser(
    eset_making_fn = Biobase::ExpressionSet,
    normalise_fn = "none",
    log_checking_fn = function(x) TRUE
  )

  identity_normaliser <- build_eset_normaliser(
    eset_making_fn = identity
  )

  withr::with_seed(
    42,
    my_matrix <- matrix(rnorm(20), nrow = 5, ncol = 4)
  )

  expect_silent(my_normaliser(my_matrix))

  expect_is(my_normaliser(my_matrix), "ExpressionSet")

  expect_equivalent(
    object = Biobase::exprs(my_normaliser(my_matrix)),
    expected = my_matrix
  )

  expect_error(
    identity_normaliser(my_matrix),
    info = paste(
      "`eset_making_fn` must produce an `ExpressionSet` in",
      "`build_eset_normaliser`"
    )
  )
})
