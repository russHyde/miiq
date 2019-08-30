###############################################################################
# 2017-07-11
#
###############################################################################

context("Tests for import functions for external data")

###############################################################################

#' @importFrom   Biobase       ExpressionSet
#'
test_that("add_gpl_to_eset: invalid inputs", {
  es <- Biobase::ExpressionSet()
  expect_error(
    object = add_gpl_to_eset("Not an ExpressionSet"),
    info = "Input to add_gpl_to_eset should include an ExpressionSet/ESet"
  )

  expect_error(
    object = add_gpl_to_eset(
      es, "Not a GPL"
    ),
    info = "Input to add_gpl_to_eset should include a GPL"
  )
})

###############################################################################
