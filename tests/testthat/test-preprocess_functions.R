###############################################################################

context("Test preprocess functions")

###############################################################################

.df <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("we can add an entrez.id and symbol to an ExpressionSet", {
  eset <- Biobase::ExpressionSet(
    matrix(1:10, nrow = 5),
    featureData = Biobase::AnnotatedDataFrame(
      .df(
        `Gene ID` = as.character(1:5),
        `Gene symbol` = c("ABC///DEF", "EEE///EEE", LETTERS[3:5])
      )
    )
  )

  expect_true(
    all(
      c("entrez.id", "symbol") %in%
        colnames(
          Biobase::featureData(
            add_entrez_columns_if_gset_has_annotgpl(eset))))
  )
  expect_equal(
    Biobase::fData(add_entrez_columns_if_gset_has_annotgpl(eset))$entrez.id,
    as.character(1:5)
  )
  expect_equal(
    Biobase::fData(add_entrez_columns_if_gset_has_annotgpl(eset))$symbol,
    c("ABC|DEF", "EEE", LETTERS[3:5])
  )
})
