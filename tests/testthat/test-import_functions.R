###############################################################################

context("Tests for import functions for external data")

###############################################################################

.df <- function(...) data.frame(..., stringsAsFactors = FALSE)

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

test_that("add_gpl_to_eset: valid inputs", {
  eset <- Biobase::ExpressionSet(
    matrix(
      1:6,
      nrow = 3,
      dimnames = list(
        paste0("G", 1:3),
        paste0("S", 1:2)
      )
    )
  )
  gpl <- new(
    "GPL",
    header = list(),
    dataTable = new(
      "GEODataTable",
      table = .df(
        ID = paste0("G", c(1, 2, 4)),
        annot = c("x", "y", "z")
      ),
      columns = .df(
        Column = c("ID", "annot"),
        Description = c("Unique identifier for the probe", "Some annotation")
      )
    )
  )

  # When all eset IDs are present in the gpl the returned eset
  # - should be the same dimension as the input;
  # - it's rownames should be in the same order;
  # - the colnames of it's featureData should match those of the GPL table
  expect_is(
    add_gpl_to_eset(eset[1:2, ], gpl), "ExpressionSet"
  )
  expect_equal(
    dim(add_gpl_to_eset(eset[1:2, ], gpl)),
    dim(eset[1:2, ]),
    info = paste(
      "If all eset IDs are present in the GPL, `add_gpl_to_eset` should not",
      "affect the eset dimensions"
    )
  )
  expect_equal(
    rownames(add_gpl_to_eset(eset[1:2, ], gpl)),
    rownames(eset[1:2, ]),
    info = paste(
      "If all eset IDs are present in the GPL, `add_gpl_to_eset` should not",
      "modify the order of the probes in the eset"
    )
  )
  expect_equal(
    colnames(featureData(add_gpl_to_eset(eset[1:2, ], gpl))),
    colnames(GEOquery::Table(gpl)),
    info = "all annotation columns in the GPL should be added to the eset"
  )

  # When some eset IDs are absent from the gpl, the returned eset
  # - should only have rows for those IDs that are present in both
  # - but have rownames in the same order as the input (modulo the missing
  # rows)
  # - and the user should be warned that not all probes have annotation data
  expect_warning(
    add_gpl_to_eset(eset, gpl),
    info = paste(
      "If the eset contains some probes that are not present in the GPL",
      "annotations, throw a warning"
    )
  )
  expect_equal(
    rownames(
      suppressWarnings(
        add_gpl_to_eset(eset, gpl)
      )
    ),
    intersect(rownames(eset), GEOquery::Table(gpl)[["ID"]]),
    info = paste(
      "If the eset contains some probes that are absent from the GPL",
      "annotations, the order of the returned probes should match that of",
      "the input eset (after removing the non-overlapping probes)"
    )
  )
})

###############################################################################
