test_that(".check_or_get_eset", {
  eset <- random_eset()

  # with neither ESet nor GLD
  expect_error(
    object = .check_or_get_eset(),
    info = "both eset and gld are NULL"
  )

  # with invalid ESet
  expect_error(
    object = .check_or_get_eset(gset = "NOT AN EXPRESSION SET"),
    info = "gset should be an expression set in .check_or_get_eset"
  )

  # with invalid GLD
  expect_error(
    object = .check_or_get_eset(geo_limma_dataset = "NOT A GLD"),
    info = paste(
      "geo_limma_dataset should be valid if provided to .check_or_get_eset"
    )
  )

  # with GLD
  expect_equal(
    object = .check_or_get_eset(
      geo_limma_dataset = eset_limma_dataset(eset = eset)
    ),
    expected = eset,
    info = ".check_or_get_eset should return the eset from a geo_limma_dataset"
  )

  # with ESet
  expect_equal(
    object = .check_or_get_eset(
      gset = eset
    ),
    expected = eset,
    info = ".check_or_get_eset should return the unaltered gset argument"
  )
})
