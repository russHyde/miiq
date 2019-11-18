###############################################################################

context("Unit tests for classes for importing external microarray datasets")

###############################################################################

#'
test_that("MicroarrayImportConfig: constructor", {

  # Using new() constuctor, no args provided
  expect_s4_class(
    object = new("MicroarrayImportConfig"),
    class = "MicroarrayImportConfig"
  )

  # Using constructor function, with all args provided
  expect_s4_class(
    object = MicroarrayImportConfig(
      acc = "GSE123",
      import_method = import_geo_affy_raw,
      normalise_method = identity
    ),
    class = "MicroarrayImportConfig"
  )

  expect_error(
    object = MicroarrayImportConfig(
      import_method = import_geo_affy_raw,
      normalise_method = identity
    ),
    info = "MicroarrayImportConfig() must have 'acc' defined"
  )
  expect_error(
    object = MicroarrayImportConfig(
      acc = "GSE123",
      normalise_method = identity
    ),
    info = "MicroarrayImportConfig() must have import_method defined"
  )
  expect_error(
    object = MicroarrayImportConfig(
      acc = "GSE123",
      import_method = import_geo_affy_raw
    ),
    info = "MicroarrayImportConfig() must have normalise_method defined"
  )

  expect_equal(
    object = MicroarrayImportConfig(
      acc = "GSE1234",
      import_method = "miiq::import_geo_processed",
      normalise_method = identity
    )@import_method,
    expected = miiq::import_geo_processed,
    info = "Import method can be specified using pkg::fn syntax"
  )

  expect_equal(
    object = MicroarrayImportConfig(
      acc = "GSE1234",
      import_method = "miiq::import_geo_processed",
      normalise_method = "base::identity"
    )@normalise_method,
    expected = identity,
    info = "Normalise method can be specified using pkg::fn syntax"
  )
})

###############################################################################
