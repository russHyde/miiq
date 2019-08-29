###############################################################################
# 2017-07-11
#
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
})

###############################################################################
