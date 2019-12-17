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

test_that("user can pass in file_list from run_download_workflow", {

  # typical values for the file_list returned by run_download_workflow for a
  # raw GEO dataset
  #
  geo_acc <- "GSE11578"
  gpl_acc <- "GPL571"
  geo_file_list <- list(
    raw_dir = file.path("~", "ext_data", "GEOquery", "GSE11578"),
    raw_archive = "GSE11578_RAW.tar",
    processed_dir = file.path("~", "ext_data", "GEOquery"),
    processed_files = "GSE11578_series_matrix.txt.gz",
    gpl_dir = file.path("~", "ext_data", "GEOquery"),
    gpl_files = "GPL571.annot.gz",
    annot_gpl = TRUE
  )

  import_config <- expect_silent(
    MicroarrayImportConfig(
      acc = geo_acc,
      import_method = "import_geo_affy_raw",
      normalise_method = "normalise_raw_affy_eset",
      file_list = geo_file_list
    )
  )

  expect_equal(import_config@raw_dir, geo_file_list$raw_dir)
  expect_equal(import_config@raw_archive, geo_file_list$raw_archive)
  expect_equal(import_config@processed_dir, geo_file_list$processed_dir)
  expect_equal(import_config@processed_files, geo_file_list$processed_files)
  expect_equal(import_config@gpl_dir, geo_file_list$gpl_dir)
  expect_equal(import_config@gpl_files, geo_file_list$gpl_files)
  expect_equal(import_config@annot_gpl, geo_file_list$annot_gpl)

  expect_error(
    MicroarrayImportConfig(
      acc = geo_acc,
      import_method = "import_geo_affy_raw",
      normalise_method = "normalise_raw_affy_eset",
      file_list = geo_file_list,
      raw_dir = "some_directory"
    ),
    info = paste(
      "when `file_list` is specified, [raw|processed|gpl]_[dir|files|archive]",
      "and annot_gpl should not be specified"
    )
  )
})
