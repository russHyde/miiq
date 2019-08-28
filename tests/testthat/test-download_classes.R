###############################################################################

context("Tests for classes for downloading external microarray datasets")

###############################################################################

test_that("MicroarrayDownloadConfig: constructor", {

  # Using new() constuctor, no args provided
  expect_s4_class(
    object = new("MicroarrayDownloadConfig"),
    class = "MicroarrayDownloadConfig"
  )

  # Using constructor function, with all args provided
  expect_s4_class(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    class = "MicroarrayDownloadConfig"
  )

  expect_equal(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    expected = new(
      "MicroarrayDownloadConfig",
      acc = "GSE123", database = "geo", dl_method = dl_geo_raw,
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "MDC returned by MicroarrayDownloadConfig constructor should be the
same as using new() constructor. Modulo function-name vs function."
  )

  # TODO: check default values - gpl_acc   = as.character(NA)
  #                            - annot_gpl = as.logical(NA)
})

###############################################################################

test_that("MicroarrayDownloadConfig: validation", {

  # If database is 'geo' then accession-number should be a GSE accession:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "not-a-gse-acc", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "If database is geo, then accession should be GSE*"
  )

  # If gpl_acc is given, it should be a GPL* accession number:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "not-a-gpl", annot_gpl = TRUE
    ),
    info = "If gpl_acc is given it should be of the form GPL*"
  )

  # If gpl_acc is given, it should be a GPL* accession number:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "not_a_function",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "dl_func_name should correspond to a defined function"
  )

  # annot_gpl should be a single logical value:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = "not_logical"
    ),
    info = "annot_gpl should be a logical value"
  )

  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = c(TRUE, FALSE)
    ),
    info = "annot_gpl should be a single logical value"
  )

  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", dl_func_name = "identity",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "args to the function refered by dl_func_name should take 'acc'
and 'dest_dir' as first two args."
  )
})

###############################################################################
