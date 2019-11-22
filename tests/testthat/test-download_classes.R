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
      acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    class = "MicroarrayDownloadConfig"
  )

  expect_equal(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
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

  valid_arguments <- list(
    acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
    dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
  )
  invalid_arguments <- list(
    # dest dir must be a real dir
    list(dest_dir = "NOT A DIRECTORY"),
    # only a single database is allowed
    list(database = c("geo", "geo")),
    # database must be either geo or aryx
    list(database = "NOT A DATABASE"),
    # only a single accession is allowed
    list(acc = c("GSE1234", "GSE9876")),
    # not every function is a valid download-function
    list(download_method = "c"),
    # for ArrayExpress, the Accession number should be "E-MTAB-..."
    list(database = "aryx", acc = "NOT AN ARRAY EXPRESS ACC")
  )

  for (replacements in invalid_arguments) {
    test_args <- valid_arguments
    test_args[names(replacements)] <- replacements

    expect_error(
      object = do.call(
        MicroarrayDownloadConfig, test_args
      ),
      info = paste(
        "invalid arguments in `MicroarrayDownloadConfig`: ", replacements
      )
    )
  }
})

###############################################################################

test_that("MicroarrayDownloadConfig: validation", {

  # If database is 'geo' then accession-number should be a GSE accession:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "not-a-gse-acc", database = "geo", download_method = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "If database is geo, then accession should be GSE*"
  )

  # If gpl_acc is given, it should be a GPL* accession number:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "not-a-gpl", annot_gpl = TRUE
    ),
    info = "If gpl_acc is given it should be of the form GPL*"
  )

  # If gpl_acc is given, it should be a GPL* accession number:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "not_a_function",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "download_method should correspond to a defined function"
  )

  # annot_gpl should be a single logical value:
  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = "not_logical"
    ),
    info = "annot_gpl should be a logical value"
  )

  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "dl_geo_raw",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = c(TRUE, FALSE)
    ),
    info = "annot_gpl should be a single logical value"
  )

  expect_error(
    object = MicroarrayDownloadConfig(
      acc = "GSE123", database = "geo", download_method = "identity",
      dest_dir = tempdir(), gpl_acc = "GPL987", annot_gpl = TRUE
    ),
    info = "args to the function refered by download_method should take 'acc'
and 'dest_dir' as first two args.",
    regexp = "`acc` and `dest_dir` should be the first args to `dl_method`"
  )
})

###############################################################################
