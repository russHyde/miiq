###############################################################################

context("Tests for Microarray download methods")

###############################################################################

mdc_lambda <- function(download_list = list(), dest_dir = tempdir()) {
  MicroarrayDownloadConfig(
    acc = "GSE12345",
    database = "geo",
    download_method = function(acc, dest_dir) download_list,
    dest_dir = dest_dir,
    gpl_acc = "GPL98765",
    annot_gpl = FALSE
  )
}

###############################################################################

test_that("input / output for `run_download_workflow()`", {

  # `run_download_workflow` runs on a MicroarrayDownloadConfig

  expect_error(
    run_download_workflow("Not a MicroarrayDownloadConfig"),
    info = "run_download_workflow requires a MicroarrayDownloadConfig"
  )

  # run_download_workflow returns the microarray-files and the GPL-annotation
  # files for the platform

  gse_downloads <- list(raw_file = "some_gse_filename")
  gpl_results <- list(
    gpl_dir = "some_dir", gpl_files = "some_gpl_filename", annot_gpl = TRUE
  )
  mdc <- mdc_lambda(download_list = gse_downloads)
  mockery::stub(
    run_download_workflow, "download_gpl_annotations", gpl_results
  )
  expect_equal(
    run_download_workflow(mdc),
    append(gse_downloads, gpl_results),
    info = "call to run_download_workflow()"
  )

})

###############################################################################

test_that("mocked tests for `download_gpl_annotations()`", {

  # the method is only defined for MicroarrayDownloadConfig objects

  expect_error(
    download_gpl_annotations("Not a MicroarrayDownloadConfig"),
    info = "download_gpl_annotations requires a MicroarrayDownloadConfig"
  )

  # expect an error if the GPL URL does not exist

  mdc <- mdc_lambda(dest_dir = tempdir())
  mockery::stub(download_gpl_annotations, "RCurl::url.exists", FALSE)
  expect_error(
    download_gpl_annotations(mdc),
    info = "URL for the GPL annotations should exist",
    regexp = "GPL URL does not exist"
  )

  # if the gpl file already exists, it's path is returned

  td <- tempfile(pattern = "gplTest")
  dir.create(td)
  gpl_file <- "GPL98765.tsv"
  file.create(file.path(td, gpl_file))
  mdc <- mdc_lambda(dest_dir = td)
  mockery::stub(download_gpl_annotations, ".make_gpl_filename", gpl_file)
  expect_equal(
    download_gpl_annotations(mdc),
    list(
      gpl_dir = td, gpl_files = gpl_file, annot_gpl = mdc@annot_gpl
    ),
    info = "If a GPL file already exists, it is returned"
  )

  # download_gpl_annotations should return the directory, and basename for
  # the downloaded files and whether an annot_gpl was accessed

  td <- tempfile(pattern = "gplTest")
  dir.create(td)
  gpl_file <- "GPL98765.tsv"
  mdc <- mdc_lambda(dest_dir = td)
  mockery::stub(download_gpl_annotations, "RCurl::url.exists", TRUE)
  mockery::stub(download_gpl_annotations, ".make_gpl_filename", gpl_file)
  mockery::stub(
    download_gpl_annotations, "GEOquery::getGEO", function(...) {
      file.create(file.path(td, gpl_file))
    }
  )
  expect_equal(
    download_gpl_annotations(mdc),
    list(
      gpl_dir = td, gpl_files = gpl_file, annot_gpl = mdc@annot_gpl
    ),
    info = "If the GPL file does not exist locally, it is obtained with getGEO"
  )

  # expect an error if, after attempting to download the data, the GPL file
  # does not exist
})

###############################################################################
