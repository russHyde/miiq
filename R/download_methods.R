###############################################################################

# --- Methods on "MicroarrayDownloadConfig"

#' 'run_download_workflow' method for the MicroarrayDownloadConfig class.
#'
#' @param        object        A MicroarrayDownloadConfig object
#'
#' @include      download_classes.R
#'
#' @export
#'

run_download_workflow <- function(object) {
  # Attempts to write this as a method on MicroarrayDownloadConfig have failed
  # and I don't know why (hence it is a function, but implemented on MDC)
  stopifnot(is(object, "MicroarrayDownloadConfig"))
  downloader <- object@dl_method

  gpl_files <- download_gpl_annotations(object)

  marray_files <- downloader(acc = object@acc, dest_dir = object@dest_dir)

  append(marray_files, gpl_files)
}

#' Download GPL-annotations from GEO for a MicroarrayDownloadConfig object.
#'
#' @param        object        A MicroarrayDownloadConfig object.
#'
#' @importFrom   GEOquery      getGEO
#' @importFrom   RCurl         url.exists
#'
#' @include      download_functions.R

download_gpl_annotations <- function(object) {
  stopifnot(is(object, "MicroarrayDownloadConfig"))
  # Attempts to write this as a method on MicroarrayDownloadConfig have failed
  # and I don't know why (hence it is a function, but implemented on MDC)

  # Note that invalid GPL id kills the function.
  # Note that invalid combination of annot.gpl and gpl.acc kill the function.
  # Wrapper around getGEO for GPLs that checks that an GPL dataset is
  #   available before attempting to download it.
  #
  # Since it is passed a MicroarrayDownloadConfig
  # - the gpl_accession will be valid
  # - the dest_dir will exist
  # - the annot_gpl will be a scalar logical

  # helpers
  .gpl_url_exists <- function(x) {
    gpl_url <- .make_gpl_url(x@gpl_acc, x@annot_gpl)
    RCurl::url.exists(gpl_url)
  }
  .make_gpl_filepath <- function(x) {
    file.path(
      x@dest_dir, .make_gpl_filename(x@gpl_acc, x@annot_gpl)
    )
  }
  .get_gpl <- function(x) {
    GEOquery::getGEO(
      GEO = x@gpl_acc, destdir = x@dest_dir, AnnotGPL = x@annot_gpl
    )
  }

  # The filepath into which getGEO would store the GPL dataset
  gpl_filepath <- .make_gpl_filepath(object)

  if (file.exists(gpl_filepath)) {
    message("GPL filename ", gpl_filepath, " already exists")
  } else {
    if (!.gpl_url_exists(object)) {
      stop("GPL URL does not exist for GPL id: ", object@gpl_acc)
    }
    # Download the GPL dataset for the current microarray platform
    .get_gpl(object)
  }

  stopifnot(file.exists(gpl_filepath))

  list(
    gpl_dir = object@dest_dir,
    gpl_files = basename(gpl_filepath),
    annot_gpl = object@annot_gpl
  )
}

###############################################################################
