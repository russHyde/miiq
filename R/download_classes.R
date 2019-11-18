### ======================================================================= ###
#   Downloaders: Classes
### ======================================================================= ###

#' @name         meta-extdata_download_classes-NULL-function
#' @title        meta-extdata_download_classes-NULL-function
#' @description   Unnecessary description for
#'   meta-extdata_download_classes-NULL
#'
#' @importFrom   methods       callNextMethod
#' @importFrom   methods       new
#' @importFrom   methods       setClass
#' @importFrom   methods       setGeneric
#' @importFrom   methods       setMethod
#' @importFrom   methods       setValidity
#' @importFrom   methods       signature
#'
NULL

###############################################################################

#' Check the validity of a Microarray-Download-Config object
#'
#' @param        object        A potential MicroarrayDownloadConfig object.
#'
#' @include      download_functions.R
#'
#' @return       Either TRUE or a character-vector of errors.
#' @noRd
#'

.check_marray_dl_config <- function(
                                    object) {
  cls <- "MicroarrayDownloadConfig"

  err_msg <- function(x) {
    paste0(cls, ": ", x)
  }

  db <- object@database
  if (length(db) != 1) {
    return(err_msg("`dl_database` should be of length 1"))
  }
  if (!(db %in% c("geo", "aryx"))) {
    return(err_msg("`dl_database` should be one of `geo` or `aryx`"))
  }

  acc <- object@acc
  if (length(acc) != 1) {
    return(err_msg("`dl_acc` should be of length 1"))
  }
  if (db == "geo" && !.is_gse_acc(acc)) {
    return(
      err_msg("`dl_acc` should be a 'GSExxx' ID is `dl_database` is `geo`")
    )
  }
  if (db == "aryx" && !.is_aryx_acc(acc)) {
    return(
      err_msg("`dl_acc` should be `E-MTAB-xxxx` if `dl_database` is `aryx`")
    )
  }

  # Tests on dl_method
  # - that it is a function is tested at construction
  # - Check that dl_method takes an accession number * directory as args
  # - acc and dest_dir should be the only args
  dl_method <- object@dl_method
  if (!all.equal(c("acc", "dest_dir"), names(formals(dl_method))[1:2])) {
    return(
      err_msg("`acc` and `dest_dir` should be the first args to `dl_method`")
    )
  }

  dest_dir <- object@dest_dir
  if (length(dest_dir) != 1 || !dir.exists(dest_dir)) {
    return(
      err_msg("`dl_dest_dir` should be of length 1 and an existing dir")
    )
  }

  # TODO: Allow NA gpl values for raw-data downloads?
  gpl <- object@gpl_acc
  if (length(gpl) != 1 || !.is_gpl_acc(gpl)) {
    return(
      err_msg("`gpl_acc` should be of length 1 and a GPLxxx accession number")
    )
  }

  annot_gpl <- object@annot_gpl
  if (length(annot_gpl) != 1) {
    return(
      err_msg("`annot_gpl` should be of length 1")
    )
  }

  TRUE
}

###############################################################################

#' Class definition for MicroarrayDownloadConfig
#'
#' @name         MicroarrayDownloadConfig-class
#' @rdname       MicroarrayDownloadConfig-class
#'
#' @exportClass   MicroarrayDownloadConfig
#'

methods::setClass(
  Class = "MicroarrayDownloadConfig",
  slots = list(
    acc = "character",
    database = "character",
    dl_method = "function",
    dest_dir = "character",
    gpl_acc = "character",
    annot_gpl = "logical"
    # TODO: raw_files, raw_archive, processed_files, samples_file
  )
)

methods::setMethod(
  f = "initialize",
  signature = methods::signature("MicroarrayDownloadConfig"),
  definition = function(.Object, ...) methods::callNextMethod()
)

# Pass-through validity checker for MicroarrayDownloadConfig, see
#   .check_marray_dl_config for details.
#
methods::setValidity(
  Class = "MicroarrayDownloadConfig",
  method = function(object) .check_marray_dl_config(object)
)

###############################################################################

#' Generic for the method 'run_download_workflow'
#'
#' @docType      methods
#' @name         run_download_workflow
#' @rdname       run_download_workflow-methods
#'
#' @export
#'
methods::setGeneric(
  name = "run_download_workflow",
  def = function(object) standardGeneric("run_download_workflow")
)

###############################################################################

#' 'run_download_workflow' method for the MicroarrayDownloadConfig class.
#'
#' @name         run_download_workflow
#' @rdname       run_download_workflow-methods
#' @aliases      run_download_workflow,MicroarrayDownloadConfig-method
#'
#' @export
#'
methods::setMethod(
  f = "run_download_workflow",
  signature = methods::signature(
    object = "MicroarrayDownloadConfig"
  ),
  definition = function(
                          object) {
    downloader <- object@dl_method

    gpl_files <- if (!is.na(object@gpl_acc)) {
      dl_gpl_annotations(
        gpl_acc = object@gpl_acc,
        dest_dir = object@dest_dir,
        annot_gpl = object@annot_gpl
      )
    } else {
      list()
    }

    marray_files <- downloader(
      acc = object@acc,
      dest_dir = object@dest_dir
    )

    append(marray_files, gpl_files)
  }
)

###############################################################################

#' Sets up a class for storing the options for downloading datasets from GEO
#' and/or ArrayExpress
#'
#' @param        acc           An accession number for use at GEO or
#'   ArrayExpress.
#' @param        database      The name of the database from which the dataset
#'   should be downloaded.
#' @param        download_method   Either a function (or the name of a
#'   function) for downloading the required dataset.
#' @param        dest_dir      The directory within which the dataset should be
#'   saved.
#' @param        gpl_acc       The accession number for any relevant GPL
#'   annotation data from GEO.
#' @param        annot_gpl     Should the NCBI-embellished platform annotation
#'   data be downloaded?
#'
#' @include      utils.R
#'
#' @export
#'

MicroarrayDownloadConfig <- function(
                                     acc,
                                     database,
                                     download_method,
                                     dest_dir = tempdir(),
                                     gpl_acc = as.character(NA),
                                     annot_gpl) {
  if (is.character(download_method)) {
    download_method <- .get_from_env(download_method)
  }

  new(
    "MicroarrayDownloadConfig",
    acc = acc,
    database = database,
    dl_method = download_method,
    dest_dir = dest_dir,
    gpl_acc = gpl_acc,
    annot_gpl = annot_gpl
  )
}

###############################################################################
