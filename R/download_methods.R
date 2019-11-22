###############################################################################

# --- Generics

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

# --- Methods on "MicroarrayDownloadConfig"

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
