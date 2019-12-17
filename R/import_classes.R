###############################################################################
#
##############################################################################

### ======================================================================= ###
#   Importers: Classes
### ======================================================================= ###

#' @name         meta-extdata_import_classes-NULL-function
#' @title        meta-extdata_import_classes-NULL-function
#' @description   Unnecessary description for meta-extdata_import_classes-NULL
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

#############################################################################

#' Class definition for MicroarrayImportConfig
#'
#' @name         MicroarrayImportConfig-class
#' @rdname       MicroarrayImportConfig-class
#'
#' @exportClass   MicroarrayImportConfig
#'
methods::setClass(
  Class = "MicroarrayImportConfig",
  slots = list(
    acc = "character",
    import_method = "function",
    normalise_method = "function",
    raw_dir = "character",
    raw_files = "character",
    raw_archive = "character",
    processed_dir = "character",
    processed_files = "character",
    processed_archive = "character",
    sdrf = "character",
    idf = "character",
    adf = "character",
    gpl_dir = "character",
    gpl_files = "character",
    annot_gpl = "logical"
  )
)

methods::setMethod(
  f = "initialize",
  signature = methods::signature("MicroarrayImportConfig"),
  definition = function(.Object, ...) {
    .Object <- methods::callNextMethod()
    .Object
  }
)

###############################################################################

# TODO: setValidity
# - 1st arg to import_method should be 'mic'

###############################################################################

#' Sets up a class for storing the options for importing datasets from GEO
#' and/or ArrayExpress
#'
#' @param        acc           An identifier for the microarray dataset.
#'   Typically this would be the GEO or ArrayExpress identifier.
#' @param        import_method   Either a function or the name of a function
#'   to be used for importing the microarray files into R. The function-name
#'   can either be specified in "function_name" (if the function is defined in
#'   the calling environment) or "pkg_name::function_name".
#' @param        normalise_method   Either a function or the name of a function
#'   to be used for normalising the imported microarray data. The function-name
#'   can be provided in "function_name" or "pkg_name::function_name" format.
#' @param        file_list     A list (as output by run_download_workflow) that
#'   defines the file paths for the data to be imported. When specified, none
#'   of the arguments
#' @param        raw_dir       GHI
#' @param        raw_files     JKL
#' @param        raw_archive   MNO
#' @param        processed_dir   PQR
#' @param        processed_files   STU
#' @param        processed_archive   987
#' @param        sdrf          VWX
#' @param        idf           YZ0
#' @param        adf           123
#' @param        gpl_dir       456
#' @param        gpl_files     789
#' @param        annot_gpl     Logical.
#'
#' @export
#'
MicroarrayImportConfig <- function(
                                   acc,
                                   import_method,
                                   normalise_method,
                                   file_list = NULL,
                                   raw_dir = as.character(NA),
                                   raw_files = as.character(NA),
                                   raw_archive = as.character(NA),
                                   processed_dir = as.character(NA),
                                   processed_files = as.character(NA),
                                   processed_archive = as.character(NA),
                                   sdrf = as.character(NA),
                                   idf = as.character(NA),
                                   adf = as.character(NA),
                                   gpl_dir = as.character(NA),
                                   gpl_files = as.character(NA),
                                   annot_gpl = as.logical(NA)) {
  # Ensure that the import-method is defined, and, if specified by name
  # convert that function name to the equivalent function.
  if (missing(import_method)) {
    stop("import_method should be defined")
  }

  # Ensure that the normalise-method is defined, and, if specified by name
  # convert that function name to the equivalent function.
  if (missing(normalise_method)) {
    stop("normalise_method should be defined")
  }

  .all_download_args_are_na <- all(vapply(
      list(
        raw_dir, raw_files, raw_archive,
        processed_dir, processed_files, processed_archive,
        sdrf, idf, adf,
        gpl_dir, gpl_files, annot_gpl
      ),
      is.na,
      logical(1)
  ))

  if (!is.null(file_list) && !.all_download_args_are_na) {
    stop(
      "When file_list is defined, none of the [raw|processed|gpl] file /",
      " archive / dir args should be defined"
    )
  }

  # Tried to do the following with do.call("new",
  # append(list("MicroarrayImportConfig"), args)) but it threw a C stack
  # problem
  lambda <- function(...) {
    new(
      "MicroarrayImportConfig",
      acc = acc,
      import_method = .get_from_env(import_method),
      normalise_method = .get_from_env(normalise_method),
      ...
    )
  }

  if (!is.null(file_list)) {
    do.call(lambda, file_list)
  } else {
    lambda(
      raw_dir = raw_dir,
      raw_files = raw_files,
      raw_archive = raw_archive,
      processed_dir = processed_dir,
      processed_files = processed_files,
      processed_archive = processed_archive,
      sdrf = sdrf, idf = idf, adf = adf, gpl_dir = gpl_dir,
      gpl_files = gpl_files, annot_gpl = annot_gpl
    )
  }
}

#############################################################################
