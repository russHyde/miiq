###############################################################################

# --- Generics

#' Generic for the method 'run_import_workflow'
#'
#' @docType      methods
#' @name         run_import_workflow
#' @rdname       run_import_workflow-methods
#' @export
#'
methods::setGeneric(
  "run_import_workflow",
  function(config) {
    standardGeneric("run_import_workflow")
  }
)

###############################################################################

# --- Methods

#' 'run_import_workflow' method for the MicroarrayImportConfig class
#'
#' @name         run_import_workflow
#' @rdname       run_import_workflow-methods
#' @aliases      run_import_workflow,MicroarrayImportConfig-method
#'
#' @importClassesFrom   AnnotationDbi   OrgDb
#'
#' @return       a Biobase::ExpressionSet
#'
#' @export
#'
methods::setMethod(
  f = "run_import_workflow",
  signature = methods::signature(
    config = "MicroarrayImportConfig"
  ),
  definition = function(config) {
    # Import the microarray data & any phenotypic data from the
    # locally-stored files.
    importer <- config@import_method
    initial_marray <- importer(config)

    # Background-Correct the microarray data (& create an ExpressionSet) if
    # necessary
    correcter <- config@normalise_method
    normalised_marray <- correcter(initial_marray)

    # Append feature annotations from GPL, if a GPL dataset is available
    final_marray <- if (
      !is.na(config@gpl_dir) && !is.na(config@gpl_files)
    ) {
      gpl <- import_gpl(config@gpl_files, config@gpl_dir)
      add_gpl_to_eset(normalised_marray, gpl)
    } else {
      normalised_marray
    }

    final_marray
  }
)

###############################################################################
