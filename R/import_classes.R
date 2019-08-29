###############################################################################
# 2017-07-11
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

# Class definition for MicroarrayImportConfig
#
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
  definition = function(
                        .Object,
                        ...) {
    .Object <- methods::callNextMethod()
    .Object
  }
)

###############################################################################

#' @importFrom   magrittr      %>%
#' @importFrom   Biobase       featureData
#' @importFrom   purrr         map
#' @importFrom   stringr       str_replace_all
#'
.add_entrez_columns_if_gset_has_annotgpl <- function(
                                                     gset,
                                                     # nolint start
                                                     entrezgene.db.id
                                                     # nolint end
                                                   ) {
  split_and_join <- function(xs) {
    stringr::str_replace_all(xs, "///", "|") %>%
      strsplit("\\|") %>%
      purrr::map(unique) %>%
      purrr::map(sort) %>%
      purrr::map(function(x) paste(x, collapse = "|")) %>%
      unlist()
  }

  featureData(gset)$entrez.id <- split_and_join(
    featureData(gset)$`Gene ID`
  )

  featureData(gset)$symbol <- split_and_join(
    featureData(gset)$`Gene symbol`
  )

  gset
}

###############################################################################

# TODO: setValidity
# - 1st arg to import_method should be 'mic'

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
  definition = function(
                        config) {
    # Import the microarray data & any phenotypic data from the
    # locally-stored files.
    importer <- config@import_method
    initial_marray <- importer(config)

    # Background-Correct the microarray data (& create an ExpressionSet) if
    # necessary
    correcter <- config@normalise_method
    correcter(initial_marray)
  }
)

###############################################################################

#' Generic for the method 'run_annotate_and_preprocess_workflow'
#'
#' @docType      methods
#' @name         run_annotate_and_preprocess_workflow
#' @rdname       run_annotate_and_preprocess_workflow-methods
#'
#' @export
#'

methods::setGeneric(
  "run_annotate_and_preprocess_workflow",
  function(
           eset, config, entrezgene_db_id) {
    standardGeneric("run_annotate_and_preprocess_workflow")
  }
)

#' run_annotate_and_preprocess_workflow for ExpressionSet and
#'   MicroarrayImportConfig
#'
#' @name         run_annotate_and_preprocess_workflow
#' @rdname       run_annotate_and_preprocess_workflow-methods
#' @aliases      run_annotate_and_preprocess_workflow,ExpressionSet,MicroarrayImportConfig,character-method
#'
#' @include      microarray_classes.R
#' @include      gld.R
#'
#' @export
#'

methods::setMethod(
  f = "run_annotate_and_preprocess_workflow",
  signature = methods::signature(
    eset = "ExpressionSet",
    config = "MicroarrayImportConfig",
    entrezgene_db_id = "character"
  ),
  definition = function(
                        eset,
                        config,
                        entrezgene_db_id) {
    # Annotate the features in the ExpressionSet
    gpl_annotated_marray <- if (
      !is.na(config@gpl_dir) && !is.na(config@gpl_files)
    ) {
      gpl <- import_gpl(config@gpl_files, config@gpl_dir)
      add_gpl_to_eset(eset, gpl)
    } else {
      eset
    }

    #
    entrez_adder <- if (config@annot_gpl) {
      # The gpl dataset contains Entrez Id and gene-symbol columns already
      # We only need to reformat them and add them back to the featureData
      # of the expressionSet
      .add_entrez_columns_if_gset_has_annotgpl
    } else {
      # Use the default esetAnnotation function defined in drug.markers
      # package
      gld_fnDefault_esetAnnotation
    }

    # Run preprocess_eset_workflow to add consistent entrez.id column and to
    # transform and median-centre the data, if necessary
    preprocess_eset_workflow(
      gset = gpl_annotated_marray,
      entrezgene.db.id = entrezgene_db_id,
      eset.annot.fn = entrez_adder,
      keep.sample.fn = get_keep_sample_function(config@acc),
      keep.probe.fn = keep_probe_fn
    )
  }
)

###############################################################################

#' Sets up a class for storing the options for importing datasets from GEO
#' and/or ArrayExpress
#'
#' @param        acc           ABC
#' @param        import_method   DEF
#' @param        normalise_method   FED
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
                                   raw_dir         = as.character(NA),
                                   raw_files       = as.character(NA),
                                   raw_archive     = as.character(NA),
                                   processed_dir   = as.character(NA),
                                   processed_files = as.character(NA),
                                   processed_archive = as.character(NA),
                                   sdrf            = as.character(NA),
                                   idf             = as.character(NA),
                                   adf             = as.character(NA),
                                   gpl_dir         = as.character(NA),
                                   gpl_files       = as.character(NA),
                                   annot_gpl       = as.logical(NA)) {
  if (missing(import_method)) {
    stop("import_method should be defined")
  }
  # Tried to do the following with do.call("new",
  # append(list("MicroarrayImportConfig"), args)) but it threw a C stack
  # problem
  new(
    "MicroarrayImportConfig",
    acc = acc, import_method = import_method,
    normalise_method = normalise_method, raw_dir = raw_dir,
    raw_files = raw_files, raw_archive = raw_archive,
    processed_dir = processed_dir, processed_files = processed_files,
    processed_archive = processed_archive,
    sdrf = sdrf, idf = idf, adf = adf, gpl_dir = gpl_dir,
    gpl_files = gpl_files, annot_gpl = annot_gpl
  )
}

#############################################################################
