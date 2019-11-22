#' Generic for the method 'run_preprocess_workflow'
#'
#' @docType      methods
#' @name         run_preprocess_workflow
#' @rdname       run_preprocess_workflow-methods
#'
#' @export
#'
methods::setGeneric(
  "run_preprocess_workflow",
  function(eset, config) {
    standardGeneric("run_preprocess_workflow")
  }
)

# nolint start

#' 'run_preprocess_workflow' for ExpressionSet and MicroarrayPreprocessConfig
#'
#' @name         run_preprocess_workflow
#' @rdname       run_preprocess_workflow-methods
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @aliases      run_preprocess_workflow,ExpressionSet,MicroarrayPreprocessConfig-method
#'
#' @return       a Biobase::ExpressionSet
#'
#' @export

# nolint end

methods::setMethod(
  f = "run_preprocess_workflow",
  signature = methods::signature(
    eset = "ExpressionSet",
    config = "MicroarrayPreprocessConfig"
  ),
  definition = function(
    eset,
    config
  ) {
    entrez_adder <- if (config@annot_gpl) {
      # The gpl dataset contains Entrez Id and gene-symbol columns already
      # We only need to reformat them and add them back to the featureData
      # of the expressionSet
      add_entrez_columns_if_gset_has_annotgpl
    } else {
      # Use the default esetAnnotation function defined in drug.markers
      # package
      gld_fnDefault_esetAnnotation
    }

    # Run preprocess_eset_workflow to add consistent entrez.id column and to
    # transform and median-centre the data, if necessary
    preprocess_eset_workflow(
      gset = eset,
      entrezgene_db = config@entrezgene_db,
      eset_annot_fn = entrez_adder,
      keep_sample_fn = config@keep_sample_fn,
      keep_probe_fn = config@keep_probe_fn
    )
  }
)
