#' Class definition for MicroarrayPreprocessConfig
#'
#' @name         MicroarrayPreprocessConfig-class
#' @rdname       MicroarrayPreprocessConfig-class
#'
#' @exportClass   MicroarrayPreprocessConfig
#'
methods::setClass(
  Class = "MicroarrayPreprocessConfig",
  slots = list(
    acc = "character",
    annot_gpl = "logical",
    entrezgene_db = "OrgDb",
    keep_sample_fn = "function",
    keep_probe_fn = "function"
  )
)

methods::setMethod(
  f = "initialize",
  signature = methods::signature("MicroarrayPreprocessConfig"),
  definition = function(
                        .Object,
                        ...) {
    .Object <- methods::callNextMethod()
    .Object
  }
)

#' Sets up a class for storing the options for preprocessing ExpressionSets
#' (that is, filtering the samples and probes; adding entrez and symbol columns)
#'
#' @param        acc           An identifier for the microarray dataset.
#'   Typically this would be the GEO or ArrayExpress identifier.
#' @param        entrezgene_db   A string or a OrgDB defining a entrezgene
#'   database.
#' @param        annot_gpl     Logical. Have the features in the dataset been
#'   annotated with a GPL?
#' @param        keep_sample_fn,keep_probe_fn   Functions that return a vector
#'   of integer indices that can be used to subset an ExpressionSet. By default,
#'   all samples and features will be kept.
#'
#' @export
#'
MicroarrayPreprocessConfig <- function(
                                   acc,
                                   entrezgene_db,
                                   keep_sample_fn = gld_fnDefault_keepSample,
                                   keep_probe_fn = gld_fnDefault_keepProbe,
                                   annot_gpl = as.logical(NA)) {
  new(
    "MicroarrayPreprocessConfig",
    acc = acc,
    entrezgene_db = entrezgene_db,
    keep_sample_fn = keep_sample_fn,
    keep_probe_fn = keep_probe_fn,
    annot_gpl = annot_gpl
  )
}

#############################################################################
