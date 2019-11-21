# These functions are used at the initial import of a microarray dataset
#

#' create a function that can normalise an ExpressionSet (optionally converting
#' the input dataset to an ExpressionSet beforehand)
#'
#' @param   eset_making_fn   A function to convert a dataset to an
#'   ExpressionSet. If the user is already working with an ExpressionSet, just
#'   use \code{identity} (the default). Otherwise, you might use
#'   \code{oligo::rma}.
#'
#' @param   normalise_fn   A function to normalise an ExpressionSet. Options are
#'   `quantile` (for quantile normalisation) and `none` (for do not modify; the
#'   default).
#'
#' @param   log_checking_fn   A function that can determine whether an
#'   ExpressionSet is already log-transformed. This helps prevent you from
#'   log-transforming an already log-transformed dataset, and from forgetting
#'   to log-transform and non-transformed dataset.
#'
#' @importFrom   methods       is
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @export

build_eset_normaliser <- function(
  eset_making_fn = identity,
  normalise_fn = c("none", "quantile"),
  log_checking_fn = check_if_pretransformed_eset
) {
  normalise_fn <- match.arg(normalise_fn)

  function(x) {
    # `x` is some dataset that can be converted to an ExpressionSet, or is
    # already an ExpressionSet
    # construct an ExpressionSet from the input, if it isn't already an eset
    eset <- eset_making_fn(x)
    stopifnot(methods::is(eset, "ExpressionSet"))

    is_logged <- log_checking_fn(eset)

    filter_and_transform_eset(
      eset = eset,
      log2_transform = !is_logged,
      drop_row_if_duplicated = TRUE,
      drop_row_if_zero_variance = TRUE,
      convert_inf_to_na = TRUE,
      normalise_method = normalise_fn,
      drop_row_na_inf_threshold = 0.25
    )
  }
}


###############################################################################

#' Take a processed ESet, eg from GEOquery, and normalise it
#'
#' Note this is a duplication of gld_fnDefault_transformAndProbeFilter from
#' drug.markers (but uses quantile normalisation)
#'
#' @param        eset          An ExpressionSet
#'
#' @export
#'
normalise_processed_eset <- function(eset) {
  build_eset_normaliser(identity, "quantile")(eset)
}

###############################################################################

#' Take newly imported affymetrix data and rma-normalise then apply some
#' filters to the normalised expression data
#'
#' @param        exprs         Newly imported Affymetrix data
#' @importFrom   oligo         rma
#'
#' @export
#'
normalise_raw_affy_eset <- function(exprs) {
  build_eset_normaliser(oligo::rma, "none")(exprs)
}

###############################################################################
