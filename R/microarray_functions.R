###############################################################################

#' Function that performs a few tests to check if a given \code{ExpressionSet}
#' has already been log-transformed
#'
#' @param        eset          an \code{ExpressionSet}.
#' @param        range_limit_if_transformed   the largest possible range
#' between high and low values if the \code{ESet} is already transformed.
#' @param        skew_threshold   if the skewness of all columns of the
#' dataset is greater than this threshold, the dataset is assumed
#' untransformed.
#'
#' @include      vector_manipulation.R
#'
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @importFrom   Biobase      exprs
#' @importFrom   moments      skewness
#'
#' @export
#'
check_if_pretransformed_eset <- function(
                                         eset = NULL,
                                         range_limit_if_transformed = 1000,
                                         skew_threshold = 3) {
  # NOTE: refactored 20180223 - argnames are now underscore-separated
  # - This function is used by rnaseq_marray_merge and 2015_prelim, but those
  # projects have not yet been refactored to use the renamed args.
  # - If this script is imported back into rnaseq_marray_merge or 2015_prelim,
  # rewrite all calls to this function

  # Perform a few tests to check if a given eset has been log transformed:

  if (!exists("na_inf_omit")) {
    # how do I ensure that na_inf_omit is imported?
    stop("na_inf_omit is required by check_if_pretransformed")
  }

  # Check that the input is a valid 'eset'
  if (is.null(eset) || !is(eset, "ExpressionSet")) {
    stop("'eset' should be a valid Biobase::ExpressionSet")
  }

  # If there is no data in the eset, or all entries are NA,
  #   we consider the eset to have been transformed
  if (nrow(eset) == 0 || ncol(eset) == 0 || all(is.na(Biobase::exprs(eset)))) {
    return(TRUE)
  }
  # Determine if any negative values are present;
  #   assume the dataset has been log-transformed if it contains any negatives
  if (any(as.vector(Biobase::exprs(eset)) < 0, na.rm = TRUE)) {
    return(TRUE)
  }
  # Determine the range of expression values in the dataset
  rng <- range(
    na_inf_omit(Biobase::exprs(eset))
  )

  # Determine the skewness of each microarray in the dataset
  skews <- apply(Biobase::exprs(eset), 2, function(x) {
    moments::skewness(na_inf_omit(x))
  })

  # If the dataset is not transformed it will have high range,
  #   and will have skewed distribution of expression values in each array
  # If either of these occur, we assume the dataset is not transformed
  untransformed <- (rng[2] - rng[1]) > range_limit_if_transformed ||
    all(skews > skew_threshold)

  return(!untransformed)
}
