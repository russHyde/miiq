###############################################################################

# Differential expression analysis functions

###############################################################################

#' limma workflow: .....
#'
#' @param        gset          An ExpressionSet
#' @param        design_fn   A function that creates a design matrix based on
#'   gset as input.
#' @param        contrast_fn   A function that creates a contrasts matrix
#'   based on an input design matrix.
#'
#' @return       A list(design, contrast, fits_init, fits).
#'
#' @importFrom   limma         lmFit   contrasts.fit   eBayes
#'
#' @export

limma_workflow <- function(
                           gset = NULL,
                           design_fn = NULL,
                           contrast_fn = NULL) {
  # Check that validity of the passed in ESet and functions
  if (
    is.null(gset) || !is(gset, "ExpressionSet")
  ) {
    stop("gset should be an ExpressionSet in limma_workflow")
  }

  if (
    is.null(design_fn) || !is.function(design_fn)
  ) {
    stop(paste(
      "design_fn should be a function(ExpressionSet => data.frame) in",
      "limma_workflow"
    ))
  }

  if (
    is.null(contrast_fn) || !is.function(contrast_fn)
  ) {
    stop(paste(
      "contrast_fn should be a function(design_df => contrasts_matrix) in",
      "limma_workflow"
    ))
  }

  # Construct design matrix and contrasts matrix using the ExpressionSet /
  #   functions
  design <- design_fn(gset)

  contrast <- contrast_fn(design)

  # Fit coefficients, contrasts and then regularise the contrasts using
  # `limma`
  fits_init <- limma::lmFit(
    object = gset,
    design = design
  )

  fits_cont <- limma::contrasts.fit(
    fit = fits_init,
    contrasts = contrast
  )

  fits <- limma::eBayes(
    fit = fits_cont
  )

  list(
    "design" = design,
    "contrast" = contrast,
    "fits_init" = fits_init,
    "fits" = fits
  )
}
