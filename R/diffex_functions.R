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

  # TODO: fix this; these introduce bugs when used with an unnamed gset (as
  # first argument) and
  # gld_fnBuilder_[exptDesign|exptContrasts|exptContrasts_fromList] because,
  # the latter expect a eset_limma_dataset as the first argument
  # - Suggest,
  #    - deprecate gld_fnBuilder_...;
  #    - add design_builder(treatment_cols, design_fn)(gset)
  #    - add contrast_builder(contrast_fn)(design)
  #    - add contrast_builder_from_list(contrast_list)(design)
  # - That is, rewrite them without any reference to geo_limma_dataset object
  # in the arg list
  design <- design_fn(gset = gset)

  contrast <- contrast_fn(design = design)

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
