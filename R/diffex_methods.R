###############################################################################

#' Generic for the method 'run_diffex_workflow'
#'
#' @docType      methods
#' @name         run_diffex_workflow
#' @rdname       run_diffex_workflow-methods
#'
#' @export
#'
methods::setGeneric(
  "run_diffex_workflow",
  function(eset, config) {
    standardGeneric("run_diffex_workflow")
  }
)

#' 'run_diffex_workflow' for ExpressionSet and DiffexConfig
#'
#' @name         run_diffex_workflow
#' @rdname       run_diffex_workflow-methods
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @aliases      run_diffex_workflow,ExpressionSet,DiffexConfig-method
#'
#' @return       a list with entries (design, contrast, fits_init, fits)
#'
#' @include      diffex_classes.R
#'
#' @importFrom   limma         lmFit   contrasts.fit   eBayes
#' @export

# nolint end

methods::setMethod(
  f = "run_diffex_workflow",
  signature = methods::signature(
    eset = "ExpressionSet",
    config = "DiffexConfig"
  ),
  definition = function(
                          eset,
                          config) {
    # Extract experimental design and required comparisons from the config
    design <- config@design_fn(gset = eset)
    contrast <- config@contrast_fn(design = design)

    # Fit coefficients, contrasts and then regularise the contrasts using
    # `limma`
    fits_init <- limma::lmFit(object = eset, design = design)

    fits_cont <- limma::contrasts.fit(fit = fits_init, contrasts = contrast)

    fits <- limma::eBayes(fit = fits_cont)

    list(
      "design" = design,
      "contrast" = contrast,
      "fits_init" = fits_init,
      "fits" = fits
    )
  }
)

###############################################################################
