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


###############################################################################

#' Builds a function that makes the design matrix for a dataset
#'
#' User passes in the column labels of the phenoData that are to be used in
#'   making the design
#' And passes in a design-function that converts these columns into the design
#' The functionBuilder checks for the validity of the input (treatment_cols /
#' design_fn) and the built function checks the validity of the ExpressionSet
#' that gets passed in to it.
#'
#' @param        treatment_cols   Which columns of the pData should be kept?
#' @param        design_fn     Function for making design matrices.
#'
#' @importFrom   Biobase       pData   sampleNames
#'
#' @include      utils.R   diffex_functions.R
#' @export
#'

design_builder <- function(treatment_cols = NULL, design_fn = NULL) {
  # check type-validity of the treatment_cols and design_fn
  stopifnot(is.character(treatment_cols) && length(treatment_cols) > 0)
  stopifnot(is.function(design_fn))

  design_function <- function(gset = NULL) {
    # check input validity of the Expressionset / geo_limma_dataset
    if (is.null(gset) || !methods::is(gset, "ExpressionSet")) {
      stop("'gset' should be an 'ExpressionSet'")
    }

    # Extract the columns of pData that contain the design-related info
    # and then construct the design from these columns
    pdata <- Biobase::pData(gset)
    if (!all(treatment_cols %in% colnames(pdata))) {
      missing_cols <- setdiff(treatment_cols, colnames(pdata))
      stop(paste(
        "Some treatment_cols were missing from pData:",
        paste(missing_cols, collapse = "/")
      ))
    }
    treatments <- pdata[, treatment_cols, drop = FALSE]
    design <- as.data.frame(design_fn(treatments))
    sn <- Biobase::sampleNames(gset)
    cn <- make.names(colnames(design))

    if (length(sn) != nrow(design)) {
      message(
        c("sampleNames (ExpressionSet): ", paste(sn, collapse = " "))
      )
      message(
        c("sampleNames (design): ", paste(rownames(design), collapse = " "))
      )
      stop("length of sampleNames (in eset) and design do not agree")
    }

    rownames(design) <- sn
    colnames(design) <- cn
    design
  }

  design_function
}

###############################################################################

#' Builds a function that makes the contrast matrix for a dataset
#'
#' Note that for a function-as-input, all this does is wrap that function and
#' check the dimensions and type of the input design against the output
#' contrast matrix.
#'
#' @param    x    Either a closure or a list of string-based contrast
#'   definitions. If a closure, this should accept a design matrix and return
#'   a contrast matrix.
#'
#' @return    A function that accepts a design matrix / data.frame and
#'   constructs a contrasts matrix from it.
#' @importFrom   limma   makeContrasts
#' @export

contrast_builder <- function(x) {

  # convert a list of string-defined contrasts into a function
  f <- if(is.function(x)) {
    x
  } else if (is.list(x)) {
    function(design) {
      cont_mat <- limma::makeContrasts(
        contrasts = x,
        levels = design
      )
      if (!is.null(names(x))) {
        colnames(cont_mat) <- names(x)
      }
      cont_mat
    }
  }

  contrast_function <- function(design) {
    stopifnot(
      (
        is(design, "data.frame") && is.numeric(as.matrix(design))
      ) || (
        is(design, "matrix") && is.numeric(design)
      )
    )
    contrast_matrix <- f(design)

    stopifnot(nrow(contrast_matrix) == ncol(design))
    contrast_matrix
  }

  contrast_function
}
