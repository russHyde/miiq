###############################################################################

# Differential expression analysis functions

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
