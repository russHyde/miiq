### ======================================================================= ###
#   Diffex analysis: Classes
### ======================================================================= ###

#' @name         meta-diffex_classes-NULL-function
#' @title        meta-diffex_classes-NULL-function
#' @description   Unnecessary description for meta-diffex_classes-NULL
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

###############################################################################

#' Class definition for DiffexConfig
#'
#' @name         DiffexConfig-class
#' @rdname       DiffexConfig-class
#'
#' @exportClass   DiffexConfig
#'

methods::setClass(
  Class = "DiffexConfig",
  slots = list(
    design_fn = "function",
    contrast_fn = "function"
  )
)

methods::setMethod(
  f = "initialize",
  signature = methods::signature("DiffexConfig"),
  definition = function(.Object, ...) {
    .Object <- methods::callNextMethod()
    .Object
  }
)

###############################################################################

#' Constructor function for DiffexConfig class
#'
#' @param   design   Either a function or a design matrix. If it is a function,
#' then design(ExpressionSet) should be a design matrix.
#' @param   contrast   Either a function or a contrast matrix. If it is a
#' function the contrast(design_matrix) should be a contrast matrix.
#'
#' @export

DiffexConfig <- function(design, contrast) {
  # Allow user to specify design and contrast as either functions or as matrix
  #
  stopifnot(is.function(design) || is.matrix(design) || is.data.frame(design))
  stopifnot(is.function(contrast) || is.matrix(contrast))

  # If a design or contrast has been specified as a matrix, make a function
  # that will return that matrix
  .fn_switch <- function(x) {
    if (is.function(x)) {
      x
    } else {
      function(...) x
    }
  }

  design_fn <- .fn_switch(design)
  contrast_fn <- .fn_switch(contrast)

  new("DiffexConfig", design_fn = design_fn, contrast_fn = contrast_fn)
}
