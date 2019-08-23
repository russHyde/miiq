###############################################################################

#' @name         nullFunction
#' @title        nullFunction
#' @importFrom   methods       setValidity
#' @importFrom   methods       setClass
#' @importFrom   methods       setGeneric
#' @importFrom   methods       setMethod
#' @importFrom   methods       setReplaceMethod
NULL

###############################################################################

#' validity function for eset_limma_dataset objects
#'
#' @importFrom   magrittr      %>%
#'
#' @noRd
#'

.validity_eld <- function(object) {
    valid <- TRUE
    msgs <- c()

    # Check that the number of probes is consistent between eset, fits and
    #   fits.init
    check_probes <- function(L) {
      non_empty_nrows <- Map(
          nrow,
          Filter(function(x) nrow(x) > 0, L)
        ) %>%
        unlist()

      if (length(non_empty_nrows) > 0) {
        return(all(non_empty_nrows == non_empty_nrows[1]))
      } else {
        return(TRUE)
      }
    }

    if (
      !check_probes(list(object@eset, object@fits.init, object@fits))
    ) {
      valid <- FALSE
      msgs <- c(
        msgs,
        paste(
          "All non-empty entries in {eset, fits, fits.init} should have",
          "the same number of probes (rows)"
        )
      )
    }

    # Check that the number of samples is consistent between eset and design
    check_samples <- function(
                              eset, design) {
      if (ncol(eset) == 0 || nrow(design) == 0) {
        return(TRUE)
      }
      return(ncol(eset) == nrow(design))
    }

    if (!check_samples(object@eset, object@design)) {
      valid <- FALSE
      msgs <- c(
        msgs,
        paste(
          "All non-empty entries in {eset, design} should have the same",
          "number of samples"
        )
      )
    }

    #
    if (!valid) {
      return(msgs)
    } else {
      return(TRUE)
    }
  }

#' @title        Class to hold ExpressionSets and limma model results
#'
#' @importClassesFrom   limma   MArrayLM
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @name         eset_limma_dataset-class
#' @rdname       eset_limma_dataset-class
#'
#' @exportClass   eset_limma_dataset
#'
methods::setClass(
  "eset_limma_dataset",
  slots = list(
    eset = "ExpressionSet",
    design = "data.frame",
    contrast = "matrix",
    fits.init = "MArrayLM",
    fits = "MArrayLM"
  ),
  prototype = list(
    eset = Biobase::ExpressionSet(),
    design = data.frame(),
    contrast = matrix(nrow = 0, ncol = 0),
    fits.init = methods::new("MArrayLM"),
    fits = methods::new("MArrayLM")
  ),
  validity = .validity_eld
)

#' wrapper function eset_limma_dataset
#'
#' @name         eset_limma_dataset
#' @rdname       eset_limma_dataset-class
#'
#' @param        ...           arguments passed to the constructor for
#'   `eset_limma_dataset`.
#'
#' @importFrom   methods       new
#'
#' @export
#'
eset_limma_dataset <- function(...) {
  methods::new("eset_limma_dataset", ...)
}

###############################################################################

# eset_limma_dataset methods:

#' Row-subsetting operator for `eset_limma_dataset` class
#'
#' @param        x             An `eset_limma_dataset` object
#' @param        i             A subset of the rows in `x`
#' @param        j,drop,...    Not used
#'
#' @name         [
#' @rdname       eset_limma_dataset
#' @aliases      [,eset_limma_dataset,numeric,missing,missing-method

methods::setMethod(
  "[",
  signature = c(
    x = "eset_limma_dataset",
    i = "numeric",
    j = "missing",
    drop = "missing"
  ),
  definition = function(
                        x, i, j, ..., drop) {
    row_subsetter <- function(
                              eset_or_fit) {
      if (nrow(eset_or_fit) > 0) {
        return(eset_or_fit[i, ])
      } else {
        return(eset_or_fit)
      }
    }

    # eset, fits and fits.init should be restrictied to the requested rows
    .eset <- row_subsetter(x@eset)
    .fits <- row_subsetter(x@fits)
    .fits_init <- row_subsetter(x@fits.init)

    # design & contrast are unaffected by subsetting the number of probes
    .design <- x@design
    .contrast <- x@contrast

    eset_limma_dataset(
      eset = .eset,
      design = .design,
      contrast = .contrast,
      fits = .fits,
      fits.init = .fits_init
    )
  }
)


###############################################################################

#' If no geo_limma_dataset is provided, or the provided one has no eset,
#' returns the gset
#'
#' Checks that the gset is an ExpressionSet
#'
#' @importFrom   methods       is
#'
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @noRd
#'

.check_or_get_eset <- function(geo_limma_dataset = NULL, gset = NULL) {
  if (is.null(gset)) {
    stopifnot(methods::is(geo_limma_dataset, "eset_limma_dataset"))
    gset <- geo_limma_dataset@eset
  }
  stopifnot(methods::is(gset, "ExpressionSet"))
  gset
}

###############################################################################

## #' @import       GEOquery
## #' @import       Biobase
## #' @import       limma

###############################################################################
