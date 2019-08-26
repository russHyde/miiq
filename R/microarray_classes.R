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
#' @noRd
#'

.validity_eld <- function(object) {
  valid <- TRUE
  msgs <- c()

  # Check that the number of probes is consistent between eset, fits and
  #   fits.init
  check_probes <- function(x) {
    non_empty_nrows <- unlist(
      Map(
        nrow,
        Filter(function(x) nrow(x) > 0, x)
      )
    )

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
  check_samples <- function(eset, design) {
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
#' @importFrom   methods       new
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
  validity = function(object) .validity_eld(object)
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
  definition = function(x, i, j, ..., drop) {
    row_subsetter <- function(eset_or_fit) {
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

.check_or_get_eset <- function(
                               geo_limma_dataset = NULL,
                               gset = NULL) {
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

#' Builds a function that returns the col indices of the desired samples
#'
#' Creates a function that can be applied to a geo-limma-dataset
#'   or to an ExpressionSet, and returns the column indices of those
#'   samples that are to kept
#' This just reduces the boilerplate involved in checking g-l-d or
#'   ExpressionSet nature of the input
#'
#' @param        column.filter.fn   Function that decides which samples to keep
#'
#' @export
#'
gld_fnBuilder_keepSample <- function(
                                     # nolint start
                                     column.filter.fn
                                     # nolint end
) {
  # nolint start
  stopifnot(is.function(column.filter.fn))
  fn <- column.filter.fn
  # nolint end

  function(
             # nolint start
             geo.limma.dataset = NULL,
             # nolint end
             gset = NULL) {
    # nolint start
    gld <- geo.limma.dataset
    # nolint end
    gset <- .check_or_get_eset(gld, gset)
    keep_cols <- fn(gset)
    keep_cols
  }
}

###############################################################################

#' Default function for choosing which samples to keep in a ESet / g-l-d
#'
#' @param        geo.limma.dataset   An eset_limma_dataset.
#' @param        gset          An ExpressionSet - used preferentially
#'   over ELD@eset.
#'
#' @importFrom   Biobase       sampleNames
#'
#' @export
#'
gld_fnDefault_keepSample <- gld_fnBuilder_keepSample(
  # keep all samples
  # crash if no samples exist
  column.filter.fn = function(gset) {
    sn <- Biobase::sampleNames(gset)
    stopifnot(length(sn) > 0)
    seq_along(sn)
  }
)

###############################################################################

#' Default function for choosing which probes to keep in an ESet / g-l-d
#'
#' @param        geo.limma.dataset   A geo_limma_dataset
#' @param        gset          An ExpressionSet
#'
#' @importFrom   Biobase       featureNames
#'
#' @export
#'
gld_fnDefault_keepProbe <- function(
                                    # nolint start
                                    geo.limma.dataset = NULL,
                                    # nolint end
                                    gset = NULL) {
  # nolint start
  # TODO: replace argnames with the lint-passing `geo_limma_dataset`
  geo_limma_dataset <- geo.limma.dataset
  # nolint end

  gset <- .check_or_get_eset(geo_limma_dataset, gset)

  features <- Biobase::featureNames(gset)

  stopifnot(length(features) > 0)

  keep_rows <- seq_along(features)
  return(keep_rows)
}

###############################################################################

#' Builds a function that makes the design matrix for a dataset
#'
#' User passes in the column labels of the phenoData that are to be used in
#'   making the design
#' And passes in a design.function that converts these columns into the design
#' The functionBuilder checks for the validity of the input (treatment.cols /
#'   design.function) and the built function checks the validity of the
#'   ExpressionSet / geo-limma-dataset that was passed in
#'
#' @param        treatment.cols   Which columns of the pData should be kept?
#' @param        design.fn     Function for making design matrices.
#'
#' @importFrom   Biobase       pData   sampleNames
#'
#' @export
#'
gld_fnBuilder_exptDesign <- function(
                                     treatment.cols = NULL,
                                     design.fn = NULL) {
  # check type-validity of the treatment.cols and design.fn
  stopifnot(is.character(treatment.cols) && length(treatment.cols) > 0)
  stopifnot(is.function(design.fn))

  design_function <- function(
                                # nolint start
                                geo.limma.dataset = NULL,
                                # nolint end
                                gset = NULL) {
    # nolint start
    geo_limma_dataset <- geo.limma.dataset
    # nolint end
    # check input validity of the Expressionset / geo_limma_dataset
    gset <- .check_or_get_eset(geo_limma_dataset, gset)

    # Extract the columns of pData that contain the design-related info
    # and then construct the design from these columns
    pdata <- Biobase::pData(gset)
    if (!all(treatment.cols %in% colnames(pdata))) {
      missing_cols <- setdiff(treatment.cols, colnames(pdata))
      stop(paste(
        "Some treatment.cols were missing from pData:",
        paste(missing_cols, collapse = "/")
      ))
    }
    treatments <- Biobase::pData(gset)[, treatment.cols]
    design <- as.data.frame(design.fn(treatments))
    sn <- Biobase::sampleNames(gset)
    cn <- make.names(colnames(design))

    if (
      length(sn) != nrow(design) ||
        length(cn) != ncol(design)
    ) {
      message(c("colnames: ", paste(cn, collapse = " ")))
      message(c("sampleNames: ", paste(sn, collapse = " ")))
      stop("length of colNames/sampleNames and design size do not agree")
    }

    rownames(design) <- sn
    colnames(design) <- cn
    design
  }

  design_function
}

###############################################################################

#' Default function for setting up an expt design from the pData of the
#'   ExpressionSet
#'
#' Just takes the "title" column of the pData and uses it as a factor.
#'
#' @param        geo.limma.dataset   A geo_limma_dataset.
#' @param        gset          An ExpressionSet.
#'
#' @export
#'
gld_fnDefault_exptDesign <- gld_fnBuilder_exptDesign(
  treatment.cols = "title",
  design.fn = function(treatments) {
    title <- treatments
    model.matrix(~title)
  }
)

###############################################################################

# TODO: add tests / classes / functions:
# - geo_config, geo_limma_dataset, limma_config
# - decideTests_gld, eset<-,
# - limma_workflow, preprocess_eset_workflow
# - gld_fnDefault_esetAnnotation
# - gld_fnDefault_transformAndProbeFilter
# - gld_fnBuilder_exptContrasts
# - gld_fnBuilder_exptContrasts_fromList

# TODO:
# - keep class definitions, setters / getters / subsetters into
#   this file
# - move functions, pipelines etc into different files
