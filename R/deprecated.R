# nolint start

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
#' @include      utils.R
#' @export
#'

# nolint start
gld_fnBuilder_keepSample <- function(
                                     column.filter.fn
                                     # nolint end
) {
  stopifnot(is.function(column.filter.fn))
  fn <- column.filter.fn

  function(
             # nolint start
             geo.limma.dataset = NULL,
             # nolint end
             gset = NULL) {
    warning("deprecation warning: recommend using a function(gset)")
    gset <- .check_or_get_eset(geo.limma.dataset, gset)
    keep_cols <- fn(gset)
    keep_cols
  }
}

###############################################################################
#' Function for filtering on microarray probes
#'
#' @param        geo.limma.dataset   An eset.limma.dataset
#' @param        gset          An ExpressionSet (overrides use of
#'   geo.limma.dataset)
#'
#' @include      utils.R   filter_functions.R
#' @export
#'
keep_probe_fn <- function(
                          # nolint start
                          geo.limma.dataset = NULL,
                          # nolint end
                          gset = NULL) {
  warning(
    "`miiq::keep_probe_fn` is deprecated, please use",
    "`miiq::keep_all_entrez_probes`"
  )

  gset <- .check_or_get_eset(geo.limma.dataset, gset)

  keep_all_entrez_probes(gset)
}

###############################################################################
#' Default function for choosing which samples to keep in a ESet / g-l-d
#'
#' @param        geo.limma.dataset   An eset_limma_dataset.
#' @param        gset          An ExpressionSet - used preferentially
#'   over ELD@eset.
#'
#' @importFrom   Biobase       sampleNames
#' @include      filter_functions.R
#' @export
#'

# nolint start
gld_fnDefault_keepSample <- gld_fnBuilder_keepSample(
  # nolint end
  column.filter.fn = keep_all_samples
)

###############################################################################
#' Default function for choosing which probes to keep in an ESet / g-l-d
#'
#' @param        geo.limma.dataset   A geo_limma_dataset
#' @param        gset          An ExpressionSet
#'
#' @importFrom   Biobase       featureNames
#'
#' @include      utils.R   filter_functions.R
#' @export
#'

# nolint start
gld_fnDefault_keepProbe <- function(
                                    geo.limma.dataset = NULL,
                                    # nolint end
                                    gset = NULL) {
  warning(
    "`miiq::gld_fnDefault_keepProbe` is deprecated, please use",
    "`miiq::keep_all_probes`"
  )
  # TODO: replace argnames with the lint-passing `geo_limma_dataset`
  gset <- .check_or_get_eset(geo.limma.dataset, gset)

  keep_all_probes(gset)
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
#' @include      utils.R   diffex_functions.R
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
  warning(
    "Creating design-makers using `miiq::gld_fnBuilder_exptDesign` is",
    "deprecated. Please use `miiq::design_builder` instead."
  )
    design_builder(
      treatment_cols = treatment.cols, design_fn = design.fn
    )(
      gset = .check_or_get_eset(geo.limma.dataset, gset)
    )
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
    title <- treatments[, 1]
    model.matrix(~title)
  }
)

# nolint end

###############################################################################
