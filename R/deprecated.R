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
#' @include      microarray_classes.R   filter_functions.R
#' @export
#'

gld_fnDefault_keepSample <- gld_fnBuilder_keepSample(
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

gld_fnDefault_keepProbe <- function(
                                    # nolint start
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
