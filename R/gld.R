#' Function for filtering on microarray probes
#'
#' @param        geo.limma.dataset   An eset.limma.dataset
#' @param        gset          An ExpressionSet (overrides use of
#'   geo.limma.dataset)
#'
#' @export
#'
keep_probe_fn <- function(
  # nolint start
                          geo.limma.dataset = NULL,
                          # nolint end
                          gset = NULL) {
  # TODO: rewrite to use `geo_limma_dataset` as argname
  # TODO: unit tests

  # nolint start
  geo_limma_dataset <- geo.limma.dataset
  # nolint end

  if (is.null(gset)) {
    stopifnot(is(geo_limma_dataset, "eset_limma_dataset"))
    gset <- geo_limma_dataset@eset
  }

  stopifnot(is(gset, "ExpressionSet"))
  stopifnot(nrow(gset) > 0)

  features <- Biobase::featureData(gset)
  stopifnot("entrez.id" %in% colnames(features))

  which(
    !is.na(features$entrez.id) &
      features$entrez.id != ""
  )
}
