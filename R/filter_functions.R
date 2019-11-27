#' keep_all_probes: return the index for each probe on an ExpressionSet
#'
#' @param   gset  An \code{ExpressionSet}.
#' @importClassesFrom   Biobase   ExpressionSet
#' @export

keep_all_probes <- function(gset) {
  stopifnot(is(gset, "ExpressionSet"))
  stopifnot(nrow(gset) > 0 && ncol(gset) > 0)
  seq_len(nrow(gset))
}

#' keep_all_entrez_probes
#'
#' @inheritParams   keep_all_probes
#' @importFrom   Biobase   featureData
#' @export

keep_all_entrez_probes <- function(gset) {
  stopifnot(is(gset, "ExpressionSet"))
  stopifnot(nrow(gset) > 0 && ncol(gset) > 0)
  stopifnot("entrez.id" %in% colnames(Biobase::featureData(gset)))
  features <- Biobase::featureData(gset)
  which(
    !is.na(features$entrez.id) & features$entrez.id != ""
  )
}

#' keep_all_samples
#'
#' @inheritParams   keep_all_probes
#' @importClassesFrom   Biobase   ExpressionSet
#' @export

keep_all_samples <- function(gset) {
  stopifnot(is(gset, "ExpressionSet"))
  stopifnot(nrow(gset) > 0 && ncol(gset) > 0)
  seq_len(ncol(gset))
}

###############################################################################

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
  warning("`keep_probe_fn` is deprecated, please use `keep_all_entrez_probes`")

  # TODO: rewrite to use `geo_limma_dataset` as argname

  # nolint start
  geo_limma_dataset <- geo.limma.dataset
  # nolint end

  if (is.null(gset)) {
    stopifnot(is(geo_limma_dataset, "eset_limma_dataset"))
    gset <- geo_limma_dataset@eset
  }

  keep_all_entrez_probes(gset)
}
