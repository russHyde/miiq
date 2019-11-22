#' Reformat Entrez ID and Gene Symbol columns
#'
#' @importFrom   magrittr      %>%
#' @importFrom   Biobase       featureData
#' @importFrom   purrr         map
#' @importFrom   stringr       str_replace_all
#'
#' @noRd
add_entrez_columns_if_gset_has_annotgpl <- function(gset, ...) {
  split_and_join <- function(xs) {
    stringr::str_replace_all(xs, "///", "|") %>%
      strsplit("\\|") %>%
      purrr::map(unique) %>%
      purrr::map(sort) %>%
      purrr::map(function(x) paste(x, collapse = "|")) %>%
      unlist()
  }

  Biobase::featureData(gset)$entrez.id <- split_and_join(
    Biobase::featureData(gset)$`Gene ID`
  )

  Biobase::featureData(gset)$symbol <- split_and_join(
    Biobase::featureData(gset)$`Gene symbol`
  )

  gset
}
