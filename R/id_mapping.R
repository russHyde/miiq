### ======================================================================= ###
#   Non-gene databases
### ======================================================================= ###

#' geo_id_to_accession
#'
#' Converts the numeric 9-mer UIDs used by the NCBI::GEO database into the
#' GSExxxxx / GPLxxxxxx / GSMxxxxx accession numbers that are also used by that
#' database.
#'
#' @param        geo_ids       A vector of numerical NCBI::GEO identifiers. May
#'   be given as a numeric vector or a character-vector.
#'
#' @return       A vector of NCBI::GEO accession numbers (eg, GPL12345,
#'   GSE98765, GSM24681). NA is returned for any input identifier that is not
#'   of the 9-digit 1-, 2- or 3-leading format used for uids by NCBI::GEO.
#'
#' @importFrom   magrittr      %>%
#' @importFrom   stringr       str_replace
#'
#' @export
#'

geo_id_to_accession <- function(
  geo_ids) {
  if (missing(geo_ids)) {
    stop("geo_ids is missing in input to geo_id_to_accession")
  }
  if (is.null(geo_ids) || length(geo_ids) == 0) {
    return(character(0))
  }

  # Work with the numeric ids as strings
  geo_ids_str <- as.character(geo_ids)

  # Only ids that contain 9 characters and that are solely numeric can be
  #   converted to GEO accession ids (for GSM, GPL, GSE)
  is_geo_numberstring <- function(xs) {
    # must be all-numeric 9-mer with a leading 1, 2, or 3
    grepl("^[1-3][0-9]{8}$", xs)
  }

  is_geo_numberstring(geo_ids_str) %>%
    ifelse(geo_ids_str, NA) %>%
    stringr::str_replace("^10*", "GPL") %>%
    stringr::str_replace("^20*", "GSE") %>%
    stringr::str_replace("^30*", "GSM")
}
