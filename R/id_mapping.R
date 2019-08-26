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

###############################################################################

is_valid_orgdb <- function(entrezgene_db) {
  # this hollow function makes it easier to test / mock within
  # `paste_gene_symbols`

  if (is.null(entrezgene_db)) {
    stop("No `entrezgene_db` was defined")
  }

  # db can either be an `OrgDb` object, or the name of an `OrgDb` object
  # in the latter case, the OrgDb should be installed

  # - this currently doesn't verify that the installed package corresponds to
  # an OrgDb
  is(entrezgene_db, "OrgDb")
}

#' Function to zip together a vector of entrez_ids with the gene symbol
#'   corresponding to each. The returned vector is in the form
#'   c("10000|AKT3", "1234|CCR5", ...)
#'
#' @param        entrez_ids    a vector of NCBI entrez gene ids.
#'
#' @param        entrezgene_db   an AnnotationDbi-style database for use in
#'   obtaining the symbols corresponding to the gene ids.
#'
#' @param        collapse_character   a single character for use in collapsing
#'   gene.id and gene.symbol into "gene.id|gene.symbol".
#'
#' @return       a vector of "gene.id|gene.symbol" entries.
#'
#' @importFrom   AnnotationDbi   keytypes   select
#'
#' @export

paste_gene_symbols <- function(entrez_ids,
                               entrezgene_db,
                               collapse_character = "|") {
  if(missing(entrez_ids) || is.null(entrez_ids) || length(entrez_ids) == 0) {
    stop("`entrez_ids` should be defined")
  }

  if(missing(entrezgene_db) || !is_valid_orgdb(entrezgene_db)) {
    stop("`entrezgene_db` should be defined and a valid `OrgDb`")
  }

  if(
    !all(c("ENTREZID", "SYMBOL") %in% AnnotationDbi::keytypes(entrezgene_db))
  ) {
    stop("ENTREZID and SYMBOL should be keytypes in `entrezgene_db`")
  }

  # Obtain mapping between entrez gene id and gene symbol
  #   then order the mapping according to the input entrez.id vector
  #   and relabel any NA symbol values as the string "---"
  gene_map <- suppressMessages(
    AnnotationDbi::select(
      entrezgene_db, keys = entrez_ids, columns = "SYMBOL", keytype = "ENTREZID"
    )
  )
  symbols <- gene_map[match(entrez_ids, gene_map[, "ENTREZID"]), "SYMBOL"]

  cleaned_map <- data.frame(
    entrez_id = entrez_ids,
    symbol = ifelse(is.na(symbols), "---", symbols),
    stringsAsFactors = FALSE
  )
  paste(cleaned_map$entrez_id, cleaned_map$symbol, sep = "|")
}
