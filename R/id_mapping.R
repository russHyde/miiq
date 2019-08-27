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

### ======================================================================= ###
#   Gene identifier databases
### ======================================================================= ###

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
  if (missing(entrez_ids) || is.null(entrez_ids) || length(entrez_ids) == 0) {
    stop("`entrez_ids` should be defined")
  }

  if (missing(entrezgene_db) || !is_valid_orgdb(entrezgene_db)) {
    stop("`entrezgene_db` should be defined and a valid `OrgDb`")
  }

  if (
    !all(c("ENTREZID", "SYMBOL") %in% AnnotationDbi::keytypes(entrezgene_db))
  ) {
    stop("ENTREZID and SYMBOL should be keytypes in `entrezgene_db`")
  }

  # Obtain mapping between entrez gene id and gene symbol
  #   then order the mapping according to the input entrez.id vector
  #   and relabel any NA symbol values as the string "---"
  gene_map <- suppressMessages(
    AnnotationDbi::select(
      entrezgene_db,
      keys = entrez_ids, columns = "SYMBOL", keytype = "ENTREZID"
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

###############################################################################

#' Function to map gene symbols to entrez gene ids. The input gene symbols can
#' be of a variety of formats (eg, REFSEQ, ENSEMBL) provided the symbol type is
#' present as a column name in the provided AnnotationDbi dataset.
#' If any symbol maps to multiple entrez gene ids, the entrez ids are collapsed
#' into a single string.
#' The function returns a vector of entrez ids, with the order of the output
#' matching exactly with the order of the input.
#' Any symbol that can"t be mapped to an entrez id returns with the value NA
#'
#' @param        gene_symbols   vector of some type of gene symbol (eg,
#'   Ensembl, Refseq, HGNC)
#' @param        symbol_type   the type of gene symbol used (should be a
#'   colname for the entrezgene_db)
#' @param        entrezgene_db   a database object for use in mapping symbol to
#'   entrezGene id.
#' @param        collapse_character   a single character for collapsing
#'   multiple ids into a single entry.
#'
#' @return       a vector of the same length as gene_symbols containing the
#'   entrez ids that correspond to the input gene_symbols.
#'
#' @importFrom   AnnotationDbi   keytypes   keys   select
#' @export

symbol_to_entrez_id <- function(
  gene_symbols = character(0),
  symbol_type  = "SYMBOL",
  entrezgene_db = NULL,
  collapse_character = "|") {

  # Import / validate the provided database and check that genesymbols/entrez
  # ids are available
  stopifnot(is_valid_orgdb(entrezgene_db))
  stopifnot(symbol_type %in% AnnotationDbi::keytypes(entrezgene_db))

  if (is.null(gene_symbols)) {
    return(character(0))
  }

  # Genbank/refseq symbols might be input in the form "NM_001206729.1"
  #   but the entrezgene databases (eg, org.Hs.eg.db) use
  #   ids in the form "NM_001206729" - therefore, we need to
  #   strip off any trailing ".XXX" characters for genbank/refseq
  #   input
  if (symbol_type == "REFSEQ") {
    gene_symbols <- unlist(Map(
      function(x) x[1],
      strsplit(gene_symbols, "\\.")
    ))
  }

  # At least one of the gene symbols must be a valid key
  key.symbols <- intersect(
    gene_symbols,
    AnnotationDbi::keys(entrezgene_db, keytype = symbol_type)
  )

  # The function returns NA (string equivalent) for any input symbol
  # that cannot be mapped to an entrezgene id
  if (is.null(key.symbols) | length(key.symbols) == 0) {
    return(
      as.character(rep(NA, length(gene_symbols)))
    )
  }

  # Explicitly use AnnotationDbi::select since dplyr masks "select"
  eg_map <- AnnotationDbi::select(
    entrezgene_db,
    keys = unique(gene_symbols),
    keytype = symbol_type,
    columns = "ENTREZID"
  )

  if (any(duplicated(eg_map[, symbol_type]))) {
    # Gene symbols that map to multiple entrez ids return a collapsed string
    #   containing the sorted values
    entrez_ids <- sapply(
      gene_symbols,
      function(symbol) {
        r <- which(eg_map[, symbol_type] == symbol)
        paste(sort(eg_map[r, "ENTREZID"]), collapse = collapse_character)
      }
    )
    entrez_ids <- as.vector(entrez_ids)
  } else {
    # If none of the gene symbols map to multiple entries, we can return
    #   without collapsing any of the entrez ids
    reorder <- match(gene_symbols, eg_map[, symbol_type])
    entrez_ids <- as.vector(
      eg_map[reorder, "ENTREZID"]
    )
  }

  # Returned values should be NA if no valid entrez.id could be found
  invalid_ids <- c("")
  entrez_ids[which(entrez_ids %in% invalid_ids)] <- as.character(NA)

  entrez_ids
}

###############################################################################

#' multisymbol_to_entrez_ids
#'
#' A function to map a vector of gene symbols to Entrez gene ids.
#' The input may contain "multiple|joined|symbols", and if it does, the output
#' for that entry will contain a "joined|string|of|ids" containing all of the
#' unique values that any of the input symbols mapped to
#' See symbol_to_entrez_id for further details
#'
#' @param        gene_symbols   vector of gene_symbols of a single type.
#' @param        symbol_type   type of the gene_symbols (eg, Refseq).
#' @param        entrezgene_db   database for mapping gene_symbols to
#'   entrez.ids.
#' @param        split_character   character string for splitting up the
#'   symbols if any of the input entries contains more than one symbol.
#' @param        collapse_character   character for collapsing entrez ids.
#' @param        optim_flag    Can't remember.
#' @param        quiet         BOOLEAN should dplyr progress bar be suppressed
#'   (sideeffect).
#'
#' @importFrom   magrittr      %>%
#' @importFrom   dplyr         bind_rows   group_by_   do
#'
#' @export

multisymbol_to_entrez_ids <- function(
                                      gene_symbols = character(0),
                                      symbol_type = "SYMBOL",
                                      entrezgene_db = NULL,
                                      split_character = "\\|",
                                      collapse_character = "|",
                                      optim_flag = 1,
                                      quiet = TRUE) {
  c()
}
