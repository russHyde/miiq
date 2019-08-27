###############################################################################

context("Tests mappings to/from Database-specific Ids")

###############################################################################

# #' @import       org.Hs.eg.db
# ensembl_fn <- function(x) symbol_to_entrez_id(
#     gene_symbols = x,
#     symbol_type = "ENSEMBL",
#     entrezgene_db = org.Hs.eg.db::org.Hs.eg.db
#   )
# optim.flag <- 1 # used for multisymbol_to_entrez_id

# ###############################################################################

### ======================================================================= ###
#   Non-gene databases
### ======================================================================= ###

test_that("geo_id_to_accession: Conversion of GEO identifiers", {
  # GSE, GDS, GSM, GPL

  expect_error(
    object = geo_id_to_accession(),
    info = "No input to geo_id_to_accession"
  )
  expect_equal(
    object = geo_id_to_accession(character(0)),
    expected = character(0),
    info = "Empty input to geo_id_to_accession"
  )
  expect_equal(
    object = geo_id_to_accession(NULL),
    expected = character(0),
    info = "NULL input to geo_id_to_accession"
  )
  expect_equal(
    object = geo_id_to_accession("200008833"),
    expected = "GSE8833",
    info = "Single string-input to geo_id_to_accession"
  )
  expect_equal(
    object = geo_id_to_accession(200001234),
    expected = "GSE1234",
    info = "Single numeric-input to geo_id_to_accession"
  )
  expect_equal(
    object = geo_id_to_accession("NOT_A_NUM"),
    expected = as.character(NA),
    info = "Strings that aren't representations of numbers should return NA"
  )
  expect_equal(
    object = geo_id_to_accession("999888777"),
    expected = as.character(NA),
    info = "9-mer numeric ids should start with 1, 2, or 3"
  )
  expect_equal(
    object = geo_id_to_accession(c("1222333444", "0122233344")),
    expected = as.character(c(NA, NA)),
    info = "Number strings must be exactly 9-long and start with 1, 2, or 3"
  )
  expect_equal(
    object = geo_id_to_accession(
      c(
        "200008833", "200014671", "100000570", "302450495", "1962", "4278",
        "NOT_A_NUMBER_STRING", "LENGTH_9_"
      )
    ),
    expected = c(
      "GSE8833", "GSE14671", "GPL570", "GSM2450495", NA, NA, NA,
      NA
    ),
    info = paste(
      "Various string inputs to geo_id_to_accession; valid GSE/GPL/GSM are",
      "mapped to accession numbers; valid GDS ids are not mapped to accession",
      "numbers; incorrectly formatted input is mapped to NA"
    )
  )
})

###############################################################################

### ======================================================================= ###
#   Gene-identifier databases
### ======================================================================= ###

test_that("Mappings from SYMBOL to Entrez ids", {
  test_df <- data.frame(
    ENTREZID = c("10000", "1234", "54407", "9876", "8765"),
    SYMBOL = c("AKT3", "CCR5", "SLC38A2", "DUPLICATED", "DUPLICATED"),
    stringsAsFactors = FALSE
  )

  symbol_fn <- function(x) {
    filtered_df <- test_df[test_df$SYMBOL %in% x, ]
    mockery::stub(symbol_to_entrez_id, "is_valid_orgdb", TRUE)
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keys", unique(test_df$SYMBOL))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keytypes", colnames(test_df))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::select", filtered_df)

    symbol_to_entrez_id(
      gene_symbols = x,
      symbol_type = "SYMBOL",
      entrezgene_db = "SOME.DATABASE"
    )
  }

  expect_equal(
    object = symbol_fn(c()),
    expected = character(0),
    info = "Empty vector input"
  )
  expect_equal(
    object = symbol_fn(NULL),
    expected = character(0),
    info = "NULL input"
  )
  expect_equal(
    object = symbol_fn(""),
    expected = as.character(NA),
    info = "Empty string input"
  )
  expect_equal(
    object = symbol_fn("AKT3"),
    expected = "10000",
    info = "Single, mappable, symbol"
  )
  expect_equal(
    object = symbol_fn(c("AKT3", "CCR5", "SLC38A2")),
    expected = c("10000", "1234", "54407"),
    info = "Multiple, mappable, symbols"
  )
  expect_equal(
    object = symbol_fn(c("CCR5", "AKT3", "SLC38A2")),
    expected = c("1234", "10000", "54407"),
    info = "Multiple, mappable, symbols - rearranged order"
  )
  expect_equal(
    object = symbol_fn("KTELC1"),
    expected = as.character(NA),
    info = "Single, unmappable, symbol"
  )
  expect_equal(
    object = symbol_fn(c("KTELC1", "AKT3")),
    expected = c(NA, "10000"),
    info = "Unmappable and mappable symbols together"
  )
  expect_equal(
    object = symbol_fn("DUPLICATED"),
    expected = "8765|9876",
    info = "Symbol that maps to multiple Entrez IDs"
  )
})

# ###############################################################################

test_that("Mappings from Genbank id to Entrez.id", {
  test_df <- data.frame(
    ENTREZID = c("10000", "10000"),
    REFSEQ = c("NM_001206729", "NM_005465"),
    stringsAsFactors = FALSE
  )

  genbank_fn <- function(x, sym = "REFSEQ") {
    trimmed_x <- gsub(pattern = "^(.*)\\.(.*)$", replacement = "\\1", x)
    filtered_df <- test_df[test_df[, sym] %in% trimmed_x, ]
    mockery::stub(symbol_to_entrez_id, "is_valid_orgdb", TRUE)
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keys", unique(test_df[, sym]))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keytypes", colnames(test_df))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::select", filtered_df)

    symbol_to_entrez_id(
      gene_symbols = x,
      symbol_type = sym,
      entrezgene_db = "SOME.DATABASE"
    )
  }

  expect_equal(
    object = genbank_fn(c()),
    expected = character(0),
    info = "Empty genbank input"
  )
  expect_equal(
    object = genbank_fn(NULL),
    expected = character(0),
    info = "NULL genbank input"
  )
  expect_equal(
    object = genbank_fn(""),
    expected = as.character(NA),
    info = "Empty-string genbank input"
  )
  expect_equal(
    object = genbank_fn("NM_001206729"),
    expected = "10000",
    info = "Single valid genbank input"
  )
  expect_equal(
    object = genbank_fn(c("NM_001206729", "NM_005465")),
    expected = c("10000", "10000"),
    info = "Two valid genbank inputs, map to the same gene"
  )
  expect_equal(
    object = genbank_fn("NM_001206729.1"),
    expected = "10000",
    info = "Genbank input with NM_mainID.subID both present"
  )
})

# ###############################################################################

# test_that("Different specifications of EntrezGene database", {
#   expect_error(
#     object = symbol_to_entrez_id(
#       gene_symbols = c("NM_001206729", "NM_005465"),
#       symbol_type = "REFSEQ",
#       entrezgene_db = "not.an.eg.db"
#     ),
#     info = "Invalid entrezgene database"
#   )
#   expect_equal(
#     object = symbol_to_entrez_id(
#       c("NM_001206729", "NM_005465"),
#       symbol_type = "REFSEQ",
#       entrezgene_db = "org.Hs.eg.db"
#     ),
#     expected = c("10000", "10000"),
#     info = "Entrezgene DB passed as a string"
#   )
# })

###############################################################################

test_that("Ensembl-Genes mapping to EntrezGene", {
  test_df <- data.frame(
    ENTREZID = c("7982", "93655", ""),
    ENSEMBL = c("ENSG00000004866", "ENSG00000004866", "ENSG00000067601"),
    stringsAsFactors = FALSE
  )

  ensembl_fn <- function(x, sym = "ENSEMBL") {
    filtered_df <- test_df[test_df[, sym] %in% x, ]
    mockery::stub(symbol_to_entrez_id, "is_valid_orgdb", TRUE)
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keys", unique(test_df[, sym]))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::keytypes", colnames(test_df))
    mockery::stub(symbol_to_entrez_id, "AnnotationDbi::select", filtered_df)

    symbol_to_entrez_id(
      gene_symbols = x,
      symbol_type = sym,
      entrezgene_db = "SOME.DATABASE"
    )
  }

  expect_equal(
    object = ensembl_fn("ENSG00000004866"),
    expected = "7982|93655",
    info = "Multiply-mapping input"
  )
  expect_equal(
    object = ensembl_fn(""),
    expected = as.character(NA),
    info = "Empty string Ensembl input"
  )
  expect_equal(
    object = ensembl_fn("ENSG00000067601"),
    expected = as.character(NA),
    info = paste(
      "Processed psedogene: initially mapped to empty string instead of <NA>"
    )
  )
})

# ###############################################################################

# test_that("Testing multisymbol mapper: Ensembl to Entrez", {
#   lambda_ensembl <- function(x) {
#     multisymbol_to_entrez_ids(
#       gene_symbols = x,
#       symbol_type = "ENSEMBL",
#       entrezgene_db = org.Hs.eg.db::org.Hs.eg.db,
#       optim_flag = optim.flag
#     )
#   }

#   expect_equal(
#     object = lambda_ensembl("ENSG00000004866"),
#     expected = c("7982|93655"),
#     info = "Multiply-mapping single gene"
#   )

#   expect_equal(
#     object = lambda_ensembl("ENSG00000004866|ENSG0000004866"),
#     expected = c("7982|93655"),
#     info = "Multiply-mapping duplicated gene"
#   )
# })

# ###############################################################################

# test_that("Testing multisymbol mapper: Genbank to Entrez", {
#   lambda_fn_genbank <- function(x) {
#     multisymbol_to_entrez_ids(
#       gene_symbols = x,
#       symbol_type = "REFSEQ",
#       entrezgene_db = org.Hs.eg.db::org.Hs.eg.db
#     )
#   }

#   expect_equal(
#     object = lambda_fn_genbank(c()),
#     expected = c(),
#     info = "Empty input(genbank)"
#   )

#   expect_equal(
#     object = lambda_fn_genbank(NULL),
#     expected = c(),
#     info = "NULL input(genbank)"
#   )

#   expect_equal(
#     object = lambda_fn_genbank(""),
#     expected = as.character(NA),
#     info = "Empty string input(genbank)"
#   )

#   expect_equal(
#     object = lambda_fn_genbank("DOESNOTMAP"),
#     expected = as.character(NA),
#     info = "Single non-mapping non-empty string(genbank)"
#   )

#   expect_equal(
#     object = lambda_fn_genbank("NM_001206729"),
#     expected = "10000",
#     info = "Singly-mapping input"
#   )
# })

###############################################################################
test_that("Testing entrez_id -> 'entrez_id|gene_symbol' mappings", {
  test_df <- data.frame(
    ENTREZID = c(
      "4", "1234", "10000", "1111", "1111", "4563", "7504", "956530"
    ),
    SYMBOL = c(
      NA, "CCR5", "AKT3", "SYMBOL1", "SYMBOL2", NA, "XK", "NA"
    ),
    stringsAsFactors = FALSE
  )

  build_test <- function(df = test_df, valid_db = TRUE) {
    function(...) {
      mockery::stub(paste_gene_symbols, "is_valid_orgdb", valid_db)
      mockery::stub(paste_gene_symbols, "AnnotationDbi::keytypes", colnames(df))
      dots <- list(...)
      if ("entrez_ids" %in% names(dots)) {
        filtered_df <- df[df$ENTREZID %in% dots$entrez_ids, ]
        mockery::stub(paste_gene_symbols, "AnnotationDbi::select", filtered_df)
      }
      paste_gene_symbols(...)
    }
  }

  expect_error(
    object = build_test()(),
    info = "No input to `paste_gene_symbols` - should fail"
  )
  expect_error(
    object = build_test()(entrez_ids = NULL),
    info = paste(
      "entrez_ids should be defined in `paste_gene_symbols` - received NULL"
    )
  )
  expect_error(
    object = build_test()(entrez_ids = character(0)),
    info = "entrez_ids should be non-empty in input to paste_gene_symbols"
  )
  expect_error(
    object = build_test()(
      entrez_ids = "1234"
    ),
    info = "entrezgene_db should be defined in `paste_gene_symbols`"
  )
  expect_error(
    object = build_test(valid_db = FALSE)(
      entrez_ids = "1234", entrezgene_db = "NOT.A.DATABASE"
    ),
    info = "The DB name / DB passed to paste_gene_symbols should be an OrgDB"
  )
  expect_equal(
    object = build_test()(
      entrez_ids = "1234", entrezgene_db = "REAL.DATABASE"
    ),
    expected = c("1234|CCR5"),
    info = "Single, correctly mapping gene id"
  )
  expect_equal(
    object = build_test()(
      entrez_ids = c("1234", "10000"), entrezgene_db = "REAL.DATABASE"
    ),
    expected = c("1234|CCR5", "10000|AKT3"),
    info = "Two, correctly mapping gene ids"
  )
  expect_equal(
    object = build_test()(
      entrez_ids = c("1111"), entrezgene_db = "REAL.DATABASE"
    ),
    expected = c("1111|SYMBOL1"),
    info = paste(
      "If several symbols are available for a gene ID, only return the first"
    )
  )
  expect_error(
    object = build_test(
      df = data.frame(
        not.entrez = letters, not.symbol = LETTERS, stringsAsFactors = FALSE
      )
    )(
      entrez_ids = letters[1:3], entrezgene_db = "REAL.DATABASE"
    ),
    regexp = "ENTREZID and SYMBOL should be keytypes in `entrezgene_db`",
    info = "The database should have colnames ENTREZID and SYMBOL"
  )
  # Note that 4663 used to be gene symbol "NA" = neuroacanthocytosis
  #   although this is now "7504|XK"
  expect_equal(
    object = build_test()(
      entrez_ids = c("4", "4663", "7504"), entrezgene_db = "REAL.DATABASE"
    ),
    expected = c("4|---", "4663|---", "7504|XK"),
    info = "Non-mapping gene ids - should return 'gene.id|---'"
  )

  expect_equal(
    object = build_test()(
      entrez_ids = "956530", entrezgene_db = "REAL.DATABASE"
    ),
    expected = "956530|NA",
    info = "It handles the string 'NA'"
  )
})

###############################################################################

test_that("Testing multisymbol mapper: Symbols to Entrez", {
  # Tests using HGNC symbols as input
  lambda_fn <- function(x) {
    multisymbol_to_entrez_ids(
      gene_symbols = x,
      symbol_type = "SYMBOL",
      entrezgene_db = org.Hs.eg.db::org.Hs.eg.db,
      optim_flag = optim.flag
    )
  }

  expect_equal(
    object = lambda_fn(c()),
    expected = c(),
    info = "Multisymbol: Empty input"
  )

  expect_equal(
    object = lambda_fn(NULL),
    expected = c(),
    info = "Multisymbol: NULL input"
  )

  # expect_equal(
  #  object = lambda_fn(""),
  #  expected = as.character(NA),
  #  info = "Multisymbol: Empty string"
  # )

  #   expect_equal(
  #     object = lambda_fn("DOESNOTMAP"),
  #     expected = as.character(NA),
  #     info = "Multisymbol: Single non-mapping input"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn("AKT3"),
  #     expected = "10000",
  #     info = "Multisymbol: Singly-mapping input"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn(c("AKT3", "CCR5", "SLC38A2")),
  #     expected = c("10000", "1234", "54407"),
  #     info = "Multisymbol: Multiple singly-mapping inputs"
  #   )
  #   expect_equal(
  #     object = lambda_fn(c("AKT3|AKT3")),
  #     expected = "10000",
  #     info = "A single many:one mapping input"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn("AKT3|CCR5"),
  #     expected = "10000|1234",
  #     info = "Multisymbol: A single many:many input"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn("AKT3|CCR5"),
  #     expected = lambda_fn("CCR5|AKT3"),
  #     info = "Lexicographic sorting of the output"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn("AKT3||DOESNOTMAP"),
  #     expected = "10000",
  #     info = "Delimited string of multiple ids; only one valid map"
  #   )
  #
  #   expect_equal(
  #     object = lambda_fn(c("AKT3|XK", "CCR5|SLC38A2")),
  #     expected = c("10000|7504", "1234|54407"),
  #     info = "Vector of multiply-mapping inputs"
  #   )
})
