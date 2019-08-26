###############################################################################
# Converted to testthat format from RUnit on 14/9/2016
#
#
###############################################################################

context("Tests for rh_microarray_functions.R")

###############################################################################

library("Biobase")
# library("org.Hs.eg.db")
# library("org.Mm.eg.db")

###############################################################################
# Commonly used test input

eset_empty <- new("ExpressionSet")

eset_NA <- new(
  "ExpressionSet",
  exprs = matrix(NA)
)

eset_no_genbank <- eset_empty
Biobase::featureData(eset_no_genbank) <- Biobase::AnnotatedDataFrame(
  data = data.frame(
    NOT.A.GENBANK.COLUMN = character(0),
    stringsAsFactors = FALSE
  )
)

eset_refseq <- eset_empty
Biobase::featureData(eset_refseq) <- Biobase::AnnotatedDataFrame(
  data = data.frame(
    "RefSeq Transcript ID" = character(0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)


###############################################################################

test_that("check_if_pretransformed_eset: validity of inputs", {

  # Input must be an EpxressionSet;
  # NULL and list should fail in particular
  expect_error(
    check_if_pretransformed_eset(eset = NULL),
    info = "NULL input should fail"
  )

  expect_error(
    check_if_pretransformed_eset(eset = list()),
    info = "empty list input should fail"
  )

  expect_error(
    check_if_pretransformed_eset(
      eset = list(eset_empty)
    ),
    info = "list of ESets should fail"
  )
})

###############################################################################

test_that("check_if_pretransformed_eset: correct outputs", {

  # Empty ESets: Assume it has been log-transformed
  expect_equal(
    check_if_pretransformed_eset(eset = eset_empty),
    TRUE,
    info = "Empty eset is transformed"
  )

  # An ESet containing all NA entries: Assume it has been log-transformed
  expect_equal(
    check_if_pretransformed_eset(eset = eset_NA),
    TRUE,
    info = "All-NA eset is transformed"
  )

  # An ESet that contains some negative values: Assume it has been
  # log-trasnformed
  expect_equal(
    check_if_pretransformed_eset(
      eset = new("ExpressionSet", exprs = matrix(-1))
    ),
    TRUE,
    info = "Negative-containing eset is transformed"
  )

  # An ESet with a high range between highest and lowest values:
  #   - Assume the eset is not log-transformed
  eset_high_range <- new(
    "ExpressionSet",
    exprs = matrix(c(1, 1002), nrow = 2)
  )
  expect_equal(
    check_if_pretransformed_eset(
      eset = eset_high_range,
      range_limit_if_transformed = 1000 # default range limit
    ),
    FALSE,
    info = "High-range datasets are assumed untransformed"
  )

  # An Eset with a high range and containing some NA values
  # - Assume the eset is not log-transformed
  eset_high_range_NA <- new(
    "ExpressionSet",
    exprs = matrix(
      # the function should determine the range as 1001
      # note - if na.omit(matrix) was used, 1st row would be dropped and
      #        the range would be 1, but na.inf.omit() converts to vector
      c(
        1, NA,
        1001, 1002
      ),
      nrow = 2, byrow = TRUE
    )
  )
  expect_equal(
    check_if_pretransformed_eset(
      eset = eset_high_range_NA,
      range_limit_if_transformed = 1000
    ),
    FALSE,
    info = paste(
      "High-range (=> untransformed) nature of the dataset depends",
      "on an NA-containing row"
    )
  )

  # Highly-skewed datasets: Assume the eset is not log-transformed
  # Default skewness threshold is 1, and test passes if all columns are skewed
  eset_skewed <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        2, 3,
        3, 4,
        10, 50
      ),
      nrow = 4, byrow = TRUE
    )
  )
  expect_equal(
    check_if_pretransformed_eset(
      eset = eset_skewed,
      skew_threshold = 1 # default skew threshold
    ),
    FALSE,
    info = "Highly-skewed low-range non-negative dataset: assume untransformed"
  )

  # IGNORE
  # If the IQR-normalised difference between mean and median is above a
  # threshold
  # ? How to test - would all such datasets also be 'highly skewed' ? #
})

###############################################################################

test_that("has_refseq_column", {
  # Empty ESet has no annotations, and hence no RefSeq column:
  expect_equal(
    has_refseq_column(eset_empty),
    FALSE,
    info = "an empty `Eset` should have no RefSeq column"
  )

  # ESet that has no GenBank / RefSeq column:
  expect_equal(
    has_refseq_column(eset_no_genbank),
    FALSE,
    info = "an `Eset` that has no RefSeq column"
  )

  # ESet that has a Refseq column
  expect_equal(
    has_refseq_column(eset_refseq),
    TRUE,
    info = "an `Eset` that has a Refseq column"
  )
})

###############################################################################

test_that("get_refseq_column: validity of inputs", {

  # Input must be an ExpressionSet
  expect_error(
    get_refseq_column(),
    info = "empty call to get_refseq_column fails"
  )

  expect_error(
    get_refseq_column(NULL),
    info = "NULL call to get_refseq_column fails"
  )

  expect_error(
    get_refseq_column(character(0)),
    info = "non-ESet input to get_refseq_column fails"
  )

  # Input should contain a genbank / refseq column
  #   and the featureData entry of the ESet should have a varLabels entry
  expect_error(
    get_refseq_column(eset = eset_empty),
    info = "the expressionSet has some valid valLabels in it's featureData"
  )

  # The genbank / refseq column should be from the set {RefSeq Transcript ID,
  #   GB_LIST, GB_ACC, RefSeq.Transcript.ID} and should be a character or
  #   factor column
  expect_error(
    get_refseq_column(eset = eset_no_genbank),
    info = paste(
      "a valid genbank column ('RefSeq Transcript ID', 'GB_LIST',",
      "'GB_ACC') should be present"
    )
  )

  eset_nonCharacter_genbank <- eset_empty
  featureData(eset_nonCharacter_genbank) <- AnnotatedDataFrame(
    data = data.frame("GB_LIST" = logical(0))
  )
  expect_error(
    get_refseq_column(eset = eset_nonCharacter_genbank),
    info = "the genbank column is a character of factor"
  )
})

###############################################################################

test_that("get_refseq_column: correct outputs", {

  # A genbank column should be returned in order of preference:
  #   'RefSeq Transcript ID' > 'GB_LIST' > 'GB_ACC'
  # The genbank column should either be a factor or contain characters
  expect_equal(
    get_refseq_column(eset_refseq),
    "RefSeq Transcript ID",
    info = "RefSeq column is found"
  )

  # Genbank column pulled out even if it's a factor:
  eset_refseq_factor <- eset_empty
  featureData(eset_refseq_factor) <- AnnotatedDataFrame(
    data = data.frame(
      "RefSeq Transcript ID" = character(0),
      check.names = FALSE,
      stringsAsFactors = TRUE
    )
  )
  expect_equal(
    get_refseq_column(eset_refseq_factor),
    "RefSeq Transcript ID",
    info = "RefSeq column is found, even if it is a factor"
  )

  # The function should be robust to application of R's make.names function to
  #   the featureData columns
  eset_refseq_makenames <- eset_empty
  featureData(eset_refseq) <- AnnotatedDataFrame(
    data = data.frame(
      "RefSeq.Transcript.ID" = character(0),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(
    get_refseq_column(eset_refseq),
    "RefSeq.Transcript.ID",
    info = paste(
      '"RefSeq Transcript ID" column was identified after',
      "make.names conversion to RefSeq.Transcript.ID"
    )
  )

  # The function should pick GB_LIST before GB_ACC
  eset_gbacc_first <- eset_empty
  featureData(eset_gbacc_first) <- AnnotatedDataFrame(
    data = data.frame(
      "GB_ACC" = character(0),
      "GB_LIST" = character(0),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(
    get_refseq_column(eset_gbacc_first),
    "GB_LIST",
    info = "GB_LIST should be picked before GB_ACC"
  )
})

###############################################################################

test_that("add_entrez_ids_to_esets: validity of inputs", {

  # - Should work using both an ExpressionSet and a list of ExpressionSets as
  # input
  # - All ExpressionSets should have a defined 'platform' and a defined
  # 'featureData' slot (the latter containing at least one of RefSeq Transcript
  # ID, GB_LIST or GB_ACC)
  # - The list of ExpressionSets need not all have the same 'platform'

  # If the input is NULL, the function should fail
  expect_error(
    add_entrez_ids_to_esets(),
    info = "add_entrez_ids_to_esets fails if no input is provided"
  )

  expect_error(
    add_entrez_ids_to_esets(esets = NULL),
    info = "add_entrez_ids_to_esets fails on NULL input"
  )

  expect_error(
    add_entrez_ids_to_esets(esets = list(NULL)),
    info = "add_entrez_ids_to_esets fails on List-of-NULL inputs"
  )

  # If any of the inputs are not ExpressionSets, the function should fail
  expect_error(
    add_entrez_ids_to_esets(esets = "NOT AN ExpressionSet"),
    info = "add_entrez_ids_to_esets fails if input is not an ExpressionSet"
  )

  # If the annotation slot is empty, the function should fail
  # - note that Biobase::ExpressionSet can only take strings in the
  # annotation(.) slot
  #
  eset_nullPlatform <- eset_empty
  annotation(eset_nullPlatform) <- character(0)
  expect_error(
    add_entrez_ids_to_esets(esets = eset_nullPlatform),
    info = paste(
      "add_entrez_ids_to_esets fails if any ESet has an empty",
      "'annotation'"
    )
  )

  eset_emptyPlatform <- eset_empty
  annotation(eset_emptyPlatform) <- ""
  expect_error(
    add_entrez_ids_to_esets(esets = eset_emptyPlatform),
    info = paste(
      "add_entrez_ids_to_esets fails if any ESet has an empty",
      "'annotation'"
    )
  )

  # If the featureData does not have a defined refseq column, the function
  # should fail
  # - implicitly tested via get_refseq_column
  eset_noRefseqCol <- eset_empty
  annotation(eset_noRefseqCol) <- "SOME PLATFORM"
  featureData(eset_noRefseqCol) <- AnnotatedDataFrame(
    data.frame(
      "NOT A REFSEQ COLUMN" = character(0),
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    add_entrez_ids_to_esets(esets = eset_noRefseqCol),
    info = paste(
      "add_entrez_ids_to_esets fails unless all ESets have a valid ",
      "refseq/genbank column"
    )
  )
})

###############################################################################

# test_that("add_entrez_ids_to_esets: correct outputs", {
#   result_helper <- function(esets, db = org.Hs.eg.db::org.Hs.eg.db) {
#     res <- add_entrez_ids_to_esets(
#       esets = esets,
#       entrezgene.db = db
#     )
#     featureData(res)[["entrez.id"]]
#   }

#   # Logic tests (low level stuff is done by multisymbol_to_entrez_ids):
#   # Blank refseq entries, NA refseq entries should map to NA
#   eset_blankRefseqs <- eset_empty
#   annotation(eset_blankRefseqs) <- "SOME_OTHER_PLATFORM"
#   featureData(eset_blankRefseqs) <- AnnotatedDataFrame(
#     data.frame("GB_ACC" = c(NA, ""), stringsAsFactors = FALSE)
#   )
#   expect_blankRefseqs <- rep(as.character(NA), 2)
#   result_blankRefseqs <- result_helper(eset_blankRefseqs)
#   expect_equal(
#     expect_blankRefseqs,
#     result_blankRefseqs,
#     info = "NA and empty string should give NA entrez.ids"
#   )

#   # A single refseq entry that maps to a single human gene
#   eset_singleRefseq <- eset_empty
#   annotation(eset_singleRefseq) <- "platform.1"
#   featureData(eset_singleRefseq) <- AnnotatedDataFrame(
#     data.frame(
#       "GB_LIST" = "NM_000579", # CCR5/'1234'
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_singleRefseq <- "1234"
#   result_singleRefseq <- result_helper(eset_singleRefseq)
#   expect_equal(
#     expect_singleRefseq,
#     result_singleRefseq,
#     info = "Single refseq id that maps to a single entrez id"
#   )

#   # A ' /// '-separated entry of refseq ids
#   eset_twoRefseq <- eset_empty
#   annotation(eset_twoRefseq) <- "GPL12345"
#   featureData(eset_twoRefseq) <- AnnotatedDataFrame(
#     data.frame(
#       "GB_LIST" = "NM_001307936 /// NM_018976",
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_twoRefseq <- "54407"
#   result_twoRefseq <- result_helper(eset_twoRefseq)
#   expect_equal(
#     expect_twoRefseq,
#     result_twoRefseq,
#     info = "Two refseqs, ///-separated, that map to a single entrez id"
#   )

#   # A ','-separated entry of refseq ids
#   eset_commaRefseq <- eset_empty
#   annotation(eset_commaRefseq) <- "GPL9876"
#   featureData(eset_commaRefseq) <- AnnotatedDataFrame(
#     data.frame(
#       "GB_LIST" = "NM_001130045,NM_153254,BC126152",
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_commaRefseq <- "254173"
#   result_commaRefseq <- result_helper(eset_commaRefseq)
#   expect_equal(
#     expect_commaRefseq,
#     result_commaRefseq,
#     info = "Comma-separated refseq list, that map to a single entrez id"
#   )

#   # refseq to entrez mapping for two different mouse probes:
#   eset_mouseRefseq <- eset_empty
#   annotation(eset_mouseRefseq) <- "GPL1261"
#   featureData(eset_mouseRefseq) <- AnnotatedDataFrame(
#     data.frame(
#       "RefSeq Transcript ID" = c("NM_017477 /// NM_201244", "NM_013477"),
#       stringsAsFactors = TRUE,
#       check.names = FALSE
#     )
#   )
#   expect_mouseRefseq <- c("54161", "11972")
#   result_mouseRefseq <- result_helper(eset_mouseRefseq,
#     org.Mm.eg.db::org.Mm.eg.db)
#   expect_equal(
#     expect_mouseRefseq,
#     result_mouseRefseq,
#     info = "refseq to entrez mapping for two mouse probes"
#   )

#   # Add entrez.ids to a swissprot-containing ESet
#   eset_swissprot <- eset_empty
#   annotation(eset_swissprot) <- "SOME GPL ID"
#   featureData(eset_swissprot) <- AnnotatedDataFrame(
#     data.frame(
#       "swissprot" = c(NA, "", "---", "ENST0000412115", "NR_046018",
#       "NM_001005221"),
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_swissprot <- c(NA, NA, NA, NA, "100287102", "729759")
#   result_swissprot <- result_helper(eset_swissprot)
#   expect_equal(
#     expect_swissprot,
#     result_swissprot,
#     info = "swissprot column can be mapped to Entrez ids"
#   )
# })

###############################################################################

test_that("swissprot_column_to_refseq", {

  # no input
  # "---" as input
  # /// and // splits - both mapping and nonmapping versions
  # factor input

  expect_error(
    swissprot_column_to_refseq(),
    info = "No input to swissprot_column_to_refseq should fail"
  )

  expect_error(
    swissprot_column_to_refseq(NULL),
    info = "NULL input to swissprot_column_to_refseq should fail"
  )

  expect_error(
    swissprot_column_to_refseq(character(0)),
    info = "Empty vector input to swissprot_column_to_refseq should fail"
  )

  expect_equal(
    swissprot_column_to_refseq(as.character(NA)),
    "",
    info = "swissprot_column_to_refseq maps NA string to empty string"
  )

  invalid_vec <- c(NA, "", "---", "ENST00000412115")
  invalid_expect <- rep("", 4)
  expect_equal(
    swissprot_column_to_refseq(invalid_vec),
    invalid_expect,
    info = paste(
      "swissprot_column_to_refseq maps ---, empty string and",
      "refseq-absent to empty string"
    )
  )

  expect_equal(
    swissprot_column_to_refseq(
      c("NR_046018", "NM_001005221", "ENST00000412115")
    ),
    c("NR_046018", "NM_001005221", ""),
    info = "refseq ids map to refseq ids and non-refseq ids do not"
  )

  expect_equal(
    swissprot_column_to_refseq(c("NR_046018.2")),
    "NR_046018",
    info = paste(
      "swissprot_column_to_refseq drops NR_XXXXX.1/2/3... to",
      "generate NR_XXXX"
    )
  )

  valid_vec <- c(
    "NR_028325 // B4DYM5",
    "NR_028325 // B4DYM5 /// NR_028325 // B4E0H4",
    "NM_001005221 // Q6IEY1 /// BC137547 // Q6IEY1",
    paste(
      "NM_001100114 // Q5VT98 /// NM_001099852 // Q5VT98 ///",
      "ENST00000327795 // Q5VT98"
    )
  )
  valid_expect <- c(
    "NR_028325",
    "NR_028325",
    "NM_001005221",
    "NM_001100114,NM_001099852"
  )
  expect_equal(
    swissprot_column_to_refseq(valid_vec),
    valid_expect,
    info = "Correct parsing of ///-separated and //-subseparated entries"
  )

  expect_equal(
    swissprot_column_to_refseq(factor(c(valid_vec, invalid_vec))),
    c(valid_expect, invalid_expect),
    info = "Correct parsing of factors by swissprot_column_to_refseq"
  )
})

###############################################################################

test_that("Unit tests for filter_and_transform_eset", {
  ft_runner <- function(
                          # As for filter_and_transform_eset, but with all
                          # booleans set to FALSE
                          eset,
                          log2_transform = FALSE,
                          drop_row_if_duplicated = FALSE,
                          drop_row_if_zero_variance = FALSE,
                          convert_inf_to_na = FALSE,
                          normalise_method = "none",
                          drop_row_na_inf_threshold = 0.25) {
    filter_and_transform_eset(
      eset = eset,
      log2_transform = log2_transform,
      drop_row_if_duplicated = drop_row_if_duplicated,
      drop_row_if_zero_variance = drop_row_if_zero_variance,
      convert_inf_to_na = convert_inf_to_na,
      normalise_method = normalise_method,
      drop_row_na_inf_threshold = drop_row_na_inf_threshold
    )
  }

  # Input must be an ExpressionSet
  expect_error(
    filter_and_transform_eset(),
    info = "No input to filter_and_transform_eset"
  )
  expect_error(
    filter_and_transform_eset(eset = NULL),
    info = "NULL input to filter_and_transform_eset"
  )
  expect_error(
    filter_and_transform_eset(eset = list()),
    info = "Non-Eset input to filter_and_transform_eset should fail"
  )
  expect_error(
    filter_and_transform_eset(eset = list(new("ExpressionSet"))),
    info = "list of ESets to filter_and_transform_eset should fail"
  )

  # Transformation of an empty eset should have no effect
  expect_equal(
    filter_and_transform_eset(
      eset = eset_empty
    ),
    eset_empty,
    info = "transformation of an empty eset should have no effect"
  )

  # Log2 transform an eset - all positive
  dnames_allPos <- list(1:2, 1:3)
  eset_allPos <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2, 4,
        1, 1 / 2, 1 / 4
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_allPos
    )
  )
  expect_allPos <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, 1, 2,
        0, -1, -2
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_allPos
    )
  )
  expect_equal(
    ft_runner(
      eset = eset_allPos,
      log2_transform = TRUE,
    ),
    expect_allPos,
    info = "Log2 transformation of an all-positive ESet"
  )

  # Log2 transform an eset - some NA present
  dnames_someNAs <- list(1:2, 1:3)
  eset_someNAs <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        NA, 2, 4,
        1, 1 / 2, NaN
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_someNAs
    )
  )
  expect_someNAs <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        NA, 1, 2,
        0, -1, NA
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_someNAs
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_someNAs,
      log2_transform = TRUE,
      drop_row_na_inf_threshold = 1,
    ),
    expect_someNAs,
    info = "transform Eset, some input entries are NA"
  )

  # Fail to Log2 transform an eset containing negatives
  dnames_someNeg <- list(1:2, 1:3)
  eset_someNeg <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        -1, 2, NA,
        1, -1 / 2, 1 / 4
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_someNeg
    )
  )
  expect_error(
    ft_runner(
      eset = eset_someNeg,
      log2_transform = TRUE
    ),
    info = "fail to log2-transform an eset containing negatives"
  )

  # Drop a duplicated row
  dnames_dupRow <- list(1:3, 1:2)
  eset_dupRow <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, 4,
        1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_dupRow
    )
  )

  expect_dupRow_dropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, 4
      ),
      nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(1:2, 1:2)
    )
  )

  expect_dupRow_notDropped <- eset_dupRow

  expect_equal(
    ft_runner(
      eset_dupRow,
      drop_row_if_duplicated = TRUE
    ),
    expect_dupRow_dropped,
    info = "Dropping a duplicated row"
  )

  expect_equal(
    ft_runner(
      eset_dupRow
    ),
    expect_dupRow_notDropped,
    info = "Not dropping a duplicated row"
  )

  # Drop a duplicated row - NA present
  dnames_dupRow_NA <- list(1:3, 1:2)
  eset_dupRow_NA <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, NA,
        3, NA
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_dupRow_NA
    )
  )
  expect_dupRow_NA_dropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, NA
      ),
      nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(1:2, 1:2)
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_dupRow_NA,
      drop_row_if_duplicated = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_dupRow_NA_dropped,
    info = "Dropping a duplicated, NA-containing row"
  )

  # Drop a duplicated row - Inf present
  dnames_dupRow_Inf <- list(1:3, 1:2)
  eset_dupRow_Inf <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf,
        1, Inf,
        3, 4
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_dupRow_Inf
    )
  )
  expect_dupRow_Int_dropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf,
        3, 4
      ),
      nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(c(1, 3), 1:2)
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_dupRow_Inf,
      drop_row_if_duplicated = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_dupRow_Int_dropped,
    info = "Dropping a duplicated, Inf-containing row"
  )

  # Convert Infs to NA during log-transformation
  dnames_Inf <- list(1:3, 1:2)
  eset_Inf <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf, # Inf --> Inf
        0, 2, # 0   --> -Inf
        1 / 2, 4
      ),
      nrow = 3, ncol = 2, byrow = TRUE,
      dimnames = dnames_Inf
    )
  )

  expect_Inf_converted <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, NA,
        NA, 1,
        -1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_Inf
    )
  )

  expect_Inf_notConverted <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, Inf,
        -Inf, 1,
        -1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_Inf
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_Inf,
      log2_transform = TRUE,
      convert_inf_to_na = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_Inf_converted,
    info = "Converting Infs to NA during log-transformation"
  )

  expect_equal(
    ft_runner(
      eset = eset_Inf,
      log2_transform = TRUE,
      convert_inf_to_na = FALSE,
      drop_row_na_inf_threshold = 1
    ),
    expect_Inf_notConverted,
    info = "Not converting Infs to NA during log-transformation"
  )

  # Drop rows with too many NA/Inf entries
  dnames_na50 <- list(1:4, 1:3)
  eset_na50 <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf, NA, # > 50% are NA or Inf
        0, 2, NA, # < 50% are NA or Inf
        1 / 2, 4, 8, # < 50% are NA or Inf
        NA, NA, NA # 100% are NA
      ),
      nrow = 4, ncol = 3, byrow = TRUE, dimnames = dnames_na50
    )
  )

  expect_na50_dropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, 2, NA,
        1 / 2, 4, 8
      ),
      nrow = 2, ncol = 3, byrow = TRUE, dimnames = list(c(2, 3), 1:3)
    )
  )

  expect_na50_notDropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf, NA,
        0, 2, NA,
        1 / 2, 4, 8,
        NA, NA, NA # ensure this line is not dropped
      ),
      nrow = 4, ncol = 3, byrow = TRUE, dimnames = dnames_na50
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_na50,
      drop_row_na_inf_threshold = 0.5
    ),
    expect_na50_dropped,
    info = "Dropping rows with <= 50% NA/Inf"
  )

  expect_equal(
    ft_runner(
      eset = eset_na50,
      drop_row_na_inf_threshold = 1
    ),
    expect_na50_notDropped,
    info = "Keeping all rows, despite NA/Inf presence"
  )

  # Test median-subtraction and IQR normalisation step:
  # - This step acts on the rows of the expression data
  #   and should act before any Infs are converted to NA
  dnames_med_iqr <- list(1:5, 1:4)
  eset_med_iqr <- new(
    "ExpressionSet",
    # col1 - all entries are different; med = 3, IQR = 4-2 = 2
    # col2 - contains some NAs; med = 2, IQR = 2.5-1.5 = 1
    # col3 - all entries are the same; med = 1, IQR=0; output should be NaN
    # col4 - includes an Inf that should be ignored; med = 2.5, IQR = 3.75 -
    # 1.75 = 2
    exprs = matrix(
      c(
        1, 1, 1, 1,
        2, 2, 1, 2,
        3, 3, 1, 3,
        4, NA, 1, Inf,
        5, NA, 1, 6
      ),
      nrow = 5, ncol = 4, byrow = TRUE, dimnames = dnames_med_iqr
    )
  )

  # Expected values on subtracting column medians
  expect_exprs_med <- matrix(
    c(
      -2, -1, 0, -1.5,
      -1, 0, 0, -0.5,
      0, 1, 0, 0.5,
      1, NA, 0, Inf,
      2, NA, 0, 3.5
    ),
    nrow = 5, ncol = 4, byrow = TRUE, dimnames = dnames_med_iqr
  )

  expect_equal(
    ft_runner(
      eset_med_iqr,
      drop_row_na_inf_threshold = 1,
      normalise_method = "median"
    ),
    new("ExpressionSet", exprs = expect_exprs_med),
    info = "Median normalisation"
  )

  # Expected values on subtracting medians and then dividing by IQRs
  expect_exprs_med_iqr <- matrix(
    c(
      -1, -1, NaN, -0.75,
      -0.5, 0, NaN, -0.25,
      0, 1, NaN, 0.25,
      0.5, NA, NaN, Inf,
      1, NA, NaN, 1.75
    ),
    nrow = 5, ncol = 4, byrow = TRUE, dimnames = dnames_med_iqr
  )

  expect_equal(
    ft_runner(
      eset_med_iqr,
      drop_row_na_inf_threshold = 1,
      normalise_method = "median_iqr"
    ),
    new("ExpressionSet", exprs = expect_exprs_med_iqr),
    info = "Median / IQR normalisation"
  )

  # Expected values on subtracting medians, dividng by IQRs and then
  #   transforming Infs to NA
  expect_exprs_med_iqr_inf2NA <- matrix(
    c(
      -1, -1, NaN, -0.75,
      -0.5, 0, NaN, -0.25,
      0, 1, NaN, 0.25,
      0.5, NA, NaN, NA,
      1, NA, NaN, 1.75
    ),
    nrow = 5, ncol = 4, byrow = TRUE, dimnames = dnames_med_iqr
  )

  expect_equal(
    ft_runner(
      eset_med_iqr,
      convert_inf_to_na = TRUE,
      drop_row_na_inf_threshold = 1,
      normalise_method = "median_iqr"
    ),
    new("ExpressionSet", exprs = expect_exprs_med_iqr_inf2NA),
    info = "Median / IQR normalisation and Inf --> NA"
  )

  # Tests for drop_row_if_zero_varianceiance
  dnames_zeroVar <- list(1:2, 1:4)
  eset_zeroVar <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 1, 1, 1,
        1, NA, 1, 1
      ),
      nrow = 2, ncol = 4, byrow = TRUE, dimnames = dnames_zeroVar
    )
  )
  expect_zeroVar <- eset_zeroVar[c(), ]

  expect_equal(
    ft_runner(
      eset_zeroVar,
      drop_row_na_inf_threshold = 1,
      drop_row_if_zero_variance = TRUE
    ),
    expect_zeroVar,
    info = "Dropping rows that have no inter-sample variance"
  )
})

###############################################################################
