###############################################################################

context("Tests for rh_microarray_functions.R")

###############################################################################

library("Biobase")

adf <- Biobase::AnnotatedDataFrame

# library("org.Hs.eg.db")
# library("org.Mm.eg.db")

###############################################################################
# Commonly used test input

eset_empty <- new("ExpressionSet")

eset_with_only_na_values <- new(
  "ExpressionSet",
  exprs = matrix(NA)
)

eset_no_genbank <- eset_empty
Biobase::featureData(eset_no_genbank) <- adf(
  data = data.frame(
    NOT.A.GENBANK.COLUMN = character(0),
    stringsAsFactors = FALSE
  )
)

eset_refseq <- eset_empty
Biobase::featureData(eset_refseq) <- adf(
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
    check_if_pretransformed_eset(eset = eset_with_only_na_values),
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
  eset_high_range_with_na_values <- new(
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
      eset = eset_high_range_with_na_values,
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

  eset_non_character_genbank <- eset_empty

  Biobase::featureData(
    eset_non_character_genbank
  ) <- adf(
    data = data.frame("GB_LIST" = logical(0))
  )
  expect_error(
    get_refseq_column(eset = eset_non_character_genbank),
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
  Biobase::featureData(eset_refseq_factor) <- adf(
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
  Biobase::featureData(eset_refseq) <- adf(
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
  Biobase::featureData(eset_gbacc_first) <- adf(
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
  eset_null_platform <- eset_empty
  annotation(eset_null_platform) <- character(0)
  expect_error(
    add_entrez_ids_to_esets(esets = eset_null_platform),
    info = paste(
      "add_entrez_ids_to_esets fails if any ESet has an empty",
      "'annotation'"
    )
  )

  eset_empty_platform <- eset_empty
  annotation(eset_empty_platform) <- ""
  expect_error(
    add_entrez_ids_to_esets(esets = eset_empty_platform),
    info = paste(
      "add_entrez_ids_to_esets fails if any ESet has an empty",
      "'annotation'"
    )
  )

  # If the featureData does not have a defined refseq column, the function
  # should fail
  # - implicitly tested via get_refseq_column
  eset_no_refseq_column <- eset_empty
  annotation(eset_no_refseq_column) <- "SOME PLATFORM"
  Biobase::featureData(eset_no_refseq_column) <- adf(
    data.frame(
      "NOT A REFSEQ COLUMN" = character(0),
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    add_entrez_ids_to_esets(esets = eset_no_refseq_column),
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
#     Biobase::featureData(res)[["entrez.id"]]
#   }

#   # Logic tests (low level stuff is done by multisymbol_to_entrez_ids):
#   # Blank refseq entries, NA refseq entries should map to NA
#   eset_blank_refseqs <- eset_empty
#   annotation(eset_blank_refseqs) <- "SOME_OTHER_PLATFORM"
#   Biobase::featureData(eset_blank_refseqs) <- adf(
#     data.frame("GB_ACC" = c(NA, ""), stringsAsFactors = FALSE)
#   )
#   expect_blank_refseqs <- rep(as.character(NA), 2)
#   result_blank_refseqs <- result_helper(eset_blank_refseqs)
#   expect_equal(
#     expect_blank_refseqs,
#     result_blank_refseqs,
#     info = "NA and empty string should give NA entrez.ids"
#   )

#   # A single refseq entry that maps to a single human gene
#   eset_single_refseq <- eset_empty
#   annotation(eset_single_refseq) <- "platform.1"
#   Biobase::featureData(eset_single_refseq) <- adf(
#     data.frame(
#       "GB_LIST" = "NM_000579", # CCR5/'1234'
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_single_refseq <- "1234"
#   result_single_refseq <- result_helper(eset_single_refseq)
#   expect_equal(
#     expect_single_refseq,
#     result_single_refseq,
#     info = "Single refseq id that maps to a single entrez id"
#   )

#   # A ' /// '-separated entry of refseq ids
#   eset_two_refseq <- eset_empty
#   annotation(eset_two_refseq) <- "GPL12345"
#   Biobase::featureData(eset_two_refseq) <- adf(
#     data.frame(
#       "GB_LIST" = "NM_001307936 /// NM_018976",
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_two_refseq <- "54407"
#   result_two_refseq <- result_helper(eset_two_refseq)
#   expect_equal(
#     expect_two_refseq,
#     result_two_refseq,
#     info = "Two refseqs, ///-separated, that map to a single entrez id"
#   )

#   # A ','-separated entry of refseq ids
#   eset_comma_refseq <- eset_empty
#   annotation(eset_comma_refseq) <- "GPL9876"
#   Biobase::featureData(eset_comma_refseq) <- adf(
#     data.frame(
#       "GB_LIST" = "NM_001130045,NM_153254,BC126152",
#       stringsAsFactors = FALSE
#     )
#   )
#   expect_comma_refseq <- "254173"
#   result_comma_refseq <- result_helper(eset_comma_refseq)
#   expect_equal(
#     expect_comma_refseq,
#     result_comma_refseq,
#     info = "Comma-separated refseq list, that map to a single entrez id"
#   )

#   # refseq to entrez mapping for two different mouse probes:
#   eset_mouse_refseq <- eset_empty
#   annotation(eset_mouse_refseq) <- "GPL1261"
#   Biobase::featureData(eset_mouse_refseq) <- adf(
#     data.frame(
#       "RefSeq Transcript ID" = c("NM_017477 /// NM_201244", "NM_013477"),
#       stringsAsFactors = TRUE,
#       check.names = FALSE
#     )
#   )
#   expect_mouse_refseq <- c("54161", "11972")
#   result_mouse_refseq <- result_helper(eset_mouse_refseq,
#     org.Mm.eg.db::org.Mm.eg.db)
#   expect_equal(
#     expect_mouse_refseq,
#     result_mouse_refseq,
#     info = "refseq to entrez mapping for two mouse probes"
#   )

#   # Add entrez.ids to a swissprot-containing ESet
#   eset_swissprot <- eset_empty
#   annotation(eset_swissprot) <- "SOME GPL ID"
#   Biobase::featureData(eset_swissprot) <- adf(
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
  dnames_all_positive <- list(1:2, 1:3)
  eset_all_positive <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2, 4,
        1, 1 / 2, 1 / 4
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_all_positive
    )
  )
  expect_all_positive <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, 1, 2,
        0, -1, -2
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_all_positive
    )
  )
  expect_equal(
    ft_runner(
      eset = eset_all_positive,
      log2_transform = TRUE
    ),
    expect_all_positive,
    info = "Log2 transformation of an all-positive ESet"
  )

  # Log2 transform an eset - some NA present
  dnames_some_na_values <- list(1:2, 1:3)
  eset_some_na_values <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        NA, 2, 4,
        1, 1 / 2, NaN
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_some_na_values
    )
  )
  expect_some_na_values <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        NA, 1, 2,
        0, -1, NA
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_some_na_values
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_some_na_values,
      log2_transform = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_some_na_values,
    info = "transform Eset, some input entries are NA"
  )

  # Fail to Log2 transform an eset containing negatives
  dnames_some_negative_values <- list(1:2, 1:3)
  eset_some_negative_values <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        -1, 2, NA,
        1, -1 / 2, 1 / 4
      ),
      nrow = 2, byrow = TRUE, dimnames = dnames_some_negative_values
    )
  )
  expect_error(
    ft_runner(
      eset = eset_some_negative_values,
      log2_transform = TRUE
    ),
    info = "fail to log2-transform an eset containing negatives"
  )

  # Drop a duplicated row
  dnames_dupd_row <- list(1:3, 1:2)
  eset_dupd_row <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, 4,
        1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_dupd_row
    )
  )

  expect_dupd_row_dropped <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, 4
      ),
      nrow = 2, ncol = 2, byrow = TRUE, dimnames = list(1:2, 1:2)
    )
  )

  expect_dupd_row_not_dropped <- eset_dupd_row

  expect_equal(
    ft_runner(
      eset_dupd_row,
      drop_row_if_duplicated = TRUE
    ),
    expect_dupd_row_dropped,
    info = "Dropping a duplicated row"
  )

  expect_equal(
    ft_runner(
      eset_dupd_row
    ),
    expect_dupd_row_not_dropped,
    info = "Not dropping a duplicated row"
  )

  # Drop a duplicated row - NA present
  dnames_dupd_row_with_na_values <- list(1:3, 1:2)
  eset_dupd_row_with_na_values <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 2,
        3, NA,
        3, NA
      ),
      nrow = 3, ncol = 2, byrow = TRUE,
      dimnames = dnames_dupd_row_with_na_values
    )
  )
  expect_dupd_row_with_na_values_dropped <- new(
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
      eset = eset_dupd_row_with_na_values,
      drop_row_if_duplicated = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_dupd_row_with_na_values_dropped,
    info = "Dropping a duplicated, NA-containing row"
  )

  # Drop a duplicated row - Inf present
  dnames_dupd_row_with_inf_values <- list(1:3, 1:2)
  eset_dupd_row_with_inf_values <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf,
        1, Inf,
        3, 4
      ),
      nrow = 3, ncol = 2, byrow = TRUE,
      dimnames = dnames_dupd_row_with_inf_values
    )
  )
  expect_dupd_row_inf_dropped <- new(
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
      eset = eset_dupd_row_with_inf_values,
      drop_row_if_duplicated = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_dupd_row_inf_dropped,
    info = "Dropping a duplicated, Inf-containing row"
  )

  # Convert Infs to NA during log-transformation
  dnames_inf <- list(1:3, 1:2)
  eset_inf <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, Inf, # Inf --> Inf
        0, 2, # 0   --> -Inf
        1 / 2, 4
      ),
      nrow = 3, ncol = 2, byrow = TRUE,
      dimnames = dnames_inf
    )
  )

  expect_inf_converted <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, NA,
        NA, 1,
        -1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_inf
    )
  )

  expect_inf_not_converted <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        0, Inf,
        -Inf, 1,
        -1, 2
      ),
      nrow = 3, ncol = 2, byrow = TRUE, dimnames = dnames_inf
    )
  )

  expect_equal(
    ft_runner(
      eset = eset_inf,
      log2_transform = TRUE,
      convert_inf_to_na = TRUE,
      drop_row_na_inf_threshold = 1
    ),
    expect_inf_converted,
    info = "Converting Infs to NA during log-transformation"
  )

  expect_equal(
    ft_runner(
      eset = eset_inf,
      log2_transform = TRUE,
      convert_inf_to_na = FALSE,
      drop_row_na_inf_threshold = 1
    ),
    expect_inf_not_converted,
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

  expect_na50_not_dropped <- new(
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
    expect_na50_not_dropped,
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
  expect_exprs_med_iqr_inf_to_na <- matrix(
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
    new("ExpressionSet", exprs = expect_exprs_med_iqr_inf_to_na),
    info = "Median / IQR normalisation and Inf --> NA"
  )

  # Tests for drop_row_if_zero_varianceiance
  dnames_zero_variance <- list(1:2, 1:4)
  eset_zero_variance <- new(
    "ExpressionSet",
    exprs = matrix(
      c(
        1, 1, 1, 1,
        1, NA, 1, 1
      ),
      nrow = 2, ncol = 4, byrow = TRUE, dimnames = dnames_zero_variance
    )
  )
  expect_zero_variance <- eset_zero_variance[c(), ]

  expect_equal(
    ft_runner(
      eset_zero_variance,
      drop_row_na_inf_threshold = 1,
      drop_row_if_zero_variance = TRUE
    ),
    expect_zero_variance,
    info = "Dropping rows that have no inter-sample variance"
  )
})

###############################################################################

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
  eset_null_platform <- eset_empty
  annotation(eset_null_platform) <- character(0)
  expect_error(
    add_entrez_ids_to_esets(esets = eset_null_platform),
    info = paste("add_entrez_ids_to_esets fails if any ESet has an empty",
                 "'annotation'")
  )

  eset_empty_platform <- eset_empty
  annotation(eset_empty_platform) <- ""
  expect_error(
    add_entrez_ids_to_esets(esets = eset_empty_platform),
    info = paste("add_entrez_ids_to_esets fails if any ESet has an empty",
                 "'annotation'")
  )


  # If the featureData does not have a defined refseq column, the function
  # should fail
  # - implicitly tested via get_refseq_column
  eset_no_refseq_col <- eset_empty
  annotation(eset_no_refseq_col) <- "SOME PLATFORM"
  Biobase::featureData(eset_no_refseq_col) <- adf(
    data.frame(
      "NOT A REFSEQ COLUMN" = character(0),
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    add_entrez_ids_to_esets(esets = eset_no_refseq_col),
    info = paste(
      "add_entrez_ids_to_esets fails unless all ESets have a valid ",
      "refseq/genbank column"
    )
  )
})

test_that("add_entrez_ids_to_esets: correct outputs", {
  # TODO: mock out the use of org.Hs.eg.db and org.Mm.eg.db

  result_helper <- function(esets, db = org.Hs.eg.db::org.Hs.eg.db) {
    res <- add_entrez_ids_to_esets(
      esets = esets,
      entrezgene.db = db
    )
    Biobase::featureData(res)[["entrez.id"]]
  }

  # Logic tests (low level stuff is done by multisymbol_to_entrez_ids):
  # Blank refseq entries, NA refseq entries should map to NA
  eset_blank_refseqs <- eset_empty
  annotation(eset_blank_refseqs) <- "SOME_OTHER_PLATFORM"
  Biobase::featureData(eset_blank_refseqs) <- adf(
    data.frame("GB_ACC" = c(NA, ""), stringsAsFactors = FALSE)
  )
  expect_blank_refseqs <- rep(as.character(NA), 2)
  result_blank_refseqs <- result_helper(eset_blank_refseqs)
  expect_equal(
    expect_blank_refseqs,
    result_blank_refseqs,
    info = "NA and empty string should give NA entrez.ids"
  )

  # A single refseq entry that maps to a single human gene
  eset_single_refseq <- eset_empty
  annotation(eset_single_refseq) <- "platform.1"
  Biobase::featureData(eset_single_refseq) <- adf(
    data.frame(
      "GB_LIST" = "NM_000579", # CCR5/'1234'
      stringsAsFactors = FALSE
    )
  )
  expect_single_refseq <- "1234"
  result_single_refseq <- result_helper(eset_single_refseq)
  expect_equal(
    expect_single_refseq,
    result_single_refseq,
    info = "Single refseq id that maps to a single entrez id"
  )

  # A ' /// '-separated entry of refseq ids
  eset_two_refseq <- eset_empty
  annotation(eset_two_refseq) <- "GPL12345"
  Biobase::featureData(eset_two_refseq) <- adf(
    data.frame(
      "GB_LIST" = "NM_001307936 /// NM_018976",
      stringsAsFactors = FALSE
    )
  )
  expect_two_refseq <- "54407"
  result_two_refseq <- result_helper(eset_two_refseq)
  expect_equal(
    expect_two_refseq,
    result_two_refseq,
    info = "Two refseqs, ///-separated, that map to a single entrez id"
  )

  # A ','-separated entry of refseq ids
  eset_comma_refseq <- eset_empty
  annotation(eset_comma_refseq) <- "GPL9876"
  Biobase::featureData(eset_comma_refseq) <- adf(
    data.frame(
      "GB_LIST" = "NM_001130045,NM_153254,BC126152",
      stringsAsFactors = FALSE
    )
  )
  expect_comma_refseq <- "254173"
  result_comma_refseq <- result_helper(eset_comma_refseq)
  expect_equal(
    expect_comma_refseq,
    result_comma_refseq,
    info = "Comma-separated refseq list, that map to a single entrez id"
  )

  # refseq to entrez mapping for two different mouse probes:
  eset_mouse_refseq <- eset_empty
  annotation(eset_mouse_refseq) <- "GPL1261"
  Biobase::featureData(eset_mouse_refseq) <- adf(
    data.frame(
      "RefSeq Transcript ID" = c("NM_017477 /// NM_201244", "NM_013477"),
      stringsAsFactors = TRUE,
      check.names = FALSE
    )
  )
  expect_mouse_refseq <- c("54161", "11972")
  result_mouse_refseq <- result_helper(eset_mouse_refseq,
                                      org.Mm.eg.db::org.Mm.eg.db)
  expect_equal(
    expect_mouse_refseq,
    result_mouse_refseq,
    info = "refseq to entrez mapping for two mouse probes"
  )
  # Add entrez.ids to a swissprot-containing ESet
  eset_swissprot <- eset_empty
  annotation(eset_swissprot) <- "SOME GPL ID"
  Biobase::featureData(eset_swissprot) <- adf(
    data.frame(
      "swissprot" = c(NA, "", "---", "ENST0000412115", "NR_046018",
                      "NM_001005221"),
      stringsAsFactors = FALSE
    )
  )
  expect_swissprot <- c(NA, NA, NA, NA, "100287102", "729759")
  result_swissprot <- result_helper(eset_swissprot)
  expect_equal(
    expect_swissprot,
    result_swissprot,
    info = "swissprot column can be mapped to Entrez ids"
  )
})
