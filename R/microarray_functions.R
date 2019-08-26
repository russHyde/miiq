###############################################################################

#' Function that performs a few tests to check if a given \code{ExpressionSet}
#' has already been log-transformed
#'
#' @param        eset          an \code{ExpressionSet}.
#' @param        range_limit_if_transformed   the largest possible range
#' between high and low values if the \code{ESet} is already transformed.
#' @param        skew_threshold   if the skewness of all columns of the
#' dataset is greater than this threshold, the dataset is assumed
#' untransformed.
#'
#' @include      vector_manipulation.R
#'
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @importFrom   Biobase      exprs
#' @importFrom   moments      skewness
#'
#' @export
#'
check_if_pretransformed_eset <- function(
                                         eset = NULL,
                                         range_limit_if_transformed = 1000,
                                         skew_threshold = 3) {
  # NOTE: refactored 20180223 - argnames are now underscore-separated
  # - This function is used by rnaseq_marray_merge and 2015_prelim, but those
  # projects have not yet been refactored to use the renamed args.
  # - If this script is imported back into rnaseq_marray_merge or 2015_prelim,
  # rewrite all calls to this function

  # Perform a few tests to check if a given eset has been log transformed:

  if (!exists("na_inf_omit")) {
    # how do I ensure that na_inf_omit is imported?
    stop("na_inf_omit is required by check_if_pretransformed")
  }

  # Check that the input is a valid 'eset'
  if (is.null(eset) || !is(eset, "ExpressionSet")) {
    stop("'eset' should be a valid Biobase::ExpressionSet")
  }

  # If there is no data in the eset, or all entries are NA,
  #   we consider the eset to have been transformed
  if (nrow(eset) == 0 || ncol(eset) == 0 || all(is.na(Biobase::exprs(eset)))) {
    return(TRUE)
  }
  # Determine if any negative values are present;
  #   assume the dataset has been log-transformed if it contains any negatives
  if (any(as.vector(Biobase::exprs(eset)) < 0, na.rm = TRUE)) {
    return(TRUE)
  }
  # Determine the range of expression values in the dataset
  rng <- range(
    na_inf_omit(Biobase::exprs(eset))
  )

  # Determine the skewness of each microarray in the dataset
  skews <- apply(Biobase::exprs(eset), 2, function(x) {
    moments::skewness(na_inf_omit(x))
  })

  # If the dataset is not transformed it will have high range,
  #   and will have skewed distribution of expression values in each array
  # If either of these occur, we assume the dataset is not transformed
  untransformed <- (rng[2] - rng[1]) > range_limit_if_transformed ||
    all(skews > skew_threshold)

  return(!untransformed)
}

###############################################################################

#' Not exported: checks if a valid refseq column is present in an
#'   \code{ExpressionSet} and returns a vector of all refseq column names if it
#'   does; otherwise returns an empty vector
#'
#' @param        gset          An \code{ExpressionSet}.
#' @param        refseq_col_types   The choice of refseq column names.
#'
#' @importClassesFrom   Biobase   ExpressionSet
#' @importFrom   Biobase       varLabels   featureData
#'
get_refseq_colnames <- function(
                                gset = NULL,
                                refseq_col_types = c(
                                  "RefSeq Transcript ID",
                                  # in case make.names has been applied:
                                  "RefSeq.Transcript.ID",
                                  "GB_LIST",
                                  "GB_ACC"
                                )) {
  # returns the names of the RefSeq columns in the gset if any exist
  # and returns an empty character vector otherwise
  stopifnot(!is.null(gset) && is(gset, "ExpressionSet"))

  refseq_cols <- intersect(
    refseq_col_types,
    Biobase::varLabels(Biobase::featureData(gset))
  )

  if (length(refseq_cols) > 0) {
    return(refseq_cols)
  } else {
    return(character(0))
  }
}

#' Not exported: checks if an ExpressionSet has a valid RefSeq column
#'
#' @inheritParams   get_refseq_colnames
#'
has_refseq_column <- function(
                              gset = NULL) {
  # Returns TRUE/FALSE depedning whether a RefSeq-esque column is present in
  #   the input gset
  refseq_cols <- get_refseq_colnames(gset = gset)
  length(refseq_cols) > 0
}

#' Returns a RefSeq-annotation-containing column name from an
#' \code{ExpressionSet}
#'
#' Extracts a single column name from an \code{ExpressionSet}. The
#' corresponding annotation column contains RefSeq IDs for the rows of that
#' \code{ExpressionSet}. Where multiple Refseq columns are present, only the
#' first is returned (order specified by \code{get_refseq_colnames}).
#'
#' Errors are thrown if no RefSeq column is present or if the RefSeq column
#' contains neither character nor factor variables.
#'
#' @inheritParams   get_refseq_colnames
#'
#' @return       A single string - the name of an annotation column that
#'   contains RefSeq data inside the \code{ExpressionSet}.
#'
#' @importFrom   Biobase       featureData
#' @export
#'

get_refseq_column <- function(
                              gset = NULL) {
  refseq_cols <- get_refseq_colnames(gset = gset)
  for (rc in refseq_cols) {
    fd <- Biobase::featureData(gset)[[rc]]
    stopifnot(is.character(fd) || is.factor(fd))
  }
  stopifnot(length(refseq_cols) >= 1)
  refseq_cols[1]
}

###############################################################################

#' Extract any NM_XXXX or NR_XXXX RefSeq/GenBank identifiers out of a SwissProt
#' column in an \code{ExpressionSet}
#'
#' Cells are delimited as follows: ABC /// DEF /// GHI1 // GHI2 /// JKL
#' So that main entries are ///-separated and subentries are //-separated.
#'
#' Returned value should be a vector of comma-separated refseq ids.
#'
#' @param        swiss        A vector of swissprot ids.
#'
#' @importFrom   magrittr      %>%
#'

swissprot_column_to_refseq <- function(swiss) {
  # TODO: decide if BC039241 - type ids should be kept as well
  if (is.null(swiss) ||
    length(swiss) == 0
  ) {
    stop("NULL input to swissprot_column_to_refseq")
  }

  grepval_refseq <- function(v) {
    gsub(
      pattern = "\\.[0-9]+$",
      replacement = "",
      x = grep("NR_|NM_", v, value = TRUE)
    )
  }

  # make a list of vectors of SwissProt IDs from the input vector
  list_of_swiss_ids <- gsub(
    pattern = " /// | // ", replacement = ",", x = swiss
  ) %>%
    strsplit(split = ",")

  # then collapse each vector into a deduplicated, comma-separated, string
  # and return a vector of comma-separated strings
  Map(grepval_refseq, list_of_swiss_ids) %>%
    lapply(function(x) paste(unique(x), collapse = ",")) %>%
    unlist()
}
