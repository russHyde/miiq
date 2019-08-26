###############################################################################

#' Function that filters and transforms the entries in the \code{exprs} entry
#' of an \code{ExpressionSet}
#'
#' @param        eset          An \code{ExpressionSet}.
#' @param        log2_transform   BOOLEAN : Should the \code{ExpressionSet} be
#'   log2-transformed?
#' @param        drop_row_if_duplicated   BOOLEAN : Should rows be dropped if
#'   they are (numerical) duplicates of another row?
#' @param        drop_row_if_zero_variance   BOOLEAN : Should rows be dropped
#'   if they are constant?
#' @param        convert_inf_to_na   BOOLEAN : Should infinite vals be
#'   converted to NA?
#' @param        normalise_method   String from "quantile" (quantile-normalise
#'   the dataset, the default), "median" (subtract the medians off),
#'   "median_iqr" (subtract off the medians and divide by the iqr) or "none"
#'   (leave the dataset untouched).
#' @param        drop_row_na_inf_threshold   What fraction of missing (NA) or
#'   infinite values are permitted in a row before that row is dropped from
#'   further analysis.
#'
#' @importFrom   Biobase       exprs   exprs<-
#' @importFrom   stats         sd   median   IQR
#' @importFrom   preprocessCore   normalize.quantiles
#'
#' @export
#'

filter_and_transform_eset <- function(
                                      eset = NULL,
                                      # -- options
                                      log2_transform = TRUE,
                                      drop_row_if_duplicated = TRUE,
                                      drop_row_if_zero_variance = TRUE,
                                      convert_inf_to_na = TRUE,
                                      normalise_method = c(
                                        "quantile", "median", "median_iqr",
                                        "none"
                                      ),
                                      drop_row_na_inf_threshold = 0.25) {
  normalise_method <- match.arg(normalise_method)

  # Function to filter / transform the entries in the exprs entry of
  #   an ExpressionSet

  # Check that the input is a valid 'eset'
  if (is.null(eset) | !is(eset, "ExpressionSet")) {
    stop("'eset' should be a valid Biobase::ExpressionSet")
  }

  if (nrow(eset) == 0 || ncol(eset) == 0) {
    return(eset)
  }

  # Log2 normalise the dataset
  if (log2_transform) {
    if (any(exprs(eset) < 0, na.rm = TRUE)) {
      stop("Attempt to log transform negatives in eset")
    }
    Biobase::exprs(eset) <- log2(exprs(eset))
  }
  # Drop duplicated rows
  if (drop_row_if_duplicated) {
    keep_rows <- rownames(unique(exprs(eset)))
    eset <- eset[keep_rows, ]
  }
  # Drop rows that have a variance of zero
  # TODO: - Should this be done after transformations have been performed?
  if (drop_row_if_zero_variance) {
    sds <- apply(exprs(eset), 1, sd)
    keep_rows <- which(sds > 0)
    eset <- eset[keep_rows, ]
  }

  # Normalisation step:
  normalisation_functions <- list(
    "none" = identity,
    "quantile" = preprocessCore::normalize.quantiles
  )
  normalisation_functions$median <- function(xs) {
    medians <- apply(xs, 2, function(x) median(na_inf_omit(x)))
    sweep(xs, MARGIN = 2, STATS = medians, FUN = "-")
  }
  normalisation_functions$median_iqr <- function(xs) {
    ys <- normalisation_functions$median(xs)
    iqrs <- apply(ys, 2, function(x) IQR(na_inf_omit(x)))
    sweep(ys, MARGIN = 2, STATS = iqrs, FUN = "/")
  }

  Biobase::exprs(eset) <- normalisation_functions[[
  normalise_method
  ]](Biobase::exprs(eset))

  # Convert infinite values to NA
  if (convert_inf_to_na) {
    inf_idx <- which(is.infinite(exprs(eset)))
    Biobase::exprs(eset)[inf_idx] <- NA
  }

  # Drop any probes that have too many NA or Inf values in the dataset
  # Each probe should have <= 100*drop.row...threshold % NA or Inf
  na_inf_count <- apply(Biobase::exprs(eset), 1, function(x) {
    length(which(is.na(x) | is.infinite(x)))
  })
  keep_probes <- which(
    na_inf_count <= (drop_row_na_inf_threshold * ncol(eset))
  )
  eset <- eset[keep_probes, ]

  # Drop rows with an NA or Inf fraction > threshold
  # Quantile normalise within each array to ensure median = 0 and iqr = 1
  return(eset)
}

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
  # Perform a few tests to check if a given eset has been log transformed:

  # Check that the input is a valid 'eset'
  if (is.null(eset) || !is(eset, "ExpressionSet")) {
    stop("'eset' should be a valid Biobase::ExpressionSet")
  }

  # If there is no data in the eset, or all entries are NA,
  #   we consider the eset to have been transformed
  if (nrow(eset) == 0 || ncol(eset) == 0 || all(is.na(Biobase::exprs(eset)))) {
    return(TRUE)
  }

  xs <- Biobase::exprs(eset)

  # Determine if any negative values are present;
  #   assume the dataset has been log-transformed if it contains any negatives
  if (any(as.vector(xs) < 0, na.rm = TRUE)) {
    return(TRUE)
  }
  # Determine the range of expression values in the dataset
  rng <- range(
    na_inf_omit(xs)
  )

  # Determine the skewness of each microarray in the dataset
  skews <- apply(xs, 2, function(x) {
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

###############################################################################
