###############################################################################

#' @name         nullFunction
#' @title        nullFunction
#' @importFrom   methods       setValidity
#' @importFrom   methods       setClass
#' @importFrom   methods       setGeneric
#' @importFrom   methods       setMethod
#' @importFrom   methods       setReplaceMethod
#'

NULL

###############################################################################

#' validity function for eset_limma_dataset objects
#'
#' @noRd
#'

.validity_eld <- function(object) {
  valid <- TRUE
  msgs <- c()

  # Check that the number of probes is consistent between eset, fits and
  #   fits.init
  check_probes <- function(x) {
    non_empty_nrows <- unlist(
      Map(
        nrow,
        Filter(function(x) nrow(x) > 0, x)
      )
    )

    if (length(non_empty_nrows) > 0) {
      return(all(non_empty_nrows == non_empty_nrows[1]))
    } else {
      return(TRUE)
    }
  }

  if (
    !check_probes(list(object@eset, object@fits.init, object@fits))
  ) {
    valid <- FALSE
    msgs <- c(
      msgs,
      paste(
        "All non-empty entries in {eset, fits, fits.init} should have",
        "the same number of probes (rows)"
      )
    )
  }

  # Check that the number of samples is consistent between eset and design
  check_samples <- function(eset, design) {
    if (ncol(eset) == 0 || nrow(design) == 0) {
      return(TRUE)
    }
    return(ncol(eset) == nrow(design))
  }

  if (!check_samples(object@eset, object@design)) {
    valid <- FALSE
    msgs <- c(
      msgs,
      paste(
        "All non-empty entries in {eset, design} should have the same",
        "number of samples"
      )
    )
  }

  #
  if (!valid) {
    return(msgs)
  } else {
    return(TRUE)
  }
}

#' @title        Class to hold ExpressionSets and limma model results
#'
#' @importClassesFrom   limma   MArrayLM
#' @importClassesFrom   Biobase   ExpressionSet
#' @importFrom   methods       new
#'
#' @name         eset_limma_dataset-class
#' @rdname       eset_limma_dataset-class
#'
#' @exportClass   eset_limma_dataset
#'

methods::setClass(
  "eset_limma_dataset",
  slots = list(
    eset = "ExpressionSet",
    design = "data.frame",
    contrast = "matrix",
    fits.init = "MArrayLM",
    fits = "MArrayLM"
  ),
  prototype = list(
    eset = Biobase::ExpressionSet(),
    design = data.frame(),
    contrast = matrix(nrow = 0, ncol = 0),
    fits.init = methods::new("MArrayLM"),
    fits = methods::new("MArrayLM")
  ),
  validity = function(object) .validity_eld(object)
)

#' wrapper function eset_limma_dataset
#'
#' @name         eset_limma_dataset
#' @rdname       eset_limma_dataset-class
#'
#' @param        ...           arguments passed to the constructor for
#'   `eset_limma_dataset`.
#'
#' @importFrom   methods       new
#'
#' @export
#'

eset_limma_dataset <- function(...) {
  methods::new("eset_limma_dataset", ...)
}

###############################################################################

# eset_limma_dataset methods:

setGeneric("feature_names", function(x) {
  standardGeneric("feature_names")
})

methods::setMethod(
  "feature_names",
  signature = c(x = "eset_limma_dataset"),
  definition = function(x) {
    if (nrow(x@eset) > 0) {
      rownames(x@eset)
    } else if (nrow(x@fits) > 0) {
      rownames(x@fits)
    } else {
      rownames(x@fits.init)
    }
  }
)

#' Row-subsetting (feature-subsetting) operator for `eset_limma_dataset` class
#'
#' @param        x             An `eset_limma_dataset` object
#' @param        i             A subset of the row numbers in `x` the eset,
#'   fits or fits.init of x (that is, a subset of the features in x)
#' @param        j,drop,...    Not used
#'
#' @name         [
#' @rdname       eset_limma_dataset
#' @aliases      [,eset_limma_dataset,numeric,missing,missing-method
#'

methods::setMethod(
  "[",
  signature = c(
    x = "eset_limma_dataset",
    i = "numeric",
    j = "missing",
    drop = "missing"
  ),
  definition = function(x, i, j, ..., drop) {
    row_subsetter <- function(eset_or_fit) {
      if (nrow(eset_or_fit) > 0) {
        return(eset_or_fit[i, ])
      } else {
        return(eset_or_fit)
      }
    }

    # eset, fits and fits.init should be restrictied to the requested rows
    .eset <- row_subsetter(x@eset)
    .fits <- row_subsetter(x@fits)
    .fits_init <- row_subsetter(x@fits.init)

    # design & contrast are unaffected by subsetting the number of probes
    .design <- x@design
    .contrast <- x@contrast

    eset_limma_dataset(
      eset = .eset,
      design = .design,
      contrast = .contrast,
      fits = .fits,
      fits.init = .fits_init
    )
  }
)

#' Row-subsetting (feature-subsetting) operator for `eset_limma_dataset` class
#'
#' @param        x             An `eset_limma_dataset` object
#' @param        i             A subset of the feature_names in `x` (the
#'   rownames of the eset, fits or fits.init entries)
#' @param        j,drop,...    Not used
#'
#' @name         [
#' @rdname       eset_limma_dataset
#' @aliases      [,eset_limma_dataset,character,missing,missing-method
#'

methods::setMethod(
  "[",
  signature = c(
    x = "eset_limma_dataset",
    i = "character",
    j = "missing",
    drop = "missing"
  ),
  definition = function(x, i, j, ..., drop) {
    keep_features <- match(i, feature_names(x))
    x[keep_features, ]
  }
)

###############################################################################

#' Default function for annotating the ExpressionSet with additional info
#'
#' This just adds entrez.ids.
#' Not exported.
#' Not sure how to import entrezgene database (entrezgene.db.id) by name
#'   when using \<AT\>import annotations.
#'
#' @param        geo_limma_dataset   A geo_limma_dataset.
#' @param        gset          An ExpressionSet.
#' @param        entrezgene_db   The AnnotationDBI (or it's name) for use in
#'   mapping to entrez ids.
#'
#' @importClassesFrom   AnnotationDbi   OrgDb
#' @include      utils.R
#'
gld_fnDefault_esetAnnotation <- function(
  geo_limma_dataset = NULL,
  gset = NULL,
  entrezgene_db = NULL
) {
  # Function can either use geo_limma_datasets or a combination of
  # ExpressionSet and an AnnotationDBI
  gset <- .check_or_get_eset(geo_limma_dataset, gset)

  if (!is(entrezgene_db, "OrgDb")
      &&
      !is.character(entrezgene_db)
  ) {
    stopifnot(is(geo_limma_dataset, "geo_limma_dataset"))
    entrezgene_db <- geo_limma_dataset@geo.config@entrezgene.db.id
  }

  # entrezgene_db is checked by symbol_to_entrez_id, so checking it"s
  #   validity is unnecessary

  # Function to annotate an ExpressionSet with entrez.ids
  eset <- add_entrez_ids_to_esets(
    gset,
    entrezgene.db = entrezgene_db
  )

  return(eset)
}

###############################################################################

#' Default function for transforming (etc) a microarray dataset
#'
#' Not exported
#'
#' @param        geo.limma.dataset   A geo_limma_dataset
#' @param        gset          An ExpressionSet - used preferentially to
#'   GLD@eset.
#'
#' @include      utils.R
#' @return       An ExpressionSet
#'

gld_fnDefault_transformAndProbeFilter <- function(
  # nolint start
  geo.limma.dataset = NULL,
  # nolint end
  gset = NULL) {
  # nolint start
  geo_limma_dataset <- geo.limma.dataset
  # nolint end

  gset <- .check_or_get_eset(geo_limma_dataset, gset)

  # Determine if the imported ExpressionSet has already been
  #   log-transformed
  is_transformed <- check_if_pretransformed_eset(
    eset = gset
  )

  # Transform (log2) and median-Normalise the dataset
  # - also do some housekeeping:
  #   - dropping identical rows;
  #   - converting Inf to NA;
  #   - dropping rows with many missing values;
  # By default the array-wide median is subtracted, but the samples
  #   are not divided by the IQR (default chagned 14/7/2016)
  #   since IQR.normalisation affects the fold-changes
  filter_and_transform_eset(
    eset = gset,
    log2_transform = !is_transformed,
    drop_row_if_duplicated = TRUE,
    convert_inf_to_na = TRUE,
    normalise_method = "median",
    drop_row_na_inf_threshold = 0.25
  )
}

###############################################################################

#' Workflow function for converting an ExpressionSet into an entrez.id
#'   annotated ExpressionSet that has been log transformed / quantile
#'   normalised and median-normalised. Probes that do not map to entrez.id are
#'   dropped by default. All annotation / filtering functions can be overridden
#'
#' @param        gset          An ExpressionSet (must have genbank ids or
#'   swissprot ids or similar for default annotation function).
#' @param        entrezgene_db   character or \code{OrgDb}. Which
#'   \code{AnnotationDbi} object to use in the eset annotation function.
#' @param        eset_annot_fn   Function that annotates an input ExpressionSet
#'   and outputs an ExpressionSet.
#' @param        keep_sample_fn   Function that returns integer indices of
#'   samples (cols) of the ExpressionSet that are to be kept.
#' @param        transform_and_filter_fn   Function that applies various
#'   filtering / transformation steps to the values in the input ExpressionSet
#'   and returns an ExpressionSet
#' @param        keep_probe_fn   Function that returns integer indices of
#'   probes (rows) of the ExpressionSet that are to be kept.
#'
#' @return       An ExpressionSet
#'
#' @importFrom   magrittr      %>%
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @include   filter_functions.R
#' @export
#'
preprocess_eset_workflow <- function(
  gset = NULL,
  entrezgene_db = NULL,
  eset_annot_fn = gld_fnDefault_esetAnnotation,
  keep_sample_fn = keep_all_samples,
  transform_and_filter_fn = gld_fnDefault_transformAndProbeFilter,
  keep_probe_fn = keep_all_probes
) {
  # validity checks
  stopifnot(is(gset, "ExpressionSet"))
  # the entrezgene.db will be checked in Entrez-annotation functions
  if (is.null(entrezgene_db)) {
    stop("No `entrezgene_db` provided in `preprocess_eset_workflow`")
  }

  lambda_sample_filter <- function(gset) {
    keep_samples <- keep_sample_fn(gset = gset)
    gset[, keep_samples]
  }

  lambda_probe_filter <- function(gset) {
    keep_probes <- keep_probe_fn(gset = gset)
    gset[keep_probes, ]
  }
  # eset annotation
  #
  annotated_eset <- eset_annot_fn(
    gset = gset,
    entrezgene_db = entrezgene_db
  )

  # sample filtering
  # transform and filter probes
  # probe filtering
  transform_and_filter_fn(
    gset = lambda_sample_filter(annotated_eset)
  ) %>%
    lambda_probe_filter()

  # return the modified eset
}

# TODO: add tests / classes / functions:
# - geo_config, geo_limma_dataset
# - decideTests_gld, eset<-,

# TODO:
# - keep class definitions, setters / getters / subsetters into
#   this file
# - move functions, pipelines etc into different files
