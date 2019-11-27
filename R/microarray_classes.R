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

#' Row-subsetting operator for `eset_limma_dataset` class
#'
#' @param        x             An `eset_limma_dataset` object
#' @param        i             A subset of the rows in `x`
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



###############################################################################

## #' @import       GEOquery
## #' @import       Biobase
## #' @import       limma

###############################################################################

#' Builds a function that returns the col indices of the desired samples
#'
#' Creates a function that can be applied to a geo-limma-dataset
#'   or to an ExpressionSet, and returns the column indices of those
#'   samples that are to kept
#' This just reduces the boilerplate involved in checking g-l-d or
#'   ExpressionSet nature of the input
#'
#' @param        column.filter.fn   Function that decides which samples to keep
#'
#' @include      utils.R
#' @export
#'

gld_fnBuilder_keepSample <- function(
                                     # nolint start
                                     column.filter.fn
                                     # nolint end
) {
  # nolint start
  stopifnot(is.function(column.filter.fn))
  fn <- column.filter.fn
  # nolint end

  function(
             # nolint start
             geo.limma.dataset = NULL,
             # nolint end
             gset = NULL) {
    warning("deprecation warning: recommend using a function(gset)")
    # nolint start
    gld <- geo.limma.dataset
    # nolint end
    gset <- .check_or_get_eset(gld, gset)
    keep_cols <- fn(gset)
    keep_cols
  }
}

###############################################################################

#' Builds a function that makes the design matrix for a dataset
#'
#' User passes in the column labels of the phenoData that are to be used in
#'   making the design
#' And passes in a design.function that converts these columns into the design
#' The functionBuilder checks for the validity of the input (treatment.cols /
#'   design.function) and the built function checks the validity of the
#'   ExpressionSet / geo-limma-dataset that was passed in
#'
#' @param        treatment.cols   Which columns of the pData should be kept?
#' @param        design.fn     Function for making design matrices.
#'
#' @importFrom   Biobase       pData   sampleNames
#'
#' @include      utils.R
#' @export
#'

gld_fnBuilder_exptDesign <- function(
                                     treatment.cols = NULL,
                                     design.fn = NULL) {
  # check type-validity of the treatment.cols and design.fn
  stopifnot(is.character(treatment.cols) && length(treatment.cols) > 0)
  stopifnot(is.function(design.fn))

  design_function <- function(
                                # nolint start
                                geo.limma.dataset = NULL,
                                # nolint end
                                gset = NULL) {
    # nolint start
    geo_limma_dataset <- geo.limma.dataset
    # nolint end
    # check input validity of the Expressionset / geo_limma_dataset
    gset <- .check_or_get_eset(geo_limma_dataset, gset)

    # Extract the columns of pData that contain the design-related info
    # and then construct the design from these columns
    pdata <- Biobase::pData(gset)
    if (!all(treatment.cols %in% colnames(pdata))) {
      missing_cols <- setdiff(treatment.cols, colnames(pdata))
      stop(paste(
        "Some treatment.cols were missing from pData:",
        paste(missing_cols, collapse = "/")
      ))
    }
    treatments <- Biobase::pData(gset)[, treatment.cols]
    design <- as.data.frame(design.fn(treatments))
    sn <- Biobase::sampleNames(gset)
    cn <- make.names(colnames(design))

    if (length(sn) != nrow(design)) {
      message(
        c("sampleNames (ExpressionSet): ", paste(sn, collapse = " "))
      )
      message(
        c("sampleNames (design): ", paste(rownames(design), collapse = " "))
      )
      print(design)
      stop("length of sampleNames (in eset) and design do not agree")
    }

    rownames(design) <- sn
    colnames(design) <- cn
    design
  }

  design_function
}

###############################################################################

#' Default function for setting up an expt design from the pData of the
#'   ExpressionSet
#'
#' Just takes the "title" column of the pData and uses it as a factor.
#'
#' @param        geo.limma.dataset   A geo_limma_dataset.
#' @param        gset          An ExpressionSet.
#'
#' @export
#'

gld_fnDefault_exptDesign <- gld_fnBuilder_exptDesign(
  treatment.cols = "title",
  design.fn = function(treatments) {
    title <- treatments
    model.matrix(~title)
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
# - geo_config, geo_limma_dataset, limma_config
# - decideTests_gld, eset<-,
# - limma_workflow
# - gld_fnBuilder_exptContrasts
# - gld_fnBuilder_exptContrasts_fromList

# TODO:
# - keep class definitions, setters / getters / subsetters into
#   this file
# - move functions, pipelines etc into different files
