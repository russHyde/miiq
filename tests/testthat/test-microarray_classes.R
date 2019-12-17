###############################################################################

context("Tests for the functions in microarray_classes.R")

###############################################################################
# Expression sets that are used multiple times:
eset_empty <- Biobase::ExpressionSet()

eset_1x2 <- Biobase::ExpressionSet(
  # 1 probe, 2 samples
  assayData = matrix(
    1:2,
    nrow = 1,
    dimnames = list("probe1", c("samp1", "samp2"))
  ),
  phenoData = Biobase::AnnotatedDataFrame(data.frame(
    title = c("A", "B"),
    row.names = c("samp1", "samp2")
  ))
)

###############################################################################
# Test the validity of eset_limma_datasets
# - Any entry can be missing
# - Otherwise, the entries must be:
#     - ExpressionSet (eset),
#     - data.frame (design),
#     - matrix (contrast)
#     - MArrayLM (fits.init)
#     - MArrayLM (fits)
test_that("Validity of eset_limma_datasets (entries and dimensions)", {
  xset <- Biobase::ExpressionSet()
  dfr <- data.frame()
  mat <- matrix()
  mlm <- new("MArrayLM", "limma")
  incorrect_class <- "Not of the right class"

  expect_error(
    object = new(
      "eset_limma_dataset",
      eset = incorrect_class, design = dfr, contrast = mat, fits.init = mlm,
      fits = mlm
    ),
    info = "If eset is defined, it should be a subclass of ExpressionSet"
  )

  expect_error(
    object = new(
      "eset_limma_dataset",
      eset = xset, design = incorrect_class, contrast = mat, fits.init = mlm,
      fits = mlm
    ),
    info = "If design is defined, it should be a data.frame"
  )

  expect_error(
    object = new(
      "eset_limma_dataset",
      eset = xset, design = dfr, contrast = incorrect_class, fits.init = mlm,
      fits = mlm
    ),
    info = "If contrast is defined, it should be a matrix"
  )

  expect_error(
    object = new(
      "eset_limma_dataset",
      eset = xset, design = dfr, contrast = mat, fits.init = incorrect_class,
      fits = mlm
    ),
    info = "If fits.init is defined, it should be a subclass of MArrayLM"
  )

  expect_error(
    object = new(
      "eset_limma_dataset",
      eset = xset, design = dfr, contrast = mat, fits.init = mlm,
      fits = incorrect_class
    ),
    info = "If fits.init is defined, it should be a subclass of MArrayLM"
  )

  # If both eset and fits are non-zero, they should have the same number of
  # rows
  eset_3x2 <- Biobase::ExpressionSet(
    assayData = matrix(nrow = 3, ncol = 2)
  )
  mlm_2x1 <- new(
    "MArrayLM",
    list(coefficients = matrix(nrow = 2, ncol = 1))
  )
  expect_error(
    object = eset_limma_dataset(eset = eset_3x2, fits = mlm_2x1),
    info = "If both eset and fits are non-zero, they should have the same
number of rows (probes)"
  )
  # If both eset and design are non-zero, they should have the same nummber of
  # samples, (ncol(eset) == nrow(design))
  design_4x1 <- data.frame(single.column = rep(0, 4))
  expect_error(
    object = eset_limma_dataset(eset = eset_3x2, design = design_4x1),
    info = paste(
      "If both eset and design are non-zero,",
      "they should have the same number of samples:",
      "2 samples vs 4 samples"
    )
  )
  design_1x1 <- data.frame(single.column = 1)
  expect_error(
    object = eset_limma_dataset(eset = eset_3x2, design = design_1x1),
    info = paste(
      "If both eset and design are non-zero,",
      "they should have the same number of samples:",
      "2 samples in eset vs 1x1 data.frame as design"
    )
  )

  # nb, setClass("MArrayLM", representation("list"), where list has entries:
  #  coefficients:    Matrix containing fitted coefficients or contrasts.
  #  stdev.unscaled:  Matrix containing unscaled standard deviations of the
  #    coefficients or contrasts.
  #  sigma:           Numeric vector containing residual standard deviations
  #    for each gene
  #  df.residual:     Numeric vector containing residual degrees of freedom for
  #    each gene.
})

test_that("Row subsetting of eset_limma_datasets works as expected", {

  # Test data: 3 probes, 2 samples, 2 coefficients, 1 contrast
  eset_3x2 <- random_eset(n_probes = 3, n_samples = 2)
  design_2x2 <- data.frame(single.column = rep(1, 2))
  contrast_2x1 <- matrix(0, nrow = 2, ncol = 1)
  fit_3x2 <- new(
    "MArrayLM",
    list(coefficients = matrix(0, nrow = 3, ncol = 2))
  )

  expect_equal(
    object = eset_limma_dataset(eset = eset_3x2)[1:2, ],
    expected = eset_limma_dataset(eset = eset_3x2[1:2, ]),
    info = "Subsetting an eset_limma_dataset that only contains eset -
affects rows of eset"
  )

  expect_equal(
    object = eset_limma_dataset(design = design_2x2)[1:2, ],
    expected = eset_limma_dataset(design = design_2x2),
    info = "Subsetting an eset_limma_dataset that only contains design,
should not affect anything"
  )

  expect_equal(
    object = eset_limma_dataset(contrast = contrast_2x1)[1:2, ],
    expected = eset_limma_dataset(contrast = contrast_2x1),
    info = "Subsetting an eset_limma_dataset that only contains contrast,
should not affect anything"
  )

  expect_equal(
    object = eset_limma_dataset(fits = fit_3x2)[1:2, ],
    expected = eset_limma_dataset(fits = fit_3x2[1:2, ]),
    info = "Subsetting an ELD that only contains fits - affects rows of
fits"
  )

  expect_equal(
    object = eset_limma_dataset(fits.init = fit_3x2)[1:2, ],
    expected = eset_limma_dataset(fits.init = fit_3x2[1:2, ]),
    info = "Subsetting an ELD that only contains fits.init - affects rows
of fits.init"
  )
})

###############################################################################

test_that("Row subsetting of eset_limma_dataset using rownames", {
  n_probes <- 20
  n_coefs <- 5
  probe_names <- paste0("P", seq_len(n_probes))
  coef_names <- paste0("C", seq_len(n_coefs))

  eset <- random_eset(n_probes = n_probes)
  rownames(eset) <- paste0("P", seq_len(n_probes))

  fit <- new(
    "MArrayLM",
    list(
      coefficients = matrix(
        0,
        nrow = n_probes, ncol = n_coefs,
        dimnames = list(probe_names, coef_names)
      )
    )
  )

  index_vec <- sample(probe_names, sample.int(n_probes, 1))

  expect_equal(
    object = eset_limma_dataset(eset = eset)[index_vec, ],
    expected = eset_limma_dataset(eset = eset[index_vec, ]),
    info = "Subsetting an ELD that only contains an ESet, by rownames"
  )

  expect_equal(
    object = eset_limma_dataset(fits = fit)[index_vec, ],
    expected = eset_limma_dataset(fits = fit[index_vec, ]),
    info = "Subsetting an ELD that only contains a fits entry, by rownames"
  )

  expect_equal(
    object = eset_limma_dataset(fits.init = fit)[index_vec, ],
    expected = eset_limma_dataset(fits.init = fit[index_vec, ]),
    info = "Subsetting an ELD that only contains a fits.init entry, by rowname"
  )
})

###############################################################################

test_that("Default keep_sample_fn - keep all samples; die if no samples", {
  # TODO: remove these since they duplicate tests in `test-filter_functions.R`
  # or replace them with a test on `run_preprocess_workflow`
  expect_error(
    object = keep_all_samples(gset = eset_empty),
    info = "Default keep.sample.fn should crash on an empty eset"
  )

  eset_1x1 <- Biobase::ExpressionSet(assayData = matrix(1))
  expect_equal(
    object = keep_all_samples(gset = eset_1x1),
    expected = c(1),
    info = "Default keep_sample_fn should keep all cols of a non-empty eset"
  )
})

###############################################################################

test_that("Default exptDesign function - use title as a factor", {
  expect_error(
    object = gld_fnDefault_exptDesign(gset = eset_empty),
    info = "Default exptDesign should crash on an empty eset"
  )

  eset_without_pheno_data <- Biobase::ExpressionSet(assayData = matrix(1))
  expect_error(
    object = gld_fnDefault_exptDesign(gset = eset_without_pheno_data),
    info = "Default exptDesign should crash if no pData entries"
  )

  eset_without_title_column <- Biobase::ExpressionSet(
    assayData = matrix(1),
    pData = Biobase::AnnotatedDataFrame(data.frame(NOT.A.TITLE = "A"))
  )
  expect_error(
    object = gld_fnDefault_exptDesign(gset = eset_without_title_column),
    info = "Default exptDesign should crash if no `title` column in pData"
  )

  eset_1level <- Biobase::ExpressionSet(
    assayData = matrix(1),
    pData = Biobase::AnnotatedDataFrame(data.frame(title = "A"))
  )
  expect_error(
    object = gld_fnDefault_exptDesign(gset = eset_1level),
    info = "Default exptDesign should crash if only a single `title` level
in pData"
  )

  expect_equal(
    object = gld_fnDefault_exptDesign(gset = eset_1x2) %>% as.data.frame(),
    expect = data.frame(
      "X.Intercept." = c(1, 1), # note use of make.names
      "titleB" = c(0, 1),
      row.names = c("samp1", "samp2")
    ),
    info = "Default exptDesign - two levels of `title`"
  )
})

###############################################################################

test_that("Default keepProbe function", {
  # TODO: remove these tests since they duplicate tests in
  # "test-filter_functions.R"; or replace them with a test on
  # `run_preprocess_workflow`
  expect_error(
    object = keep_all_probes(gset = eset_empty),
    info = "Default keepProbe function: crash if there are no probes"
  )

  expect_equal(
    object = keep_all_probes(gset = eset_1x2),
    expected = 1,
    info = "Default keepProbe function: keep all probes"
  )
})

###############################################################################
