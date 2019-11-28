###############################################################################

# -- Utils for setting up datasets

random_eset <- function(
                        n_probes = sample.int(50, 1),
                        n_samples = sample.int(50, 1),
                        features = NULL) {
  values <- matrix(
    rnorm(n_probes * n_samples), nrow = n_probes, ncol = n_samples
  )
  features <- if(is.null(features)) {
    Biobase::annotatedDataFrameFrom(values, byrow = TRUE)
  } else {
    Biobase::AnnotatedDataFrame(features)
  }
  Biobase::ExpressionSet(values, featureData = features)
}

###############################################################################
