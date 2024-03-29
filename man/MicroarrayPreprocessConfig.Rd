% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_classes.R
\name{MicroarrayPreprocessConfig}
\alias{MicroarrayPreprocessConfig}
\title{Sets up a class for storing the options for preprocessing ExpressionSets
(that is, filtering the samples and probes; adding entrez and symbol
columns)}
\usage{
MicroarrayPreprocessConfig(acc, entrezgene_db, keep_sample_fn = NULL,
  keep_probe_fn = NULL, annot_gpl = as.logical(NA))
}
\arguments{
\item{acc}{An identifier for the microarray dataset.
Typically this would be the GEO or ArrayExpress identifier.}

\item{entrezgene_db}{A string or an \code{OrgDB} defining an
entrezgene database.}

\item{keep_sample_fn, keep_probe_fn}{Functions that return a vector
of integer indices that can be used to subset an ExpressionSet. By
default, all samples and features will be kept. Can be specified as a
function literal or a "function_name" or a "pkg_name::function_name"
string. If NULL, the default functions will keep all samples or keep all
probes.}

\item{annot_gpl}{Logical. Have the features in the dataset been
annotated with a GPL?}
}
\description{
Sets up a class for storing the options for preprocessing ExpressionSets
(that is, filtering the samples and probes; adding entrez and symbol
columns)
}
