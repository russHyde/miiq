% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalisation_functions.R
\name{normalise_processed_eset}
\alias{normalise_processed_eset}
\title{Take a processed ESet, eg from GEOquery, and normalise it}
\usage{
normalise_processed_eset(eset)
}
\arguments{
\item{eset}{An ExpressionSet}
}
\description{
Note this is a duplication of gld_fnDefault_transformAndProbeFilter from
drug.markers (but uses quantile normalisation)
}
