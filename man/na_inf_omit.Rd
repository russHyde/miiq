% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vector_manipulation.R
\name{na_inf_omit}
\alias{na_inf_omit}
\title{Function to remove any NA or infinite values from a vector}
\usage{
na_inf_omit(x)
}
\arguments{
\item{x}{Some vector - only makes sense if \code{x} is
\code{numeric}.}
}
\description{
The returned vector will be shorter than the input vector if any values
are dropped.
}
