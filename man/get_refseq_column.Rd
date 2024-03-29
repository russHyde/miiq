% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microarray_functions.R
\name{get_refseq_column}
\alias{get_refseq_column}
\title{Returns a RefSeq-annotation-containing column name from an
\code{ExpressionSet}}
\usage{
get_refseq_column(gset = NULL)
}
\arguments{
\item{gset}{An \code{ExpressionSet}.}
}
\value{
A single string - the name of an annotation column that
  contains RefSeq data inside the \code{ExpressionSet}.
}
\description{
Extracts a single column name from an \code{ExpressionSet}. The
corresponding annotation column contains RefSeq IDs for the rows of that
\code{ExpressionSet}. Where multiple Refseq columns are present, only the
first is returned (order specified by \code{get_refseq_colnames}).
}
\details{
Errors are thrown if no RefSeq column is present or if the RefSeq column
contains neither character nor factor variables.
}
