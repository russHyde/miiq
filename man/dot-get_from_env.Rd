% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.get_from_env}
\alias{.get_from_env}
\title{Obtain an object from a specific package based on a string}
\usage{
.get_from_env(x)
}
\arguments{
\item{x}{a string in the form "some_object_name" or
"some_package::some_object_name". Or, a variable (in which case the
variable is just returned; this allows the user to not need to know
whether they are working with a variable or the variable name).}
}
\description{
Suppose you provide x = "abc::def", then this will return the `def` function
or object from the package `abc`. If x = "ghi", then the `ghi` object from
the enclosing environment will be returned. If the actual variable `ghi`,
rather than it's name, is passed in it will be returned unmodified.
}
