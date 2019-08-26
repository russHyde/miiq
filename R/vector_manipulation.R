#' Function to remove any infinite values from a vector
#'
#' The returned vector will be shorter than the input vector if any infinite
#' values are found.
#'
#' @param        x             Some vector - only makes sense if \code{x} is
#'   \code{numeric}.
#'
#' @export

inf_omit <- function(x) {
  x[!is.infinite(x)]
}

#' Function to remove any NA or infinite values from a vector
#'
#' The returned vector will be shorter than the input vector if any values
#' are dropped.
#'
#' @inheritParams   inf_omit
#'
#' @importFrom   stats         na.omit
#'
#' @export

na_inf_omit <- function(x) {
  without_infs <- inf_omit(x)
  if (is.null(without_infs)) {
    without_infs
  } else {
    as.vector(na.omit(without_infs))
  }
}

###############################################################################
