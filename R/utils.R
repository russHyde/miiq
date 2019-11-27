###############################################################################

#' Obtain an object from a specific package based on a string
#'
#' Suppose you provide x = "abc::def", then this will return the `def` function
#' or object from the package `abc`. If x = "ghi", then the `ghi` object from
#' the enclosing environment will be returned. If the actual variable `ghi`,
#' rather than it's name, is passed in it will be returned unmodified.
#'
#' @param    x        a string in the form "some_object_name" or
#'   "some_package::some_object_name". Or, a variable (in which case the
#'   variable is just returned; this allows the user to not need to know
#'   whether they are working with a variable or the variable name).
#'

.get_from_env <- function(x) {
  if (is.function(x) || is(x, "OrgDb")) {
    return(x)
  }

  stopifnot(is.character(x) & length(x) == 1)

  env_fn <- strsplit(x, "::")[[1]]

  if (length(env_fn) == 1) {
    get(env_fn)
  } else {
    getExportedValue(env_fn[1], env_fn[2])
  }
}

###############################################################################

#' If no geo_limma_dataset is provided, or the provided one has no eset,
#' returns the gset
#'
#' Checks that the gset is an ExpressionSet
#'
#' @importFrom   methods       is
#'
#' @importClassesFrom   Biobase   ExpressionSet
#'
#' @noRd
#'

.check_or_get_eset <- function(
                               geo_limma_dataset = NULL,
                               gset = NULL) {
  if (is.null(gset)) {
    stopifnot(methods::is(geo_limma_dataset, "eset_limma_dataset"))
    gset <- geo_limma_dataset@eset
  }
  stopifnot(methods::is(gset, "ExpressionSet"))
  gset
}

###############################################################################
