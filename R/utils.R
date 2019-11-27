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
