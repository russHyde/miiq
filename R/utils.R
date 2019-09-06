#' Obtain an object from a specific package based on a string
#'
#' Suppose you provide x = "abc::def", then this will return the `def` function
#' or variable from the package `abc`. If x = "ghi", then the `ghi` function
#' or variable from the enclosing environment will be returned
#'
#' @param    x        a string in the form "some_object_name" or
#'   "some_package::some_object_name".
#'

.get_from_env <- function(x) {
  stopifnot(is.character(x) && length(x) == 1)

  env_fn <- strsplit(x, "::")[[1]]

  if (length(env_fn) == 1) {
    env <- -1
    fn <- env_fn
  } else {
    env <- paste0("package:", env_fn[1])
    fn <- env_fn[2]
  }

  get(fn, pos = env)
}
