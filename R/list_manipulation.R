#' get_or_else
#'
#' Given a list .x, this function looks to see if the field .field is defined
#' within that list. If it is defined, the value for that entry is returned. If
#' .field is not defined, then the function either returns the value
#' .default_val or the value .default_fn(.field). The latter is lazily
#' evaluated, ie, it's only computed if .field is not in .x.
#'
#' @param        .x            A list
#'
#' @param        .field        A string. This may or may not be one of the
#' names in the list .x; if it is the associated value is returned. If not
#' either .default_val or .default_fn(.field) is returned.
#'
#' @param        .default_val   A value that is to be returned if .field is not
#' in .x. Only one of .default_val and .default_fn may be defined.
#'
#' @param        .default_fn   A single argument function that is evaluated if
#' .field is not in .x. The value .field is passed into this function as an
#' argument.  Only one of .default_val and .default_fn may be defined.
#'
#' @export
#'
get_or_else <- function(
                        .x,
                        .field,
                        .default_val,
                        .default_fn) {
  # If .field is not in the names of .x, return .default, when it is defined,
  # or .default_fn(.field) when .default_fn is defined

  # Numeric field names are not used at present

  # .x should be a list or a NULL node (in which case it's treated as an empty
  # list)
  if (missing(.x) ||
    !(is.list(.x) || is.null(.x))
  ) {
    stop(".x should be a defined list or NULL in get_or_else")
  }
  if (missing(.field) ||
    !(is.character(.field))
  ) {
    stop(".field should be defined in get_or_else")
  }
  if (missing(.default_val) && missing(.default_fn)) {
    stop(".default_val or .default_fn should be defined in get_or_else")
  }
  if (!missing(.default_val) && !missing(.default_fn)) {
    stop(".default_val and .default_fn should not both be specified")
  }

  if (.field %in% names(.x)) {
    return(.x[[.field]])
  }

  if (!missing(.default_val)) {
    return(.default_val)
  } else {
    return(.default_fn(.field))
  }
}
