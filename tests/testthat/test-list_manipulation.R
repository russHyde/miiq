context("Tests for list manipulation functions")

test_that("get_or_else", {
  expect_error(
    object = get_or_else(),
    info = "No input to get_or_else: should be a list"
  )
  expect_error(
    object = get_or_else("Not a list"),
    info = "Non-list input to get_or_else: should be a list"
  )
  expect_error(
    object = get_or_else(.x = list(a = 123)),
    info = "No fieldname to get_or_else: should be a string or index"
  )
  expect_error(
    object = get_or_else(.x = list(a = 123), .field = data.frame()),
    info = "Inappropriate fieldname to get_or_else: should be a string or
index"
  )
  expect_error(
    object = get_or_else(.x = list(), .field = "a"),
    info = "No .default_value or .default_fn in get_or_else"
  )

  expect_error(
    object = get_or_else(
      .x = list(a = 123, b = NA),
      .field = 1,
      .default_val = NA
    ),
    info = "Numeric fieldnames are not used yet in get_or_else"
  )

  expect_equal(
    object = get_or_else(
      .x = list(a = 123, b = NA),
      .field = "b",
      .default_val = NA
    ),
    expected = NA,
    info = "String fieldname to get_or_else"
  )
  expect_equal(
    object = get_or_else(
      .x = NULL,
      .field = "b",
      .default_val = "some.default"
    ),
    expected = "some.default",
    info = "NULL nodes can be used in get_or_else as if empty lists"
  )

  expect_equal(
    object = get_or_else(
      .x = list(a = 123, b = NA),
      .field = "b",
      .default_fn = identity
    ),
    expected = NA,
    info = "String fieldname with.default_fn to get_or_else"
  )

  expect_equal(
    object = get_or_else(
      .x = list(a = 123, b = NA),
      .field = "c",
      .default_val = "some.default"
    ),
    expected = "some.default",
    info = "String Field is missing in get_or_else; using .default_val"
  )
  expect_equal(
    object = get_or_else(
      .x = list(a = c(1:3), b = letters[3:4]),
      .field = "c",
      .default_val = c("Some", "default")
    ),
    expected = c("Some", "default"),
    info = "Vector valued default in get_or_else"
  )

  expect_error(
    object = get_or_else(
      .x = list(),
      .field = "a",
      .default_val = NA,
      .default_fn = identity
    ),
    info = "User should specify either .default_val or .default_fn but not
both"
  )
  expect_equal(
    object = get_or_else(
      .x = list(a = 123, b = NA),
      .field = "c",
      .default_fn = function(x) rep(x, 3)
    ),
    expected = c("c", "c", "c"),
    info = "String Field is missing in get_or_else; using .default_fn"
  )
})
