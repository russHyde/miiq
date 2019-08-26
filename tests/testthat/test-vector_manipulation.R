###############################################################################

context("Tests for vector manipulation functions")

###############################################################################

test_that("inf_omit", {
  expect_equal(
    inf_omit(NULL),
    NULL,
    info = "NULL input --> NULL output"
  )

  expect_equal(
    inf_omit(numeric(0)),
    numeric(0),
    info = "Empty input --> empty output"
  )

  expect_equal(
    inf_omit(Inf),
    numeric(0),
    info = "Single Inf --> empty output"
  )

  expect_equal(
    inf_omit(c(1, 2, Inf, 3, -Inf, 1)),
    c(1, 2, 3, 1),
    info = "mix of numbers and Infs"
  )

  expect_equal(
    inf_omit(as.numeric(NA)),
    as.numeric(NA),
    info = "numeric NA --> unchanged in output"
  )

  expect_equal(
    inf_omit(NaN),
    NaN,
    info = "NaN input --> unchanged in output"
  )

  expect_equal(
    inf_omit(c("abc", "def", "Inf", "-Inf", "NA", NA)),
    c("abc", "def", "Inf", "-Inf", "NA", NA),
    info = "Inf omit shouldn't affect strings"
  )
})



###############################################################################
test_that("na_inf_omit", {
  input_null <- NULL
  input_empty_number <- numeric(0)
  input_single_inf <- c(Inf)
  input_number_and_inf <- c(1, 2, Inf, 3, -Inf, 1)
  input_na <- as.numeric(NA)
  input_nan <- NaN
  input_string <- c("abc", "def", "Inf", "-Inf", "NA", NA)

  expect_equal(
    na_inf_omit(input_null),
    NULL,
    info = "na_inf_omit: NULL input"
  )

  expect_equal(
    na_inf_omit(input_empty_number),
    numeric(0),
    info = "na_inf_omit: Empty input"
  )

  expect_equal(
    na_inf_omit(input_single_inf),
    numeric(0),
    info = "na_inf_omit: Single Inf as input"
  )

  expect_equal(
    na_inf_omit(input_number_and_inf),
    c(1, 2, 3, 1),
    info = "na_inf_omit: mix of numbers and infs"
  )

  expect_equal(
    na_inf_omit(input_na),
    numeric(0),
    info = "na_inf_omit: numeric NA input"
  )

  expect_equal(
    na_inf_omit(input_nan),
    numeric(0),
    info = "na_inf_omit: NaN input"
  )
  expect_equal(
    na_inf_omit(input_string),
    c("abc", "def", "Inf", "-Inf", "NA"),
    info = "na_inf_omit: collection of strings"
  )
})
