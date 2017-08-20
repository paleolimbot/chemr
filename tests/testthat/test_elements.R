
context("elements and the periodic tibble")

test_that("element functions return the desired values", {

  expect_true(is_element("H"))
  expect_false(is_element("H2O"))

  expect_identical(elsymbol(6), "C")
  expect_identical(elz("C"), c(C=6L))
  expect_is(elmass("C"), "numeric")

  expect_identical(elname("H"), c(H="Hydrogen"))
  expect_identical(elgroup("H"), c(H=1L))
  expect_identical(elperiod("H"), c(H=1L))

  expect_identical(elvalence("H")$H, c(-1L, 0L, 1L))

  # NA handling
  expect_identical(elsymbol(NA), NA_character_)
  expect_true(is.na(elz(NA)))
  expect_true(is.na(elmass(NA)))
  expect_true(is.na(elname(NA)))
  expect_true(is.na(elperiod(NA)))
  expect_true(is.na(elgroup(NA)))
  expect_is(elvalence(NA), "list")
  expect_null(elvalence(NA)[[1]])

  # invalid value handling
  expect_identical(elsymbol(-1), NA_character_)
  expect_true(is.na(elz("NA")))
  expect_true(is.na(elmass("NA")))
  expect_true(is.na(elname("NA")))
  expect_true(is.na(elperiod("NA")))
  expect_true(is.na(elgroup("NA")))
  expect_is(elvalence("NA"), "list")
  expect_null(elvalence("NA")[[1]])

  # vectorization
  expect_length(elsymbol(1:4), 4)
  expect_length(elsymbol(-10:10), 21)
  expect_length(elz(c("H", "He", "Li")), 3)
  expect_length(elmass(c("H", "He", "Li")), 3)
  expect_length(elname(c("H", "He", "Li")), 3)
  expect_length(elperiod(c("H", "He", "Li")), 3)
  expect_length(elgroup(c("H", "He", "Li")), 3)
  expect_length(elvalence(c("H", "He", "Li")), 3)
})

test_that("elements as lists work as expected", {

  expect_is(with(elzs, H:O), "integer")
  expect_length(with(elzs, H:O), 8)
  expect_equal(with(elmasses, 2*H + O), mass("H2O"))
  expect_identical(with(elvalences, H), c(-1L, 0L, 1L))

})
