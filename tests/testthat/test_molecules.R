
context("molecule and molecules objects")


test_that("molecule character parsing works as intended", {
  # good
  expect_silent(as_molecule("H2O"))
  expect_equivalent(as_molecule("H2O")[], c(H=2, O=1))
  expect_equal(charge(as_molecule("H2O")), 0)
  expect_silent(as_molecule("H3O+"))
  expect_silent(as_molecule("NH3-"))
  expect_silent(as_molecule("Mg+2"))
  expect_silent(as_molecule("Cl-"))
  expect_silent(as_molecule("O-2"))

  # bad
  expect_warning(as_molecule("H2o"), "Bad molecule text:.*")
  expect_error(as_molecule("Hh2O"), "names\\(x\\) contained the following bad symbols:.*")
  expect_warning(as_molecule("H2O++"), "Bad molecule text:.*")
  expect_warning(as_molecule("Hhhh2O+"), "Bad molecule text:.*")
})
