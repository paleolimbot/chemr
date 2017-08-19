
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

  # make sure is vectorized
  expect_length(as_molecules(c("H2O", "H+", "H3O+")), 3)
  expect_is(as_molecules(c("H2O", "H+", "H3O+")), "molecules")
  expect_is(as_molecules("H2O"), "molecules")
  expect_is(as_molecule("H2O"), "molecule")

  # NA handling
  expect_identical(as_molecule(NA_character_), NA_molecule_)
  expect_identical(as_molecule("NA"), NA_molecule_)
  expect_identical(as_molecule("<NA_molecule_>"), NA_molecule_)

  # bad
  expect_warning(as_molecule("H2o"), "Bad molecule text:.*")
  expect_error(as_molecule("Hh2O"), "names\\(x\\) contained the following bad symbols:.*")
  expect_warning(as_molecule("H2O++"), "Bad molecule text:.*")
  expect_warning(as_molecule("Hhhh2O+"), "Bad molecule text:.*")
})

test_that("formula works for creating molecule(s)", {
  expect_identical(as_molecule(~H2O), as_molecule('H2O'))
  expect_identical(as_molecules(~H2O), as_molecules('H2O'))
  expect_identical(as_molecules(~H2O + H2), as_molecules(c('H2O', "H2")))
})

test_that("coersion rules work as expected", {
  water <- as_molecule(~H2O)
  waters <- molecules(water)
  expect_false(identical(water, waters))
  expect_identical(waters[[1]], water)
  expect_identical(as_molecules(water), waters)
  expect_identical(as_molecule(waters), water)
  expect_warning(as_molecule(molecules(water, water)), "Using first of 2 molecules in x")
})

test_that("molecule(s) constructors work as expected", {
  expect_identical(molecule(H=1), as_molecule("H"))
  expect_identical(molecule(H=1, charge = 1), as_molecule("H+"))
  expect_identical(molecule(H=NA), molecule(H=1))
})

test_that("the print method works for molecule and molecules", {
  water <- as_molecule(~H2O)
  waters <- molecules(water, water)
  expect_output(print(water), "<molecule>.*")
  expect_output(print(waters, "<molecules>.*?[1].*"))
  expect_identical(water, print(water))
  expect_identical(waters, print(waters))
})

test_that("converting molecules to strings works as expected", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", "<NA_molecule>")
  expect_identical(as.character(mols), mols)
})
