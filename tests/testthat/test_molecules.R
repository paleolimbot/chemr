
context("molecule and molecules objects")


test_that("molecule character parsing works as intended", {
  # good
  expect_silent(as_molecule_single("H2O"))
  expect_is(as_molecule_single("H2O"), "molecule_single")
  expect_equivalent(as_molecule_single("H2O")[], c(H=2, O=1))
  expect_equal(charge(as_molecule_single("H2O")), 0)
  expect_silent(as_molecule_single("H3O+"))
  expect_silent(as_molecule_single("NH3-"))
  expect_silent(as_molecule_single("Mg+2"))
  expect_silent(as_molecule_single("Cl-"))
  expect_silent(as_molecule_single("O-2"))

  # make sure is vectorized
  expect_length(as_molecules(c("H2O", "H+", "H3O+")), 3)
  expect_is(as_molecules(c("H2O", "H+", "H3O+")), "molecules")
  expect_is(as_molecules("H2O"), "molecules")
  expect_is(as_molecule_single("H2O"), "molecule_single")

  # NA handling
  expect_identical(as_molecule_single(NA_character_), NA_molecule_)
  expect_identical(as_molecule_single("NA"), NA_molecule_)
  expect_identical(as_molecule_single("<NA_molecule_>"), NA_molecule_)
  expect_silent(as_molecules(c("H2O", "<NA_molecule_>")))

  # bad
  expect_warning(as_molecule_single("H2o"), "Bad molecule text:.*")
  expect_error(as_molecule_single("Hh2O"), "names\\(x\\) contained the following bad symbols:.*")
  expect_warning(as_molecule_single("H2O++"), "Bad molecule text:.*")
  expect_warning(as_molecule_single("Hhhh2O+"), "Bad molecule text:.*")
})

test_that("formula works for creating molecule(s)", {
  expect_identical(as_molecule_single(~H2O), as_molecule_single('H2O'))
  expect_identical(as_molecules(~H2O), as_molecules('H2O'))
  expect_identical(as_molecules(~H2O + H2), as_molecules(c('H2O', "H2")))
})

test_that("coersion rules work as expected", {
  water <- as_molecule_single(~H2O)
  waters <- molecules(water)
  expect_false(identical(water, waters))
  expect_identical(waters[[1]], water)
  expect_identical(as_molecules(water), waters)
  expect_identical(as_molecule_single(waters), water)
  expect_warning(as_molecule_single(molecules(water, water)), "Using first of 2 molecules in x")
})

test_that("molecule(s) constructors work as expected", {
  expect_identical(molecule_single(H=1), as_molecule_single("H"))
  expect_identical(molecule_single(H=1, charge = 1), as_molecule_single("H+"))
  expect_identical(molecule_single(H=NA), molecule_single(H=1))
})

test_that("the print method works for molecule and molecules", {
  water <- as_molecule_single(~H2O)
  waters <- molecules(water, water)
  expect_output(print(water), "<molecule_single>.*")
  expect_output(print(waters, "<molecules>.*?[1].*"))
  expect_identical(water, print(water))
  expect_identical(waters, print(waters))
})

test_that("converting molecules to strings works as expected", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", "<NA_molecule_>")
  expect_identical(as.character(mols), mols)
})

test_that("subsetting a molecules object returns a molecules object", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", "<NA_molecule_>")
  mol_objs <- as_molecules(mols)
  expect_is(mol_objs[1], "molecules")
  expect_is(mol_objs[1:4], "molecules")
  expect_is(mol_objs[mols == "H2O"], "molecules")
})

test_that("combining molecule(s) returns a molecules object", {
  mols1 <- c("H2O", "H3O+", "H+", "Mg+2")
  mols2 <- c("O-2", "O2", "CH3COOH", "Cl-", "<NA_molecule_>")
  mol_objs1 <- as_molecules(mols1)
  mol_objs2 <- as_molecules(mols2)
  silica <- as_molecule_single(~SiO2)
  water <- as_molecule_single(~H2O)

  expect_is(c(silica, water), "molecules")
  expect_equal(as.character(c(silica, water)), c("SiO2", "H2O"))
  expect_is(c(silica, mols1), "molecules")
  expect_equal(as.character(c(silica, mol_objs1)), c("SiO2", mols1))
  expect_is(c(mol_objs1, mol_objs2), "molecules")
  expect_equal(as.character(c(mol_objs1, mol_objs2)), c(mols1, mols2))
})
