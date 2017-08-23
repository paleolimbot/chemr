
context("molecule and mol objects")


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
  expect_length(as_mol(c("H2O", "H+", "H3O+")), 3)
  expect_is(as_mol(c("H2O", "H+", "H3O+")), "mol")
  expect_is(as_mol("H2O"), "mol")
  expect_is(as_molecule_single("H2O"), "molecule_single")

  # NA handling
  expect_identical(as_molecule_single(NA_character_), NA_molecule_)
  expect_identical(as_molecule_single("NA"), NA_molecule_)
  expect_identical(as_molecule_single(NA_character_), NA_molecule_)
  expect_silent(as_mol(c("H2O", NA_character_)))

  # zero count handling
  expect_identical(as_molecule_single("H2O0"), molecule_single(H=2))

  # bad
  expect_warning(as_molecule_single("H2o"), "Bad molecule text:.*")
  expect_error(as_molecule_single("Hh2O"), "names\\(x\\) contained the following bad symbols:.*")
  expect_warning(as_molecule_single("H2O++"), "Bad molecule text:.*")
  expect_warning(as_molecule_single("Hhhh2O+"), "Bad molecule text:.*")
})

test_that("formula works for creating molecule(s)", {
  expect_identical(as_molecule_single(~H2O), as_molecule_single('H2O'))
  expect_identical(as_mol(~H2O), as_mol('H2O'))
  expect_identical(as_mol(~H2O + H2), as_mol(c('H2O', "H2")))
})

test_that("coersion rules work as expected", {
  water <- as_molecule_single(~H2O)
  waters <- mol(water)
  expect_false(identical(water, waters))
  expect_identical(waters[[1]], water)
  expect_identical(as_mol(water), waters)
  expect_identical(as_molecule_single(waters), water)
  expect_warning(as_molecule_single(mol(water, water)), "Using first of 2 mol in x")
})

test_that("molecule(s) constructors work as expected", {
  expect_identical(molecule_single(H=1), as_molecule_single("H"))
  expect_identical(molecule_single(H=1, charge = 1), as_molecule_single("H+"))
  expect_identical(molecule_single(H=NA), molecule_single(H=1)) # NA <- 1
  expect_identical(molecule_single(H=0, O=2), molecule_single(O = 2)) # 0 counts removed
})

test_that("the print method works for molecule and mol", {
  water <- as_molecule_single(~H2O)
  waters <- mol(water, water)
  expect_output(print(water), "<molecule_single>.*")
  expect_output(print(waters, "<mol>.*?[1].*"))
  expect_identical(water, print(water))
  expect_identical(waters, print(waters))
})

test_that("as.mol and as_mol are identical", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  expect_identical(as.mol(mols), as_mol(mols))
})

test_that("converting mol to strings works as expected", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_objs <- as_mol(mols)
  expect_identical(as.character(mol_objs), mols)
})

test_that("subsetting a mol object returns a mol object", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_objs <- as_mol(mols)
  expect_is(mol_objs[1], "mol")
  expect_is(mol_objs[1:4], "mol")
  expect_is(mol_objs[mols == "H2O"], "mol")
  expect_true(is.na(mol_objs[NA_integer_]))
})

test_that("combining molecule(s) returns a mol object", {
  mols1 <- c("H2O", "H3O+", "H+", "Mg+2")
  mols2 <- c("O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_objs1 <- as_mol(mols1)
  mol_objs2 <- as_mol(mols2)
  silica <- as_molecule_single(~SiO2)
  water <- as_molecule_single(~H2O)

  expect_is(c(silica, water), "mol")
  expect_equal(as.character(c(silica, water)), c("SiO2", "H2O"))
  expect_is(c(silica, mols1), "mol")
  expect_equal(as.character(c(silica, mol_objs1)), c("SiO2", mols1))
  expect_is(c(mol_objs1, mol_objs2), "mol")
  expect_equal(as.character(c(mol_objs1, mol_objs2)), c(mols1, mols2))
})

test_that("rep works for molecule objects", {
  m1 <- as_molecule_single(~H2O)
  mols <- as_mol(c("O-2", "O2", "CH3COOH", "Cl-", NA_character_))

  # check class
  expect_is(rep(m1, 3), "mol")
  expect_is(rep(mols, 3), "mol")

  # check length
  expect_length(rep(m1, 3), 3)
  expect_length(rep(mols, 3), length(mols) * 3)

  # check values
  expect_equal(as.character(rep(m1, 3)), c("H2O", "H2O", "H2O"))
  expect_equal(as.character(rep(mols, 3)),
               rep(c("O-2", "O2", "CH3COOH", "Cl-", NA_character_), 3))
})

test_that("unique() works for mol vectors", {
  mols <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_objs <- as_mol(mols)
  expect_is(unique(mol_objs), "mol")
  expect_identical(unique(rep(mol_objs, 2)), mol_objs)
})

test_that("is.na() works for molecule_single and mol objects", {
  expect_true(is.na(NA_molecule_))
  expect_false(is.na(as_molecule_single(~H2O)))
  mols <- mol(~H2O, ~`H+`, NA_molecule_, ~NH2)
  expect_identical(is.na(mols), c(FALSE, FALSE, TRUE, FALSE))
})

test_that("molecule_single arithmetic works as intended", {
  m1 <- as_molecule_single(~H2O)
  m2 <- as_molecule_single(~`H+`)

  # multiplication, division
  expect_is(m1 * 2, "molecule_single")
  expect_equal(as.character(m1*2), "H4O2")
  expect_identical(2 * m1, m1 * 2) # communicativity
  expect_true(is.na(NA * m1)) # NA handling
  expect_true(is.na(NA_molecule_ * 4))
  # default error for non-numerics is fine
  expect_error(m1 * "fish", "non-numeric argument to binary operator")
  expect_is(m1 / 2, "molecule_single")
  expect_equal(as.character(m1 / 2), "HO0.5")
  expect_error(2 / m1, "Can't divide by a molecule_single") # can't divide by a mol
  expect_error(m1 / "fish", "Can't divide a molecule_single by an object of type character")
  expect_true(is.na(m1 / NA)) # NA handling
  expect_true(is.na(NA_molecule_ / 4))
  # divide by zero
  expect_identical(m1 / 0, NA_molecule_)
  expect_warning(m1 / 0, "Divide by zero")

  # order matters in addition
  expect_is(m1 + m2, "molecule_single")
  expect_equal(as.character(m1 + m2), "H2OH+")
  expect_equal(as.character(m2 + m1), "HH2O+")
  # type coersion
  expect_identical(m1 + "H+", m1 + m2)
  expect_identical("H+" + m1, m2 + m1)
  # NA handling
  expect_true(is.na(m1 + NA_molecule_))
  expect_true(is.na(NA_molecule_ + m1))
  expect_true(is.na(NA_molecule_ + NA_molecule_))
  # unary operator
  expect_identical(+m1, m1)



  # equality operator
  expect_length(m1 == m2, 1)
  expect_false(m1 == m2)
  expect_true(m1 == m1)
  expect_true(m2 == m2)
  # type coersion
  expect_true(m1 == "H2O")
  expect_true("H2O" == m1)
  expect_false(m1 == "H+")

  # make sure charge is considered
  m3 <- as_molecule_single(~H)
  expect_false(m2 == m3)

  # NA handling
  expect_true(is.na(m1 == NA_molecule_))
  expect_true(is.na(NA_molecule_ == m1))
  expect_true(is.na(NA_molecule_ == NA_molecule_))
})

test_that("mol arithmetic works as intended", {
  mol_text <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_single <- as_molecule_single(~H2O)
  mol_len_1 <- as_mol(~H2O)
  mols <- as_mol(mol_text)

  # check dimensions for *, /
  expect_is(mol_len_1 * 2, "mol")
  expect_length(mol_len_1, 1)
  expect_is(mols * 2, "mol")
  expect_length(mols * 2, length(mols))

  # check dimensions for +
  expect_length(mols + mol_len_1, length(mols))
  expect_equal(which(is.na(mols)), which(is.na(mols + mol_len_1))) # NA propogation
  expect_length(mols + mols, length(mols))
  # communicativity
  expect_equal(length(mol_len_1 + mols), length(mols + mol_len_1))
  # unary operator (doesn't work for mols!)
  # expect_identical(+mols, mols)

  # mols + mol_single # doesn't work, citing 'incompatible methods'

  # check ==
  expect_is(mols == mols, "logical")
  expect_length(mols == mols, length(mols))
  expect_true(all(na.omit(mols == mols)))
  expect_equal(which(is.na(mols)), which(is.na(mols == mols)))
})

test_that("combine_molecules works for both molecule_single and mol objects", {
  mol_text <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", "Cl-", NA_character_)
  mol_single <- as_molecule_single(~H2O)
  mol_len_1 <- as_mol(~H2O)
  mols <- as_mol(mol_text)

  # test with mol_single objects
  expect_is(combine_molecules(mol_single, mol_single), "molecule_single")
  expect_identical(combine_molecules(mol_single, mol_single), mol_single + mol_single)

  # test with mol objects
  expect_is(combine_molecules(mols, mols), "mol")
  expect_identical(combine_molecules(mols, mols), mols + mols)

  # test vectorization
  expect_is(combine_molecules(mols, mol_single), "mol")
  expect_length(combine_molecules(mols, mol_single), length(mols))
  expect_identical(combine_molecules(mols, mol_single),
                   combine_molecules(mols, mol_len_1))
  expect_identical(combine_molecules(mols, mol_len_1),
                   mols + mol_len_1)

  # test NA handling
  expect_true(is.na(combine_molecules(mol_single, NA_molecule_)))
  expect_true(is.na(combine_molecules(NA_molecule_, NA_molecule_)))
})

test_that("remove zero values, simplify works as intended", {

  m1 <- as_molecule_single("CH3COOH")
  m2 <- as_molecule_single("C2H4O2")
  expect_identical(simplify(m1), m2)

  expect_identical(remove_zero_counts(new_molecule_single(c(H=2, O=0), mass = unname(2*elmass("H")))),
                   as_molecule_single("H2"))

  # NA handling
  expect_identical(simplify(NA_molecule_), NA_molecule_)
  expect_identical(remove_zero_counts(NA_molecule_), NA_molecule_)
})

test_that("data frame, matrix representations are correct", {
  mol_text <- c("H2O", "H3O+", "H+", "Mg+2", "O-2", "O2", "CH3COOH", NA_character_,
                "Cl-", NA_character_)
  mol_single <- as_molecule_single(~H2O)
  mol_len_1 <- as_mol(~H2O)
  mols <- as_mol(mol_text)

  # create all
  mat <- as.matrix(mols)
  tbl <- tibble::as_tibble(mols)
  df <- as.data.frame(mols)

  # check class
  expect_is(mat, "matrix")
  expect_is(tbl, "tbl")
  expect_is(df, "data.frame")

  # check dimensions
  expect_equal(length(mols), nrow(mat))
  expect_equal(length(mols), nrow(tbl))
  expect_equal(length(mols), nrow(df))
  expect_equal(ncol(mat), 5)
  required_names <- c("mol", "mol_text", "mass", "charge")
  expect_true(all(required_names %in% names(tbl)))
  expect_true(all(required_names %in% names(df)))
  expect_equal(length(required_names) + ncol(mat), ncol(tbl))
  expect_equal(length(required_names) + ncol(mat), ncol(df))

  # recreate mols from df
  mat2 <- df[-(1:3)] %>%
    plyr::mlply(molecule_single) %>%
    new_mol() %>%
    as.matrix()
  # equal except for rownames
  rownames(mat2) <- NULL
  rownames(mat) <- NULL
  expect_equal(mat, mat2)

  # one length input
  expect_identical(as.matrix(mol_single), as.matrix(mol_len_1))
  expect_identical(as.data.frame(mol_single), as.data.frame(mol_len_1))
  expect_identical(tibble::as_tibble(mol_single), tibble::as_tibble(mol_len_1))

  # zero-length input
  expect_is(as.matrix(mol()), "matrix")
  expect_equal(nrow(as.matrix(mol())), 0)
  expect_equal(ncol(as.matrix(mol())), 0)
  expect_is(as.data.frame(mol()), "data.frame")
  expect_equal(nrow(as.data.frame(mol())), 0)
  expect_equal(ncol(as.data.frame(mol())), 4)
  expect_is(tibble::as_tibble(mol()), "data.frame")
  expect_equal(nrow(tibble::as_tibble(mol())), 0)
  expect_equal(ncol(tibble::as_tibble(mol())), 4)
})
