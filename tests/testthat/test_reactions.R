
context("reaction objects")

test_that("character parsing of reactions works as intended", {
  robj <- as_reaction("2H+ + 2OH- = 2H2O")
  expect_equal(as.character(robj$mol), c("H+", "OH-", "H2O"))
  expect_equal(robj$coefficient, c(2, 2, -2))
})

test_that("formula parsing works as intended", {
  robj <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  expect_equal(as.character(robj$mol), c("O2", "H2", "H2O"))
  expect_equal(robj$coefficient, c(1, 2, -2))
})

test_that("reaction constructor works as intended", {
  robj <- reaction(lhs = c("H2", "O2"), rhs = "H2O",
                   counts_lhs = c(2, 1), counts_rhs = 2)
  expect_equal(as.character(robj$mol), c("H2", "O2", "H2O"))
  expect_equal(robj$coefficient, c(2, 1, -2))

  robj2 <- as_reaction(as_mol(c("H2", "O2", "H2O")),
                       coefficient = c(2, 1, -2))
  expect_identical(robj, robj2)

})

test_that("reaction subsetting works as intended", {
  robj <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  expect_is(lhs(robj), "reaction")
  expect_equal(as.character(lhs(robj)$mol), c("O2", "H2"))
  expect_equal(lhs(robj)$coefficient, c(1, 2))
  expect_is(rhs(robj), "reaction")
  expect_equal(as.character(rhs(robj)$mol), "H2O")
  expect_equal(rhs(robj)$coefficient, -2)
})

test_that("reaction charater printing works as intended", {
  rtxt <- "O2 + 2H2 = 2H2O"
  robj <- as_reaction(rtxt)
  expect_identical(as.character(robj), rtxt)
  expect_output(print(robj), "<reaction>.*")
  expect_identical(print(robj), robj)
})

test_that("reaction simplification works as intended", {
  r <- as_reaction("NH3 + H+ + H2O = H3O+ + NH3")
  expect_equal(simplify(r)$mol, unique(r$mol))
  expect_is(simplify(r), "reaction")

  r2 <- as_reaction("0NH3 + H+ + H2O = H3O+")
  expect_is(remove_zero_counts(r2), "reaction")
  expect_equal(as.character(remove_zero_counts(r2)$mol), c("H+", "H2O", "H3O+"))

  r3 <- as_reaction(as_mol(c("H2O", "H+", NA_character_)),
                           coefficient = c(1, 2, -3))
  expect_is(remove_zero_counts(r3), "reaction")
  expect_equal(as.character(remove_zero_counts(r3)$mol), c("H2O", "H+"))
})

test_that("is_balanced() checks for a balanced reaction", {
  r1 <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  expect_true(is_balanced(r1))
  r2 <- as_reaction(O2 + H2 ~ 2*H2O)
  expect_false(is_balanced(r2))
  r3 <- as_reaction(`O2-4` + 2*H2 ~ 2*H2O)
  expect_true(is_balanced(r3, charge = FALSE))
  expect_false(is_balanced(r3, charge = TRUE))
})

test_that("data frame representations are correct", {
  r1 <- as_reaction(O2 + 2*H2 ~ 2*H2O)

  # generate representations
  mat <- as.matrix(r1)
  tbl <- tibble::as_tibble(r1)
  df <- as.data.frame(r1)

  # matrix output
  expect_true(is_balanced(r1))
  expect_is(mat, "matrix")
  expect_equal(sum(mat), 0)
  expect_equal(ncol(mat), 2)
  expect_equal(nrow(mat), 3)

  # tbl output
  expect_is(tbl, "tbl")
  expect_is(tbl, "data.frame")
  expect_identical(as.data.frame(tbl), df)
  expect_equal(nrow(tbl), 3)
  expect_equal(names(tbl), c("mol", "mol_text", "charge", "mass", "coefficient",
                             "O", "H"))
  expect_equal(tbl$coefficient, r1$coefficient)

  # empty intput
  r0 <- reaction(mol(), mol())
  expect_is(as.matrix(r0), "matrix")
  expect_equal(nrow(as.matrix(r0)), 0)
  expect_equal(ncol(as.matrix(r0)), 0)
  expect_is(tibble::as_tibble(r0), "tbl")
  expect_equal(nrow(tibble::as_tibble(r0)), 0)
  expect_equal(ncol(tibble::as_tibble(r0)), 5)
  expect_is(as.data.frame(r0), "data.frame")
  expect_identical(as.data.frame(r0), as.data.frame(tibble::as_tibble(r0)))
})

test_that("balance() balances reaction", {
  r1 <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  r2 <- as_reaction(O2 + H2 ~ H2O)
  r3 <- as_reaction(`O2-4` + 2*H2 ~ 2*H2O)

  # r1 is balanced
  expect_identical(balance(r1), r1)
  # r2 is not
  expect_true(is_balanced(balance(r2)))
  expect_identical(balance(r2), r1)

  # check geochem examples
  # NaAlSi3O8 <=> KAl3Si3O10(OH)2
  r4 <- as_reaction("NaAlSi3O8 = KAl3Si3O10O2H2 + SiO2 + K+ + Na+ + H+")
  expect_true(is_balanced(balance(r4)))
  # K1.5Al5.5Si6.5O20(OH)4 <=> Al2Si2O5(OH)4 + SiO2
  r5 <- as_reaction("K1.5Al5.5Si6.5O20O4H4 = Al2Si2O5O4H4 + SiO2 + K+ + H2O + H+")
  expect_true(is_balanced(balance(r5)))
  # Fe+2 = Fe(OH)3
  r6 <- as_reaction("Fe+2 = FeO3H3 + H2O + H+")
  expect_true(is_balanced(balance(r6, charge = F), charge = F))
  # H2O = O2 + H+
  r7 <- as_reaction("H2O = O2 + H+")
  expect_true(is_balanced(balance(r7, charge = F), charge = F))
  # Fe+2 + O2 = Fe(OH)3 (doesn't work!)
  #r8_balanced <- as_reaction("4Fe+2 + O2 + 10H2O = 4FeO3H3 + 8H+")
  #r8 <- as_reaction("Fe+2 + O2 + H2O = FeO3H3 + H+")
  #is_balanced(r8_balanced)
  #balance(r8)

  # Clinochlore [Mg10Al4Si6O20(OH)16] + Muscovite + Calcite =>
  # Clinozoisite [Ca2Al3Si3O12(OH)] + Phlogopite

  # Clinochlore [Mg10Al4Si6O20(OH)16] + Muscovite =>
  # Pyrope [Mg3Al2Si3O12] + K-Spar

  # Phlogopite =>
  # Montmorillonite [Na2/3Mg2/3Al10/3Si8O20(OH)4]

  # Anorthite + Diopside =>
  # Clinochlore [Mg10Al4Si6O20(OH)16] + Calcite

  # Diopside + Hornblende [NaCa2Al3Mg4Si6O22(OH)2] =>
  # Mg-Chlorite [Mg10Al4Si6O20(OH)16] + Calcite + Albite

  # Clinochlore [Mg10Al4Si6O20(OH)16] + Calcite =>
  # Kaolinite + Dolomite

  # Albite => Illite using Si, Al

  # Magnetite + H2S => Pyrite

  # Mg-Chlorite [Mg7Fe3Al4Si6O20(OH)16] + Calcite =>
  # Muscovite + Pyrite + Dolomite
  # on Al, Fe

  # Ilmenite => Rutile + Hematite

  # KAl3Si3O10(OH)2 = Mg10Al4Si6O20(OH)16
  # using SiO2, H2O, K+, H+, Mg+2

  # final:
  # 4 KAl3Si3O10(OH)2 + 6 SiO2 + 48 H2O + 30 Mg+2 =
  # 3 Mg10Al4Si6O20(OH)16 + 56 H+ + 4 K+

})
