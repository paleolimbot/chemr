
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

