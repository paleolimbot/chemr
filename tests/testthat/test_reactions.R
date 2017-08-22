
context("reaction objects")

test_that("character parsing of reactions works as intended", {
  robj <- as_reaction("2H+ + 2OH- = 2H2O")
  expect_equal(as.character(robj$mol), c("H+", "OH-", "H2O"))
  expect_equal(robj$coefficients, c(2, 2, -2))
})

test_that("formula parsing works as intended", {
  robj <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  expect_equal(as.character(robj$mol), c("O2", "H2", "H2O"))
  expect_equal(robj$coefficients, c(1, 2, -2))
})

test_that("reaction constructor works as intended", {
  robj <- reaction(lhs = c("H2", "O2"), rhs = "H2O",
                   counts_lhs = c(2, 1), counts_rhs = 2)
  expect_equal(as.character(robj$mol), c("H2", "O2", "H2O"))
  expect_equal(robj$coefficients, c(2, 1, -2))
})
