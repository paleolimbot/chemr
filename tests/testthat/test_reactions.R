
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

  robj2 <- as_reaction(as_mol(c("H2", "O2", "H2O")),
                       coefficients = c(2, 1, -2))
  expect_identical(robj, robj2)

})

test_that("reaction subsetting works as intended", {
  robj <- as_reaction(O2 + 2*H2 ~ 2*H2O)
  expect_is(lhs(robj), "reaction")
  expect_equal(as.character(lhs(robj)$mol), c("O2", "H2"))
  expect_equal(lhs(robj)$coefficients, c(1, 2))
  expect_is(rhs(robj), "reaction")
  expect_equal(as.character(rhs(robj)$mol), "H2O")
  expect_equal(rhs(robj)$coefficients, -2)
})

test_that("reaction charater printing works as intended", {
  rtxt <- "O2 + 2H2 = 2H2O"
  robj <- as_reaction(rtxt)
  expect_identical(as.character(robj), rtxt)
  expect_output(print(robj), "<reaction>.*")
  expect_identical(print(robj), robj)
})
