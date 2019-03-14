context("test-formatting")

test_that("formatting works", {
  chem_exprs <- c(
    chem_expression(as_reaction("2H+ + 2OH- = 2H2O")),
    chem_expression(as_reaction_list("2H+ + 2OH- = 2H2O")),
    chem_expression(as_mol(c("SO4-2", "Ca+2", "H2O", "Pb(OH)2", "Pb(SO4)2-2")))
  )

  expect_is(parse(text = chem_exprs), "expression")

  plot(seq_along(chem_exprs), type = "n")
  text(x = 1, y = seq_along(chem_exprs), labels = parse(text = chem_exprs), adj = c(0, 0.5))
})
