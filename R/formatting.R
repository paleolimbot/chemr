
#' Format objects as markdown
#'
#' @param x A mol, reaction, or reaction list object
#' @param equals_sign An equals sign to use for a reaction object
#' @param ... Passed to as.character()
#'
#' @return A character vector formatted for markdown
#' @export
#'
#' @examples
#' chem_markdown(as_reaction("2H+ + 2OH- = 2H2O"))
#' chem_markdown(as_reaction_list("2H+ + 2OH- = 2H2O"))
#' chem_markdown(as_mol(c("SO4-2", "Ca+2")))
#'
chem_markdown <- function(x, ...) {
  UseMethod("chem_markdown")
}

#' @rdname chem_markdown
#' @export
chem_markdown.default <- function(x, ...) {
  as.character(
    x,
    wrap_super = function(x) paste0("^", x, "^"),
    wrap_sub = function(x) paste0("~", x, "~"),
    ...
  )
}

#' @rdname chem_markdown
#' @export
chem_markdown.reaction <- function(x, equals_sign = NULL, ...) {
  if(is.null(equals_sign)) {
    equals_sign <- "\u21CC"
  }

  as.character(
    x,
    equals_sign = equals_sign,
    wrap_super = function(x) paste0("^", x, "^"),
    wrap_sub = function(x) paste0("~", x, "~"),
    wrap_coeff = function(x) paste0("**", x, "**"),
    ...
  )
}

#' @rdname chem_markdown
#' @export
chem_markdown.reaction_list <- function(x, ...) {
  chem_markdown.reaction(x, ...)
}


#' Format objects as plotmath expressions
#'
#' For use on plots using parse = TRUE
#'
#' @inheritParams chem_markdown
#' @export
#'
#' @examples
#' chem_expression(as_reaction("2H+ + 2OH- = 2H2O"))
#' chem_expression(as_reaction_list("2H+ + 2OH- = 2H2O"))
#' chem_expression(as_mol(c("SO4-2", "Ca+2")))
#'
chem_expression <- function(x, ...) {
  UseMethod("chem_expression")
}

#' @rdname chem_expression
#' @export
chem_expression.default <- function(x, ...) {
  paste0(
    "paste(",
    as.character(
      x,
      element_sep = ", ",
      wrap_super = function(x) paste0("{}^'", x, "'"),
      wrap_sub = function(x) paste0("[", x, "]"),
      wrap_sub_molecule = function(x) paste0("paste('(', ", x, ", ')')"),
      ...
    ),
    ")"
  )
}

#' @rdname chem_expression
#' @export
chem_expression.reaction <- function(x, equals_sign = NULL, ...) {
  if(is.null(equals_sign)) {
    equals_sign <- ", ' = ', "
  }

  paste0(
    "paste(",
    as.character(
      x,
      equals_sign = equals_sign,
      mol_sep = ", ' + ', ",
      element_sep = ", ",
      wrap_super = function(x) paste0("{}^'", x, "'"),
      wrap_sub = function(x) paste0("[", x, "]"),
      wrap_coeff = function(x) ifelse(x == "", "", paste0("'", x, "', ")),
      ...
    ),
    ")"
  )
}

#' @rdname chem_expression
#' @export
chem_expression.reaction_list <- function(x, ...) {
  chem_expression.reaction(x, ...)
}
