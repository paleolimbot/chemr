
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
#' as_markdown(as_reaction("2H+ + 2OH- = 2H2O"))
#' as_markdown(as_reaction_list("2H+ + 2OH- = 2H2O"))
#' as_markdown(as_mol(c("SO4-2", "Ca+2")))
#'
as_markdown <- function(x, ...) {
  UseMethod("as_markdown")
}

#' @rdname as_markdown
#' @export
as_markdown.default <- function(x, ...) {
  as.character(
    x,
    wrap_super = function(x) paste0("^", x, "^"),
    wrap_sub = function(x) paste0("~", x, "~"),
    ...
  )
}

#' @rdname as_markdown
#' @export
as_markdown.reaction <- function(x, equals_sign = NULL, ...) {
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

#' @rdname as_markdown
#' @export
as_markdown.reaction_list <- function(x, ...) {
  as_markdown.reaction(x, ...)
}
