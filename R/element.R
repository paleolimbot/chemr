
#' Element properties
#'
#' Based on the periodic tibble (see \link{pt}).
#'
#' @param x An element symbol (e.g. H, He)
#' @param z The atomic number Z
#' @param locale Currently not supported, but will eventually return "Al" as "Aluminum"
#'   for unsavory Americans.
#'
#' @return A named vector of the relevant value(s)
#' @export
#' @rdname element
#'
#' @examples
#' is_element("H")
#' is_element("Hg")
#' is_element("H2O")
#'
#' elsymbol(1)
#'
#' elmass("H")
#' elmass("O")
#'
#' elz("Pb")
#' elname("Al")
#' elgroup("Al")
#' elperiod("Al")
#' elvalence("Al")
#'
is_element <- function(x) {
  x %in% as.character(pt$symbol)
}

#' @rdname element
#' @export
elsymbol <- function(z) {
  z <- ifelse(is.na(z) | (z < 1) | (z > nrow(pt)), NA_integer_, z)
  pt$symbol[z]
}

#' @rdname element
#' @export
elmass <- function(x) {
  stats::setNames(pt$mass[match(x, pt$symbol)], x)
}

#' @rdname element
#' @export
elz <- function(x) {
  stats::setNames(pt$z[match(x, pt$symbol)], x)
}

#' @rdname element
#' @export
elname <- function(x, locale = "en_UK") {
  stats::setNames(pt$name[match(x, pt$symbol)], x)
}

#' @rdname element
#' @export
elgroup <- function(x) {
  stats::setNames(pt$group[match(x, pt$symbol)], x)
}

#' @rdname element
#' @export
elperiod <- function(x) {
  stats::setNames(pt$period[match(x, pt$symbol)], x)
}

#' @rdname element
#' @export
elvalence <- function(x) {
  stats::setNames(pt$valence[match(x, pt$symbol)], x)
}

#' Element properties as lists
#'
#' Based on the periodic tibble (see \link{pt}).
#'
#' @export
#'
#' @examples
#' with(elmasses, 2*H + O)
#' with(elzs, H:O)
#'
elsymbols <- stats::setNames(as.list(pt$symbol), pt$symbol)

#' @rdname elsymbols
#' @export
elmasses <- stats::setNames(as.list(pt$mass), pt$symbol)

#' @rdname elsymbols
#' @export
elzs <- stats::setNames(as.list(pt$z), pt$symbol)

#' @rdname elsymbols
#' @export
elnames <- stats::setNames(as.list(pt$name), pt$symbol)

#' @rdname elsymbols
#' @export
elgroups <- stats::setNames(as.list(pt$group), pt$symbol)

#' @rdname elsymbols
#' @export
elperiods <- stats::setNames(as.list(pt$period), pt$symbol)

#' @rdname elsymbols
#' @export
elvalences <- stats::setNames(pt$valence, pt$symbol)
