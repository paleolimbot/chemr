
# load the periodic tibble
data("pt", envir = environment())

is_symbol <- function(x) {
  x %in% as.character(pt$symbol)
}

elmass <- function(x) {
  stats::setNames(pt$mass[match(x, pt$symbol)], x)
}

