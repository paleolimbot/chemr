
#' The Periodic Tibble
#'
#' @format A tibble with 118 observations of 6 variables
#' \describe{
#'   \item{z}{The atomic number of the element}
#'   \item{N}{The symbol of the element}
#'   \item{name}{The (UK) English name of the element}
#'   \item{group}{The group (column) of the element}
#'   \item{period}{The period (row) of the element}
#'   \item{mass}{The atomic mass of the element (g/mol)}
#' }
#'
#' @source \url{https://en.wikipedia.org/wiki/List_of_chemical_elements}
"pt"

# load the periodic tibble in the package namespace
data("pt", envir = environment())
