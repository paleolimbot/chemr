
#' Create a molecule object
#'
#' @param ... A named list of symbols and counts
#' @param charge The charge of the molecule
#' @param validate Validate elements in molecule
#' @param x An objec to be coerced to type molecule(s)
#'
#' @return An object of class molecule
#' @export
#'
#' @examples
#' m1 <- molecule(H=2, O=1)
#' m2 <- molecule(H=1, charge = 1)
#' molecules(m1, m2)
#'
molecule <- function(..., charge = 0L, validate = TRUE) {
  x <- c(...)
  x <- stats::setNames(as.double(x), names(x))
  x[is.na(x)] <- 1
  if(is.null(names(x))) stop("Arguments to molecule must be named")
  m <- new_molecule(x, charge = as.integer(charge),
                    mass = sum(x * elmass(names(x))))
  if(validate) validate_molecule(m)
  m
}

#' @rdname molecule
#' @export
molecules <- function(..., validate = TRUE) {
  ml <- new_molecules(lapply(list(...), as_molecule))
  if(validate) validate_molecules(ml)
  ml
}

#' @rdname molecule
#' @export
as_molecule <- function(x, ...) UseMethod("as_molecule")

#' @rdname molecule
#' @export
as_molecules <- function(x, ...) UseMethod("as_molecules")

#' @rdname molecule
#' @export
as_molecule.molecule <- function(x, ...) {
  x
}

#' @rdname molecule
#' @export
as_molecule.molecules <- function(x, ...) {
  if(length(x) > 1) warning("Using first of ", length(x), " molecules in x")
  x[[1]]
}

#' @rdname molecule
#' @export
as_molecule.character <- function(x, validate = TRUE, ...) {
  if(length(x) > 1) warning("More than one molecule in x. Did you mean as_molecules()?")
  m <- parse_molecules(x[1], validate = FALSE)[[1]]
  if(validate) validate_molecule(m)
  m
}

#' @rdname molecule
#' @export
as_molecule.formula <- function(x, validate = TRUE, ...) {
  vars <- all.vars(x)
  if(length(vars) > 1) warning("More than one molecule in formula. Did you mean as_molecules()?")
  m <- parse_molecules(vars[1], validate = FALSE)[[1]]
  if(validate) validate_molecule(m)
  m
}

#' @rdname molecule
#' @export
as_molecules.molecules <- function(x, ...) {
  x
}

#' @rdname molecule
#' @export
as_molecules.molecule <- function(x, validate = TRUE, ...) {
  molecules(x, validate = validate)
}

#' @rdname molecule
#' @export
as_molecules.character <- function(x, validate = TRUE, ...) {
  parse_molecules(x)
}

#' @rdname molecule
#' @export
as_molecules.formula <- function(x, validate = TRUE, ...) {
  parse_molecules(all.vars(x), validate = validate)
}


#' Create, validate molecule objects
#'
#' @param x A double vector or molecule object
#' @param charge An optional charge to assign to the molecule
#' @param mass A mass to assign the molecule
#'
#' @return An double vector with class "molecule"
#' @export
#'
#' @examples
#' m <- new_molecule(c(H=2, O=1), charge = 0L, mass = 18.01528)
#' validate_molecule(m)
#' is_molecule(m)
#'
new_molecule <- function(x, charge = 0L, mass = NA_real_) {
  if(!is.double(x)) stop("x must be a double vector")
  structure(x, charge = charge, mass = mass, class = "molecule")
}

#' @rdname molecule
#' @export
new_molecules <- function(x) {
  if(!is.list(x)) stop("x must be of type list")
  structure(x, class = "molecules")
}

#' @rdname molecule
#' @export
NA_molecule_ <- new_molecule(stats::setNames(numeric(0), character(0)),
                             charge = NA_integer_, mass = NA_real_)

#' @rdname molecule
#' @export
validate_molecule <- function(x) {
  # check class
  if(!is_molecule(x)) stop("x must be of type molecule")
  # check base type
  if(!is.double(x)) stop("x must be a double vector")
  # check names
  if(is.null(names(x))) stop("x must have names")
  # check attributes
  if(is.null(attr(x, "charge"))) stop("x is missing attr 'charge'")
  if(!is.integer(attr(x, "charge"))) stop("attr(x, 'charge') is not an integer")
  if(is.null(attr(x, "mass"))) stop("x is missing attr 'mass'")
  if(!is.double(attr(x, "mass"))) stop("attr(x, 'mass') is not a double")

  # check symbols
  bad_symbols <- names(x)[!is_symbol(names(x))]
  if(length(bad_symbols) > 0) stop("names(x) contained the following bad symbols: ",
                                   paste(bad_symbols, collpase = ", "))

  # return x, invisibly
  invisible(x)
}

#' @rdname molecule
#' @export
validate_molecules <- function(x) {
  if(!is.list(x)) stop("x must be a list")
  validation <- lapply(x, function(x) try(validate_molecule(x), silent = TRUE))
  val_error <- vapply(validation, inherits, "try-error", FUN.VALUE = logical(1))
  if(any(val_error)) stop("Molecules at positions ",
                          paste(which(val_error), collapse = ", "),
                          " are invalid molecule objects")
  invisible(x)
}

#' @rdname molecule
#' @export
is_molecule <- function(x) {
  inherits(x, "molecule")
}

#' @rdname molecule
#' @export
is_molecules <- function(x) {
  inherits(x, "molecules")
}

#' Coerce molecule(s) to character
#'
#' @param x A molecule(s) object
#' @param ... Ignored
#'
#' @return A character vector
#' @export
#'
#' @examples
#' print(as_molecule(~H2O))
#' print(as_molecules(~H2O))
#' as.character(NA_molecule_)
#' as.character(molecules(NA_molecule_))
#'
print.molecule <- function(x, ...) {
  cat("<molecule>", as.character(x))
  invisible(x)
}

#' @rdname print.molecule
#' @export
print.molecules <- function(x, ...) {
  cat("<molecules>\n")
  print(as.character(x), quote = FALSE)
  invisible(x)
}

#' @rdname print.molecule
#' @export
as.character.molecule <- function(x, ...) {
  if(identical(x, NA_molecule_)) return("<NA_molecule_>")
  counts <- ifelse(x == 1, "", unclass(x)) # will be character
  charge <- charge(x)
  charge <- ifelse(charge == 1, "+",
                   ifelse(charge == -1, "-",
                          ifelse(charge > 0, paste("+", charge),
                                 ifelse(charge == 0, "", charge))))
  mat <- rbind(names(x), counts)
  dim(mat) <- NULL
  paste0(paste(mat, collapse = ""), charge)
}

#' @rdname print.molecule
#' @export
as.character.molecules <- function(x, ...) {
  vapply(x, as.character.molecule, character(1))
}


#' Access properties of a molecule object
#'
#' @param x A molecule object or symbol string
#'
#' @return The mass of the molecule or element in g/mol
#' @export
#'
#' @examples
#' mass("H")
#' m <- as_molecule("H2O")
#' mass(m)
#'
mass <- function(x) UseMethod("mass")

#' @rdname mass
#' @export
mass.default <- function(x) {
  mass(as_molecules(x))
}

#' @rdname mass
#' @export
mass.molecule <- function(x) {
  attr(x, "mass")
}

#' @rdname mass
#' @export
mass.molecules <- function(x) {
  vapply(x, attr, "mass", FUN.VALUE = double(1))
}

#' @rdname mass
#' @export
charge <- function(x) UseMethod("charge")

#' @rdname mass
#' @export
charge.default <- function(x) {
  charge(as_molecules(x))
}

#' @rdname mass
#' @export
charge.molecule <- function(x) {
  attr(x, "charge")
}

#' @rdname mass
#' @export
charge.molecules <- function(x) {
  vapply(x, attr, "charge", FUN.VALUE = integer(1))
}

# internal function to parse molecule text
.el_regex <- "([A-Z][a-z]{0,2})([0-9]*)"
.charge_regex <- "([-+][0-9]*)$"

parse_molecules <- function(txt, validate = TRUE, na = c("NA", "<NA_molecule_>", "")) {
  mol_na <- txt %in% na
  matches <- stringr::str_match_all(txt, .el_regex)
  # process charges
  charge_matches <- stringr::str_extract(txt, .charge_regex)
  # charge of "-" is -1, "+" is "+1"
  charges <- ifelse(charge_matches == "-", "-1", charge_matches)
  charges <- ifelse(charges == "+", 1L, suppressWarnings(as.integer(charges)))
  # no charge is zero
  charges[is.na(charges)] <- 0L

  mol_list <- lapply(seq_along(matches), function(i) {
    mol_string <- txt[i]
    # check for NA
    if(is.na(mol_string) || mol_na) return(NA_molecule_)
    match <- matches[[i]]

    # check for full match
    all_matches <- c(match[, 1, drop = TRUE], stats::na.omit(charge_matches[i]))
    if(nchar(paste(all_matches, collapse = "")) != nchar(mol_string)) {
      warning("Bad molecule text: ", mol_string)
      return(NA_molecule_)
    }

    # extract symbols
    symbols <- match[, 2, drop = TRUE]

    # if counts is NA it should be 1
    counts <- as.double(match[, 3, drop = TRUE])
    counts[is.na(counts)] <- 1

    # els is a named version of counts
    els <- stats::setNames(counts, symbols)

    # return els, classed as a molecule
    new_molecule(els, charge = charges[i], mass = sum(elmass(symbols) * counts))
  })

  # return mol_list classed as molecules
  ml <- new_molecules(mol_list)
  if(validate) validate_molecules(ml)
  ml
}
