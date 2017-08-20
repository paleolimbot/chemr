
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
#' m1 <- molecule_single(H=2, O=1)
#' m2 <- molecule_single(H=1, charge = 1)
#' mol(m1, m2)
#'
mol <- function(..., validate = TRUE) {
  ml <- new_mol(lapply(list(...), as_molecule_single))
  if(validate) validate_mol(ml)
  ml
}

#' @rdname mol
#' @export
molecule_single <- function(..., charge = 0L, validate = TRUE) {
  x <- c(...)
  x <- stats::setNames(as.double(x), names(x))
  x[is.na(x)] <- 1
  if(is.null(names(x))) stop("Arguments to molecule must be named")
  m <- new_molecule_single(x, charge = as.integer(charge),
                    mass = sum(x * elmass(names(x))))
  if(validate) validate_molecule_single(m)
  m
}

#' @rdname mol
#' @export
as_molecule_single <- function(x, ...) UseMethod("as_molecule_single")

#' @rdname mol
#' @export
as_mol <- function(x, ...) UseMethod("as_mol")

#' @rdname mol
#' @export
as.mol <- function(x, ...) as_mol(x, ...)

#' @rdname mol
#' @export
as_molecule_single.molecule_single <- function(x, ...) {
  x
}

#' @rdname mol
#' @export
as_molecule_single.mol <- function(x, ...) {
  if(length(x) > 1) warning("Using first of ", length(x), " mol in x")
  x[[1]]
}

#' @rdname mol
#' @export
as_molecule_single.character <- function(x, validate = TRUE, ...) {
  if(length(x) > 1) warning("More than one molecule in x. Did you mean as_mol()?")
  m <- parse_mol(x[1], validate = FALSE)[[1]]
  if(validate) validate_molecule_single(m)
  m
}

#' @rdname mol
#' @export
as_molecule_single.formula <- function(x, validate = TRUE, ...) {
  vars <- all.vars(x)
  if(length(vars) > 1) warning("More than one molecule in formula. Did you mean as_mol()?")
  m <- parse_mol(vars[1], validate = FALSE)[[1]]
  if(validate) validate_molecule_single(m)
  m
}

#' @rdname mol
#' @export
as_mol.mol <- function(x, ...) {
  x
}

#' @rdname mol
#' @export
as_mol.molecule_single <- function(x, validate = TRUE, ...) {
  mol(x, validate = validate)
}

#' @rdname mol
#' @export
as_mol.character <- function(x, validate = TRUE, ...) {
  parse_mol(x)
}

#' @rdname mol
#' @export
as_mol.formula <- function(x, validate = TRUE, ...) {
  parse_mol(all.vars(x), validate = validate)
}

#' Combine, subset molecule(s) objects
#'
#' @param x A mol object
#' @param i The index to extract
#' @param ... Objects to combine
#'
#' @return A mol object
#' @export
#'
#' @examples
#' c(as_mol(~H2O), as_mol(~`NH3`))
#' mols <- as_mol(c("H2O", "NH3"))
#' mols[1]
#'
c.molecule_single <- function(...) {
  do.call(c.mol, lapply(list(...), as_mol))
}

#' @rdname c.molecule_single
#' @export
c.mol <- function(...) {
  args <- lapply(list(...), unclass)
  new_mol(do.call(c, args))
}

#' @rdname c.molecule_single
#' @export
`[.mol` <- function(x, i, ...) {
  new_mol(unclass(x)[i, ...])
}

#' Create, validate molecule objects
#'
#' @param x A double vector or molecule object
#' @param charge An optional charge to assign to the molecule
#' @param mass A mass to assign the molecule
#'
#' @return An double vector with class "molecule_single"
#' @export
#'
#' @examples
#' m <- new_molecule_single(c(H=2, O=1), charge = 0L, mass = 18.01528)
#' validate_molecule_single(m)
#' is_molecule_single(m)
#'
new_molecule_single <- function(x, charge = 0L, mass = NA_real_) {
  if(!is.double(x)) stop("x must be a double vector")
  structure(x, charge = charge, mass = mass, class = "molecule_single")
}

#' @rdname mol
#' @export
new_mol <- function(x) {
  if(!is.list(x)) stop("x must be of type list")
  structure(x, class = "mol")
}

#' @rdname mol
#' @export
NA_molecule_ <- new_molecule_single(stats::setNames(numeric(0), character(0)),
                                    charge = NA_integer_, mass = NA_real_)

#' @rdname mol
#' @export
validate_molecule_single <- function(x) {
  # check class
  if(!is_molecule_single(x)) stop("x must be of type molecule")
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

#' @rdname mol
#' @export
validate_mol <- function(x) {
  if(!is.list(x)) stop("x must be a list")
  validation <- lapply(x, function(x) try(validate_molecule_single(x), silent = TRUE))
  val_error <- vapply(validation, inherits, "try-error", FUN.VALUE = logical(1))
  if(any(val_error)) stop("mol at positions ",
                          paste(which(val_error), collapse = ", "),
                          " are invalid molecule objects")
  invisible(x)
}

#' @rdname mol
#' @export
is_molecule_single <- function(x) {
  inherits(x, "molecule_single")
}

#' @rdname mol
#' @export
is_mol <- function(x) {
  inherits(x, "mol")
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
#' print(as_molecule_single(~H2O))
#' print(as_mol(~H2O))
#' as.character(NA_molecule_)
#' as.character(mol(NA_molecule_))
#'
print.molecule_single <- function(x, ...) {
  cat("<molecule_single>", as.character(x))
  invisible(x)
}

#' @rdname print.molecule_single
#' @export
print.mol <- function(x, ...) {
  cat("<mol>\n")
  print(as.character(x), quote = FALSE)
  invisible(x)
}

#' @rdname print.molecule_single
#' @export
as.character.molecule_single <- function(x, ...) {
  if(identical(x, NA_molecule_)) return("<NA_molecule_>")
  counts <- ifelse(x == 1, "", unclass(x)) # will be character
  charge <- charge(x)
  charge <- ifelse(charge == 1, "+",
                   ifelse(charge == -1, "-",
                          ifelse(charge > 0, paste0("+", charge),
                                 ifelse(charge == 0, "", charge))))
  mat <- rbind(names(x), counts)
  dim(mat) <- NULL
  paste0(paste(mat, collapse = ""), charge)
}

#' @rdname print.molecule_single
#' @export
as.character.mol <- function(x, ...) {
  vapply(x, as.character.molecule_single, character(1))
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
#' m <- as_molecule_single("H2O")
#' mass(m)
#'
mass <- function(x) UseMethod("mass")

#' @rdname mass
#' @export
mass.default <- function(x) {
  mass(as_mol(x))
}

#' @rdname mass
#' @export
mass.molecule_single <- function(x) {
  attr(x, "mass")
}

#' @rdname mass
#' @export
mass.mol <- function(x) {
  vapply(x, attr, "mass", FUN.VALUE = double(1))
}

#' @rdname mass
#' @export
charge <- function(x) UseMethod("charge")

#' @rdname mass
#' @export
charge.default <- function(x) {
  charge(as_mol(x))
}

#' @rdname mass
#' @export
charge.molecule_single <- function(x) {
  attr(x, "charge")
}

#' @rdname mass
#' @export
charge.mol <- function(x) {
  vapply(x, attr, "charge", FUN.VALUE = integer(1))
}

# internal function to parse molecule text
.el_regex <- "([A-Z][a-z]{0,2})([0-9]*)"
.charge_regex <- "([-+][0-9]*)$"

parse_mol <- function(txt, validate = TRUE, na = c("NA", "<NA_molecule_>", "")) {
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
    if(is.na(mol_string) || mol_na[i]) return(NA_molecule_)
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
    new_molecule_single(els, charge = charges[i], mass = sum(elmass(symbols) * counts))
  })

  # return mol_list classed as mol
  ml <- new_mol(mol_list)
  if(validate) validate_mol(ml)
  ml
}
