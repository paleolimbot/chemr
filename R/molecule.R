
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
  # NA counts are 1
  x[is.na(x)] <- 1
  # zero counts are removed
  x <- x[x != 0]
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
  if(!is.logical(validate)) stop("Invalid value for validate: ", validate)
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
  if(!is.logical(validate)) stop("Invalid value for validate: ", validate)
  parse_mol(x, validate = validate)
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
#' @param times The number of times to repeat the molecule
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
  l <- unclass(x)[i, ...]
  # NULLs should be NA_molecule
  l[vapply(l, is.null, logical(1))] <- list(NA_molecule_)
  new_mol(l)
}

#' @rdname c.molecule_single
#' @export
rep.molecule_single <- function(x, times, ...) {
  new_mol(rep(list(x), times))
}

#' @rdname c.molecule_single
#' @export
rep.mol <- function(x, times, ...) {
  new_mol(rep(unclass(x), times))
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
is.na.molecule_single <- function(x) {
  identical(x, NA_molecule_)
}

#' @rdname mol
#' @export
is.na.mol <- function(x) {
  vapply(x, identical, NA_molecule_, FUN.VALUE = logical(1))
}

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
  bad_symbols <- names(x)[!is_element(names(x))]
  if(length(bad_symbols) > 0) stop("names(x) contained the following bad symbols: ",
                                   paste(bad_symbols, collpase = ", "))

  # check counts
  if(any(x <= 0)) stop("All counts must be > 0")

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
  counts <- ifelse(unclass(x) == 1, "", unclass(x)) # will be character
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


#' Molecule arithmetic
#'
#' @param x A mol object
#' @param y The object with which to perform arithmetic
#' @param ... Objects to combine
#'
#' @return A mol object
#' @export
#' @rdname arithmetic
#'
#' @examples
#' m <- as_mol("H2O")
#' m*2
#' m + as_mol("H+")
#'
`*.molecule_single` <- function(x, y) {
  # NA x or y results in NA_molecule_ output
  if(is.na(x) || is.na(y)) return(NA_molecule_)
  # multiply coefficients, charge, mass by y
  # if x is a mol_single:
  if(is_molecule_single(x)) {
    m <- new_molecule_single(unclass(x) * y,
                             charge = charge(x) * y, mass = mass(x) * y)
  } else {
    m <- new_molecule_single(unclass(y) * x,
                             charge = charge(y) * x, mass = mass(y) * x)
  }

  # remove zero counts
  remove_zero_counts.molecule_single(m)
}

#' @rdname arithmetic
#' @export
`/.molecule_single` <- function(x, y) {
  # NA x or y results in NA_molecule_ output
  if(is.na(x) || is.na(y)) return(NA_molecule_)
  # divide coefficients, charge, mass by y
  if(!is_molecule_single(x)) stop("Can't divide by a molecule_single")
  if(!is.numeric(y)) stop("Can't divide a molecule_single by an object of type ", class(y)[1])
  # zero division results in NA molecule
  if(y == 0) {
    warning("Divide by zero in /.molecule_single -> NA_molecule")
    return(NA_molecule_)
  }
  new_molecule_single(unclass(x) / y,
                      charge = charge(x) / y, mass = mass(x) / y)
}

#' @rdname arithmetic
#' @export
`+.molecule_single` <- function(x, y) {
  # turn X,Y into a molecule_single
  x <- as_molecule_single(x)
  y <- as_molecule_single(y)
  # NA x or y results in NA_molecule_ output
  if(is.na(x) || is.na(y)) return(NA_molecule_)
  # combine raw coefficient vectors, add charges, masses
  m <- new_molecule_single(c(unclass(x), unclass(y)),
                           charge = charge(x) + charge(y),
                           mass = mass(x) + mass(y))
  # remove zero counts
  remove_zero_counts.molecule_single(m)
}

#' @rdname arithmetic
#' @export
combine_molecules <- function(...) {
  vals <- list(...)
  if(any(vapply(vals, is_mol, logical(1)))) {
    # vectorize using mapply
    margs <- c(list(combine_molecules), lapply(vals, as_mol), list(SIMPLIFY = FALSE))
    return(new_mol(do.call(mapply, margs)))
  }
  mols <- mol(..., validate = FALSE)
  if(any(is.na(mols))) return(NA_molecule_)
  unclassed <- lapply(mols, unclass)
  m <- new_molecule_single(do.call(c, unclassed),
                           charge = sum(charge(mols)),
                           mass = sum(mass(mols)))
  # remove zero counts
  remove_zero_counts.molecule_single(m)
}

#' @rdname arithmetic
#' @export
`==.molecule_single` <- function(x, y) {
  # turn X,Y into a molecule_single
  x <- as_molecule_single(x)
  y <- as_molecule_single(y)
  # combine raw coefficient vectors, add charges, masses
  bare_x <- unname(unclass(x))
  bare_y <- unname(unclass(y))
  all(names(x) == names(y)) && all(bare_x == bare_y) && (charge(x) == charge(y))
}

#' @rdname arithmetic
#' @export
`*.mol` <- function(x, y) {
  result <- mapply(`*.molecule_single`, x, y, SIMPLIFY = FALSE)
  new_mol(result)
}

#' @rdname arithmetic
#' @export
`/.mol` <- function(x, y) {
  result <- mapply(`/.molecule_single`, x, y, SIMPLIFY = FALSE)
  new_mol(result)
}

#' @rdname arithmetic
#' @export
`+.mol` <- function(x, y) {
  result <- mapply(`+.molecule_single`, x, y, SIMPLIFY = FALSE)
  new_mol(result)
}

#' @rdname arithmetic
#' @export
`==.mol` <- function(x, y) {
  mapply(`==.molecule_single`, x, y, SIMPLIFY = TRUE)
}


#' Simplify chemr objects
#'
#' @param x A mol, molecule_single, or reaction object
#' @param ... Unused
#'
#' @return An object with the same class as x
#' @export
#'
#' @examples
#' simplify(as_mol("CHOOOH"))
#'
simplify <- function(x, ...) UseMethod("simplify")

#' @rdname simplify
#' @export
remove_zero_counts <- function(x, ...) UseMethod("remove_zero_counts")

#' @rdname simplify
#' @export
simplify.molecule_single <- function(x, ...) {
  unique_names <- unique(names(x))
  obj <- tapply(x, names(x), sum)[unique_names]
  if(length(obj) == 0) {
    obj <- stats::setNames(numeric(0), character(0))
  } else {
    obj <- stats::setNames(as.double(obj), names(obj))
  }
  new_molecule_single(obj, charge = charge(x), mass = mass(x))
}

#' @rdname simplify
#' @export
remove_zero_counts.molecule_single <- function(x, ...) {
  new_molecule_single(x[x != 0], charge = charge(x), mass = mass(x))
}

#' @rdname simplify
#' @export
simplify.mol <- function(x, ...) {
  new_mol(lapply(x, simplify.molecule_single))
}

#' @rdname simplify
#' @export
remove_zero_counts.mol <- function(x, ...) {
  new_mol(lapply(x, remove_zero_counts.molecule_single))
}

# internal function to parse molecule text
.el_regex <- "([A-Z][a-z]{0,2})([0-9.]*)"
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

    # if counts == 0, it should be removed
    non_zero_counts <- counts != 0
    counts <- counts[non_zero_counts]
    symbols <- symbols[non_zero_counts]

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


