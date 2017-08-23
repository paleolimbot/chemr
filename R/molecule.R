
#' Create a molecule object
#'
#' @param ... A named list of symbols and counts
#' @param charge The charge of the molecule
#' @param count Optional count to associate with the object
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
molecule_single <- function(..., charge = 0, count = NA_real_, validate = TRUE) {
  mol_list <- list(...)
  # check for empty molecule
  if(length(mol_list) == 0) return(new_molecule_single(list()))
  if(is.null(names(mol_list))) stop("At least one argument to molecule_single must be named")

  # make recursive list of x
  mol_list_parsed <- lapply(mol_list, function(x) {
    if(is.list(x)) {
      do.call(molecule_single, c(x, list(validate = validate)))
    } else if(is.null(x)) {
      NULL
    } else if(is.na(x)) {
      stop("Invalid value in molecule_single: NA")
    } else if(is.numeric(x)) {
      as.double(x)
    } else {
      stop("Invalid type in molcule_single: ", class(x)[1])
    }
  })
  # parse non-element names as molecules
  sub_mols <- vapply(mol_list_parsed, is_molecule_single, logical(1))
  non_elements <- !stringr::str_detect(names(mol_list_parsed), "^[A-Z][a-z]?$") & !sub_mols
  non_element_counts <- mol_list_parsed[non_elements]
  mol_list_parsed[non_elements] <- mapply(as_molecule_single,
                                          names(mol_list_parsed)[non_elements],
                                          count = mol_list_parsed[non_elements],
                                          SIMPLIFY = FALSE)

  # sub molecules don't have names
  sub_mols <- vapply(mol_list_parsed, is_molecule_single, logical(1))
  names(mol_list_parsed) <- ifelse(sub_mols, "", names(mol_list))

  # create, validate molecule object
  mol <- new_molecule_single(mol_list_parsed, count = count, charge = charge)
  if(validate) validate_molecule_single(mol)
  mol
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
  parse_single(x[1], ..., validate = validate)
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

#' @rdname c.molecule_single
#' @export
unique.mol <- function(x, ...) {
  new_mol(NextMethod())
}

#' Create, validate molecule objects
#'
#' @param x A double vector or molecule object
#' @param charge An optional charge to assign to the molecule
#' @param count An optional count to assign to the molecule
#'
#' @return An double vector with class "molecule_single"
#' @export
#'
#' @examples
#' m <- new_molecule_single(list(H=2, O=1), charge = 0)
#' validate_molecule_single(m)
#' is_molecule_single(m)
#'
new_molecule_single <- function(x, charge = 0, count = NA_real_) {
  if(!is.list(x)) stop("x must be a list")
  structure(x, charge = charge, count = count, class = "molecule_single")
}

#' @rdname mol
#' @export
new_mol <- function(x) {
  if(!is.list(x)) stop("x must be of type list")
  structure(x, class = "mol")
}

#' @rdname mol
#' @export
NA_molecule_ <- new_molecule_single(stats::setNames(list(), character(0)),
                                    charge = NA_integer_, count = NA_real_)

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
  if(!is_molecule_single(x)) stop("x must be of type molecule_single")
  # check base type
  if(!is.list(x)) stop("x must be a list")
  # check names
  if(is.null(names(x))) stop("x must have names")
  # check attributes
  if(is.null(attr(x, "charge"))) stop("x is missing attr 'charge'")
  if(!is.numeric(attr(x, "charge"))) stop("attr(x, 'charge') is not numeric")
  if(is.null(attr(x, "count"))) stop("x is missing attr 'count'")
  if(!is.double(attr(x, "count"))) stop("attr(x, 'count') is not a double")

  # check symbols
  bad_symbols <- names(x)[!is_element(names(x)) & (names(x) != "")]
  if(length(bad_symbols) > 0) stop("names(x) contained the following bad symbols: ",
                                   paste(bad_symbols, collpase = ", "))

  # check sub molecules
  lapply(x[names(x) == ""], validate_molecule_single)

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
  cat("<molecule_single>", as.character(x, ...))
  invisible(x)
}

#' @rdname print.molecule_single
#' @export
print.mol <- function(x, ...) {
  cat("<mol>\n")
  print(as.character(x, ...), quote = FALSE)
  invisible(x)
}

#' @rdname print.molecule_single
#' @export
as.character.molecule_single <- function(x, ...) {
  if(identical(x, NA_molecule_)) return(NA_character_)

  counts <- vapply(x, function(el) {
    if(is.numeric(el)) {
      el
    } else if(is_molecule_single(el)) {
      attr(el, "count")
    } else {
      stop("Invalid type in molecule_single: ", class(el)[1])
    }
  }, numeric(1))

  symbols <- names(x)
  symbols[symbols == ""] <- paste0(
    "(",
    vapply(x[symbols == ""], as.character, trim = TRUE, ..., FUN.VALUE = character(1)),
    ")"
  )

  counts <- ifelse(counts == 1, "", format(counts, trim = TRUE, ...)) # will be character
  charge <- charge(x)
  charge <- ifelse(charge == 1, "+",
                   ifelse(charge == -1, "-",
                          ifelse(charge > 0, paste0("+", format(charge, trim = TRUE, ...)),
                                 ifelse(charge == 0, "", format(charge, trim = TRUE, ...)))))
  mat <- rbind(symbols, counts)
  dim(mat) <- NULL
  paste0(paste(mat, collapse = ""), charge)
}

#' @rdname print.molecule_single
#' @export
as.character.mol <- function(x, ...) {
  vapply(x, as.character.molecule_single, ..., FUN.VALUE = character(1))
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
  counts <- vapply(x, function(el) {
    if(is.list(el)) {
      attr(el, "count")
    } else {
      el
    }
  }, numeric(1))

  masses <- elmass(names(x))
  masses[names(x) == ""] <- vapply(x[names(x) == ""], mass.molecule_single, numeric(1))

  sum(counts * masses)
}

#' @rdname mass
#' @export
mass.mol <- function(x) {
  vapply(x, mass, FUN.VALUE = double(1))
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
  vapply(x, attr, "charge", FUN.VALUE = numeric(1))
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
    m <- lapply(x, function(el) {
      if(is.list(el)) {
        attr(el, "count") <- attr(el, "count") * y
        el
      } else {
        el * y
      }
    })

    new_molecule_single(m, charge = charge(x) * y, count = attr(x, "count"))
  } else {
    y * x
  }
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
  # use mul operator
  x * (1/y)
}

#' @rdname arithmetic
#' @export
`+.molecule_single` <- function(x, y) {
  if(missing(y)) return(x) # unary operator
  # turn X,Y into a molecule_single
  x <- as_molecule_single(x)
  y <- as_molecule_single(y)
  # NA x or y results in NA_molecule_ output
  if(is.na(x) || is.na(y)) return(NA_molecule_)
  # combine raw coefficient vectors, add charges, masses
  m <- new_molecule_single(c(unclass(x), unclass(y)),
                           charge = charge(x) + charge(y))
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
                           count = NA_real_)
  # remove zero counts
  remove_zero_counts.molecule_single(m)
}

#' @rdname arithmetic
#' @export
`==.molecule_single` <- function(x, y) {
  # turn X,Y into a molecule_single
  x <- as_molecule_single(x)
  y <- as_molecule_single(y)
  # compare character representations
  as.character(x) == as.character(y)
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
#' @param tol Tolerance for zero detection
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

  # simplify sub mols
  sub_mols <- vapply(x, is_molecule_single, logical(1))
  x[sub_mols] <- lapply(x[sub_mols], function(el) {
    simplify(el) * attr(el, "count")
  })

  # turn x into a named vector
  if(length(x) > 0) {
    x_simple <- unlist(x)
  } else {
    x_simple <- x
  }

  unique_names <- unique(names(x_simple))
  obj <- tapply(x_simple, names(x_simple), sum)[unique_names]
  if(length(obj) == 0) {
    obj <- stats::setNames(numeric(0), character(0))
  } else {
    obj <- stats::setNames(as.double(obj), names(obj))
  }
  new_molecule_single(as.list(obj), charge = charge(x), count = attr(x, "count"))
}

#' @rdname simplify
#' @export
remove_zero_counts.molecule_single <- function(x, tol = .Machine$double.eps^0.5, ...) {
  counts <- vapply(x, function(el) {
    if(is.list(el)) {
      attr(el, "count")
    } else {
      el
    }
  }, numeric(1))
  new_molecule_single(x[abs(counts) >= tol], charge = charge(x), count = attr(x, "count"))
}

#' @rdname simplify
#' @export
simplify.mol <- function(x, ...) {
  new_mol(lapply(x, simplify.molecule_single))
}

#' @rdname simplify
#' @export
remove_zero_counts.mol <- function(x, tol = .Machine$double.eps^0.5, ...) {
  new_mol(lapply(x, remove_zero_counts.molecule_single, tol = tol))
}

# internal function to parse molecule text
.el_regex <- "([A-Z][a-z]?|\\(.+?\\))(-?[0-9.]*)"
# .sub_mol_regex <- "\\(.+?\\)(-?[0-9.]*)"
.mol_regex <- "^(.+?)([-+][0-9]*)?$"

parse_single <- function(txt, validate = TRUE, na = c("NA", ""), count = NA_real_) {
  if(txt %in% na || is.na(txt)) return(NA_molecule_)

  # match molecule
  full_match <- stringr::str_match(txt, .mol_regex)
  if(is.na(full_match[, 1, drop = TRUE])) warning("Bad molecule text: ", txt)
  charge_str <- full_match[, 3, drop = TRUE]
  mol_string <- full_match[, 2, drop = TRUE]
  charge <- ifelse(is.na(charge_str) || (charge_str == ""), 0,
                   ifelse(charge_str == "-", -1,
                          ifelse(charge_str == "+", 1,
                                 suppressWarnings(as.numeric(charge_str)))))

  # match elements, sub-molecules
  match <- stringr::str_match_all(mol_string, .el_regex)[[1]]

  # check for full match
  all_matches <- c(match[, 1, drop = TRUE])
  if(nchar(paste(all_matches, collapse = "")) != nchar(mol_string)) {
    warning("Bad molecule text: ", mol_string)
    return(NA_molecule_)
  }

  # extract symbols
  symbols <- stringr::str_replace_all(match[, 2, drop = TRUE], "[\\(\\)]", "")

  # if counts is NA it should be 1
  counts <- as.double(match[, 3, drop = TRUE])
  counts[is.na(counts)] <- 1

  # make list of names, counts, call molecule_single
  do.call(molecule_single,
          c(as.list(stats::setNames(counts, symbols)),
            list(count = count, charge = charge, validate = validate)))
}

parse_mol <- function(txt, validate = TRUE, na = c("NA", "")) {

  mol_list <- lapply(txt, parse_single, validate = FALSE, na = na)

  # return mol_list classed as mol
  ml <- new_mol(mol_list)
  if(validate) validate_mol(ml)
  ml
}

#' Data frame reprsentation of a molecule object
#'
#' @param x A molecule object
#' @param ... Ignored
#'
#' @return A tibble with one row per molecule
#' @export
#' @rdname molecule_tibble
#'
#' @importFrom tibble as_tibble
#'
#' @examples
#' as.data.frame(as_mol(~H2O))
#' library(tibble)
#' as_tibble(as_mol(~H2O))
#' as_tibble(mol("H2O", "NH3"))
#'
as.data.frame.molecule_single <- function(x, ...) {
  as.data.frame(as_tibble.molecule_single(x, ...))
}

#' @rdname molecule_tibble
#' @export
as_tibble.molecule_single <- function(x, ...) {
  as_tibble.mol(as_mol(x), ...)
}

#' @rdname molecule_tibble
#' @export
as_tibble.mol <- function(x, ...) {
  df <- cbind(
    tibble::tibble(
      mol = x,
      mol_text = as.character(x),
      mass = mass(x),
      charge = charge(x)
    ),
    element_tbl_mol(x)
  )
  tibble::as_tibble(df)
}

#' @rdname molecule_tibble
#' @export
as.data.frame.mol <- function(x, ...) {
  as.data.frame(as_tibble.mol(x, ...))
}

#' @rdname molecule_tibble
#' @export
as.matrix.molecule_single <- function(x, ...) {
  as.matrix.mol(as_mol(x), ...)
}

#' @rdname molecule_tibble
#' @export
as.matrix.mol <- function(x, ...) {
  matrix <- as.matrix(element_tbl_mol(x))
  rownames(matrix) <- make.unique(as.character(x), sep = "_")
  matrix
}

# at the base of many of the above functions
element_tbl <- function(mol_single) {
  m_simple <- unclass(remove_zero_counts(simplify(mol_single)))
  # make sure there is always one row
  m_simple$.dummy <- 1
  tibble::as_tibble(m_simple)
}

# multi version of above
element_tbl_mol <- function(mol, fill_empty = 0) {
  if(length(mol) == 0) return(tibble::tibble())
  # make elements frame
  elements <- purrr::map_df(mol, element_tbl)
  # remove .dummy column
  elements <- elements[setdiff(names(elements), ".dummy")]
  # make NA values 0 in elements frame
  elements[] <- lapply(elements, function(x) {
    replace(x, is.na(x), fill_empty)
  })
  # return elements
  elements
}
