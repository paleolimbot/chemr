
#' Create reaction objects
#'
#' @param lhs The lefthand side
#' @param counts_lhs The counts for molecules on the lefthand side
#' @param rhs The righthand side
#' @param counts_rhs The counts for molecules on the righthand side
#' @param log_k The base 10 logarithm of the reaction quotient
#' @param validate Flag to validate molecules in the reaction
#' @param x An object to convert to a reaction object
#' @param coefficient coefficient corresponding to as_reaction input
#' @param ... Passed to/from methods
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' reaction(c("H2O", "H+"), "H3O+", log_k = 14)
#'
reaction <- function(lhs, rhs, counts_lhs = 1, counts_rhs = 1, log_k = NA_real_, validate = TRUE) {
  lhs <- as_mol(lhs, validate = validate)
  # check length of counts_lhs
  if((length(counts_lhs) != length(lhs)) && (length(counts_lhs) != 1)) {
    stop("counts_lhs must be 1 or length(rhs)")
  }
  counts_lhs <- rep_len(counts_lhs, length(lhs))
  rhs <- as_mol(rhs, validate = validate)
  # check length of counts_rhs
  if((length(counts_rhs) != length(rhs)) && (length(counts_rhs) != 1)) {
    stop("counts_lhs must be 1 or length(rhs)")
  }
  counts_rhs <- rep_len(counts_rhs, length(rhs)) * -1

  log_k <- as.numeric(log_k)

  r <- new_reaction(list(mol = c(lhs, rhs), coefficient = c(counts_lhs, counts_rhs)), log_k = log_k)
  # need to validate reaction counts regardless of validate arg
  validate_reaction(r)
  # return r
  r
}

#' @rdname reaction
#' @export
reaction_list <- function(..., validate = TRUE) {
  rl <- new_reaction_list(lapply(list(...), as_reaction))
  if(validate) validate_reaction_list(rl)
  rl
}

#' @rdname reaction
#' @export
as_reaction <- function(x, ...) UseMethod("as_reaction")

#' @rdname reaction
#' @export
as.reaction <- function(x, ...) as_reaction(x, ...)

#' @rdname reaction
#' @export
as_reaction.reaction <- function(x, ...) {
  x
}

#' @rdname reaction
#' @export
as_reaction.mol <- function(x, coefficient = rep_len(1, length(x)), log_k = NA_real_, validate = TRUE, ...) {
  if(length(coefficient) != length(x)) {
    stop("length(coefficient) is not equal to length(x)")
  }
  r <- new_reaction(list(mol = x, coefficient = coefficient), log_k = log_k)
  if(validate) validate_reaction(r)
  r
}

#' @rdname reaction
#' @export
as_reaction.formula <- function(x, log_k = NA_real_, validate = TRUE, ...) {

  counts_lhs <- count_names(x[[2]])
  counts_rhs <- count_names(x[[3]])

  mols_lhs <- as_mol(names(counts_lhs), validate = validate, ...)
  mols_rhs <- as_mol(names(counts_rhs), validate = validate, ...)

  reaction(lhs = mols_lhs, rhs = mols_rhs,
           counts_lhs =  stats::setNames(counts_lhs, NULL),
           counts_rhs = stats::setNames(counts_rhs, NULL),
           log_k = log_k, validate = FALSE)
}

#' @rdname reaction
#' @export
as_reaction.character <- function(x, log_k = NA_real_, validate = TRUE, ...) {
  # character representation is the same the phreeqc reaction format
  # CaCl2 = Ca+2 + 2Cl
  sides <- stringr::str_split(x, "\\s*=\\s*")[[1]]
  if(length(sides) != 2) stop("Reaction must have exactly two sides")
  side_obj <- lapply(sides, parse_side, validate = validate)
  # create reaction from mols and counts
  reaction(lhs = side_obj[[1]]$mols,
           rhs = side_obj[[2]]$mols,
           counts_lhs = side_obj[[1]]$counts,
           counts_rhs = side_obj[[2]]$counts,
           log_k = log_k, validate = FALSE)

}

#' @rdname reaction
#' @export
as_reaction_list <- function(x, ...) UseMethod("as_reaction_list")

#' @rdname reaction
#' @export
as.reaction_list <- function(x, ...) UseMethod("as_reaction_list")

#' @rdname reaction
#' @export
as_reaction_list.reaction_list <- function(x, ...) {
  x
}

#' @rdname reaction
#' @export
as_reaction_list.list <- function(x, validate = TRUE, ...) {
  do.call(reaction_list, c(x, list(validate = validate, ...)))
}

#' @rdname reaction
#' @export
as_reaction_list.reaction <- function(x, validate = TRUE, ...) {
  reaction_list(x, validate = validate)
}

#' @rdname reaction
#' @export
as_reaction_list.character <- function(x, validate = TRUE, ...) {
  do.call(reaction_list, c(x, list(validate = validate, ...)))
}

parse_side <- function(side, validate = TRUE) {
  components <- stringr::str_split(side, "\\s+\\+\\s*")[[1]]
  # component is a count plus a molecule id
  component_match <- stringr::str_match(components, "^\\s*([0-9.]*)\\s*(.*?)\\s*$")
  bad_components <- components[is.na(component_match[, 1, drop = TRUE])]
  if(length(bad_components) > 0) stop("Bad components in reaction: ", paste(bad_components, collapse = ", "))
  # extract counts
  counts <- component_match[, 2, drop = TRUE]
  counts <- ifelse(counts == "", "1", counts)
  counts <- as.numeric(counts)

  # parse molecules
  mols <- as_mol(component_match[, 3, drop = TRUE], validate = validate)

  # return mols and counts
  list(mols = mols, counts = counts)
}

count_names <- function(x, env = new.env(emptyenv())) {
  if(is.atomic(x)) {
    stop("Misplaced atomic in formula")
  } else if(is.name(x)) {
    count_names_add(x, env)
  } else if(is.call(x)) {
    if(identical(x[[1]], quote(`+`))) {
      for(i in setdiff(seq_along(x), 1)) {
        count_names(x[[i]], env)
      }
    } else if(identical(x[[1]], quote(`*`))) {
      # * should have one atomic and one name
      if(length(x) != 3) stop("Invalid syntax around `*`: too many arguments")
      if(is.atomic(x[[2]])) {
        if(is.name(x[[3]])) {
          count_names_add(x[[3]], env, n = x[[2]])
        } else {
          stop("Invalid `*` syntax: need one name and one coefficient")
        }
      } else if(is.atomic(x[[3]])) {
        if(is.name(x[[2]])) {
          count_names_add(x[[2]], env, n = x[[3]])
        } else {
          stop("Invalid `*` syntax: need one name and one coefficient")
        }
      } else {
        stop("Invalid `*` syntax: need one name and one coefficient")
      }
    } else {
      stop("Invalid operator in formula: ", x[[1]])
    }
  } else if (is.pairlist(x)) {
    stop("Misplaced pairlist in formula")
  } else {
    stop("Don't know how to handle type ", typeof(x),
         call. = FALSE)
  }
  # return env as a named vector
  unlist(as.list(env))
}

count_names_add <- function(name, env, n=1L) {
  name <- as.character(name)
  if(!is.numeric(n)) stop("Invalid coefficient: ", n)
  if(name %in% names(env)) {
    env[[name]] <- env[[name]] + as.integer(n)
  } else {
    env[[name]] <- as.integer(n)
  }
}

#' Create, validate reaction objects
#'
#' @param x A reaction object
#' @param log_k,value The base 10 log of the reaction quotient
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' r <- new_reaction(list(mol = as_mol(c("H2O", "H+", "H3O+")),
#'                        coefficient = c(1, 1, -1)))
#' validate_reaction(r)
#'
new_reaction <- function(x, log_k = NA_real_) {
  # check base type
  if(!is.list(x)) stop("x must be a list")
  # return structure
  structure(x, class = "reaction", log_k = log_k)
}

#' @rdname new_reaction
#' @export
validate_reaction <- function(x) {
  # check type
  if(!is.list(x)) stop("x must be a list")
  if(!inherits(x, "reaction")) stop("x must inherit 'reaction'")

  # check required names
  if(!("mol" %in% names(x))) stop("Required component 'mol' missing from x")
  if(!("coefficient" %in% names(x))) stop("Required component 'coefficient' missing from x")
  # check types
  if(!is_mol(x$mol)) stop("x$mol is not a mol vector")
  if(!is.numeric(x$coefficient)) stop("x$coefficient is not numeric")
  # check lengths
  if(length(x$mol) != length(x$coefficient)) stop("length(x$mol) != length(x$coefficient)")
  # check log_k
  if(!is.double(attr(x, "log_k"))) stop('attr(x, "log_k") is not a double')

  # return x, invisibly
  invisible(x)
}

#' @rdname new_reaction
#' @export
#' @importFrom purrr %||%
log_k <- function(x) {
  if(inherits(x, "reaction_list")) {
    new_reaction_list(lapply(x, log_k))
  } else if(inherits(x, "reaction")) {
    attr(x, "log_k") %||% NA_real_
  } else {
    log_k(as_reaction(x))
  }
}

#' @rdname new_reaction
#' @export
`log_k<-` <- function(x, value) {
  if(inherits(x, "reaction_list")) {
    new_reaction_list(purrr::map2(x, value, `log_k<-`))
  } else if(inherits(x, "reaction")) {
    attr(x, "log_k") <- as.numeric(value)
  } else {
    stop("Can't assign log_k to object of class ", paste(class(x), collapse = "/"))
  }

  x
}

#' @rdname new_reaction
#' @export
is_reaction <- function(x) {
  inherits(x, "reaction")
}

#' @rdname new_reaction
#' @export
new_reaction_list <- function(x) {
  # check base type
  if(!is.list(x)) stop("x must be a list")
  # return structure
  structure(x, class = "reaction_list")
}

#' @rdname new_reaction
#' @export
is_reaction_list <- function(x) {
  inherits(x, "reaction_list")
}

#' @rdname new_reaction
#' @export
validate_reaction_list <- function(x) {
  # check type
  if(!is.list(x)) stop("x must be a list")
  if(!inherits(x, "reaction_list")) stop("x must inherit 'reaction_list'")

  validation <- lapply(x, function(x) try(validate_reaction(x), silent = TRUE))
  val_error <- vapply(validation, inherits, "try-error", FUN.VALUE = logical(1))
  if(any(val_error)) stop("reaction at positions ",
                          paste(which(val_error), collapse = ", "),
                          " are invalid reaction objects")

  # return x, invisibly
  invisible(x)
}

#' Subset, combine reaction objects
#'
#' @param x A reaction object
#' @param i The index
#' @param ... Ignored
#'
#' @return A reaction object
#' @export
#' @rdname reactionsubset
#'
#'
#' @examples
#' r <- as_reaction("2H2 + O2 = 2H2O")
#' lhs(r)
#' rhs(r)
#'
`[.reaction` <- function(x, i, ...) {
  new_reaction(list(mol = x$mol[i], coefficient = x$coefficient[i]))
}

#' @export
#' @rdname reactionsubset
`[.reaction_list` <- function(x, i, ...) {
  new_reaction_list(unclass(x)[i])
}

#' @export
#' @rdname reactionsubset
lhs <- function(x) {
  if(is_reaction_list(x)) return(new_reaction_list(lapply(x, lhs)))
  x[x$coefficient >= 0]
}

#' @export
#' @rdname reactionsubset
rhs <- function(x) {
  if(is_reaction_list(x)) return(new_reaction_list(lapply(x, rhs)))
  -(x[x$coefficient < 0])
}

#' Coerce reactions to a character vector
#'
#' @param x A reaction object
#' @param equals_sign,mol_sep An equals/plus sign to use for the reaction
#' @param wrap_coeff Wrap coefficients for fancy formatting
#' @param wrap_super,wrap_sub Wrap super/subscript for fancy formatting
#' @param ... Ignored
#'
#' @return A character vector
#' @export
#'
#' @examples
#' r <- as_reaction("2H2 + O2 = 2H2O")
#' as.character(r)
#' print(r)
#'
as.character.reaction <- function(x, equals_sign = "=", wrap_coeff = identity,
                                  wrap_super = identity, wrap_sub = identity,
                                  mol_sep = " + ", ...) {
  sides <- vapply(list(lhs(x), rhs(x)), function(r) {
    coeffs <- abs(r$coefficient)
    coeffs <- ifelse(abs(coeffs - 1) < .Machine$double.eps^0.5,
                     "",
                     wrap_coeff(format(coeffs, ...)))
    paste0(coeffs, as.character(r$mol, wrap_super = wrap_super, wrap_sub = wrap_sub, ...), collapse = mol_sep)
  }, character(1))
  paste(sides, collapse = paste0(" ", equals_sign, " "))
}

#' @rdname as.character.reaction
#' @export
as.character.reaction_list <- function(x, ...) {
  vapply(x, as.character, character(1), ...)
}

#' @rdname as.character.reaction
#' @export
print.reaction <- function(x, ...) {
  cat("<reaction>", as.character(x, ...), "\tlog_k =", log_k(x))
  invisible(x)
}

#' @rdname as.character.reaction
#' @export
print.reaction_list <- function(x, ...) {
  cat("<reaction_list>\n")
  for(reaction in x) {
    print(reaction, ...)
    cat("\n")
  }
  invisible(x)
}

#' Reaction arithmetic
#'
#' @param x A reaction object
#' @param y A reaction object or operand
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' r1 <- as_reaction("O2 + 2H2 = 2H2O", log_k = -46.62)
#' r2 <- as_reaction("H3O+ + OH- = 2H2O", log_k = 14)
#'
#' -r1
#' -r2
#' r1 + r2
#' r1 - r2
#' r1 * 2
#' 2 * r1
#' r1 * -2
#'
`*.reaction` <- function(x, y) {
  if(is_reaction(x) && is.numeric(y)) {
    x$coefficient <- x$coefficient * y
    log_k(x) <- log_k(x) * y
    x
  } else if(is_reaction(y) && is.numeric(x)) {
    y$coefficient <- y$coefficient * x
    log_k(y) <- log_k(y) * x
    y
  } else {
    stop("* operator not defined for types ", class(x)[1], ", ", class(y)[1])
  }
}

#' @export
`/.reaction` <- function(x, y) {
  if(is_reaction(x) && is.numeric(y)) {
    x$coefficient <- x$coefficient / y
    log_k(x) <- log_k(x) / y
    x
  } else {
    stop("/ operator not defined for types ", class(x)[1], ", ", class(y)[1])
  }
}

#' @export
`+.reaction` <- function(x, y) {
  if(missing(y)) return(x)
  x <- as_reaction(x)
  y <- as_reaction(y)
  new_reaction(
    list(
      mol = c(x$mol, y$mol),
      coefficient = c(x$coefficient, y$coefficient)
    ),
    log_k = log_k(x) + log_k(y)
  )
}

#' @export
`-.reaction` <- function(x, y) {
  if(missing(y)) return(x * -1)
  x <- as_reaction(x)
  y <- as_reaction(y)
  x + (y * -1)
}

#' @export
`*.reaction_list` <- function(x, y) {
  if((is_reaction_list(x) && is.numeric(y)) || (is_reaction_list(y) && is.numeric(x))) {
    xy <- vctrs::vec_recycle_common(unclass(x), unclass(y))
    new_reaction_list(purrr::pmap(xy, `*.reaction`))
  } else {
    stop("* operator not defined for types ", class(x)[1],
         ", ", class(y)[1])
  }
}

#' @export
`/.reaction_list` <- function(x, y) {
  if(is_reaction_list(x) && is.numeric(y)) {
    xy <- vctrs::vec_recycle_common(unclass(x), unclass(y))
    new_reaction_list(purrr::pmap(xy, `/.reaction`))
  } else {
    stop("/ operator not defined for types ", class(x)[1],
         ", ", class(y)[1])
  }
}

#' @export
`+.reaction_list` <- function(x, y) {
  if(missing(y)) return(x)
  x <- as_reaction_list(x)
  y <- as_reaction_list(y)
  xy <- vctrs::vec_recycle_common(unclass(x), unclass(y))
  new_reaction_list(purrr::pmap(xy, `+.reaction`))
}

#' @export
`-.reaction_list` <- function(x, y) {
  if(missing(y)) return(x * -1)
  x <- as_reaction_list(x)
  y <- as_reaction_list(y)
  x + (y * -1)
}

#' Simplify reaction objects
#'
#' @param x A reaction object
#' @param ... ignored
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' r <- as_reaction("NH3 + H+ + H2O = H3O+ + NH3")
#' simplify_reaction(r)
#'
#' r2 <- as_reaction("0NH3 + H+ + H2O = H3O+")
#' remove_zero_counts(r2)
#'
simplify_reaction <- function(x, ...) {
  UseMethod("simplify_reaction")
}

#' @rdname simplify_reaction
#' @export
simplify_reaction.reaction <- function(x, ...) {
  unique_mols <- unique(x$mol)
  unique_coefficient <- vapply(unique_mols, function(m) {
    sum(x$coefficient[x$mol == as_mol(m)])
  }, numeric(1))
  new_reaction(list(mol = unique_mols, coefficient = unique_coefficient), log_k = log_k(x))
}

#' @rdname simplify_reaction
#' @export
remove_zero_counts.reaction <- function(x, ...) {
  y <- x[!is.na(x$mol) & (x$coefficient != 0)]
  log_k(y) <- log_k(x)
  y
}

#' @rdname simplify_reaction
#' @export
remove_zero_counts.reaction_list <- function(x, ...) {
  new_reaction_list(lapply(x, remove_zero_counts.reaction))
}

#' @rdname simplify_reaction
#' @export
simplify_reaction.reaction_list <- function(x, ...) {
  new_reaction_list(lapply(x, simplify_reaction.reaction))
}

#' Data frame reprsentation of a reaction object
#'
#' @param x A reaction object
#' @param ... Ignored
#'
#' @return A tibble with one row per molecule in reaction
#' @export
#' @rdname reaction_tibble
#'
#' @importFrom tibble as_tibble
#'
#' @examples
#' as.data.frame(as_reaction("O2 + 2H2 = 2 H2O"))
#' library(tibble)
#' as_tibble(as_reaction("O2 + 2H2 = 2 H2O"))
#'
#' as.matrix(as_reaction("O2 + 2H2 = 2 H2O"))
#'
as.data.frame.reaction <- function(x, ...) {
  as.data.frame(as_tibble.reaction(x, ...))
}

#' @rdname reaction_tibble
#' @export
as_tibble.reaction <- function(x, ...) {
  tbl_mol <- as_tibble.mol(x$mol)
  tbl_mol$coefficient <- x$coefficient
  if(nrow(tbl_mol) > 0) {
    tbl_mol$reaction <- as.character(x)
    tbl_mol$log_k <- log_k(x)
  } else {
    tbl_mol$reaction <- character(0)
    tbl_mol$log_k <- numeric(0)
  }
  # reorder columns
  cols_first <- c("reaction", "log_k", "mol", "coefficient", "charge", "mass")
  tbl_mol[c(cols_first, setdiff(names(tbl_mol), cols_first))]
}

#' @rdname reaction_tibble
#' @export
as.matrix.reaction <- function(x, ...) {
  as.matrix(x$mol) * x$coefficient
}

#' @rdname reaction_tibble
#' @export
as_tibble.reaction_list <- function(x, ...) {
  if(is.null(names(x))) {
    names(x) <- as.character(seq_along(x))
  }

  purrr::map_df(x, as_tibble.reaction, .id = "which")
}

#' @rdname reaction_tibble
#' @export
as.data.frame.reaction_list <- function(x, ...) {
  as.data.frame(as_tibble.reaction_list(x, ...))
}

#' Balance reactions
#'
#' @param x A reaction object
#' @param charge Balance charge?
#' @param tol Tolerance for zero detection
#'
#' @return A balanced reaction or logical describing if the reaction is balanced
#' @export
#'
#' @examples
#' is_balanced("O2 + 2H2 = 2H2O")
#' is_balanced("O2 + H2 = 2H2O")
#' is_balanced("O2-4 + 2H2 = 2H2O", charge = FALSE)
#' is_balanced("O2-4 + 2H2 = 2H2O", charge = TRUE)
#'
is_balanced <- function(x, charge = TRUE, tol = .Machine$double.eps^0.5) {
  if(is_reaction_list(x)) {
    return(vapply(x, is_balanced, charge = charge, tol = tol, FUN.VALUE = logical(1)))
  }

  x <- as_reaction(x)
  all_mols <- x$mol * x$coefficient
  result <- remove_zero_counts(simplify_mol(do.call(combine_molecules, all_mols)), tol = tol)
  if(charge) {
    (length(result) == 0) && (abs(charge(result)) <= tol)
  } else {
    length(result) == 0
  }
}

#' @rdname is_balanced
#' @export
balance <- function(x, charge = TRUE, tol = .Machine$double.eps^0.5) {
  if(is_reaction_list(x)) {
    return(new_reaction_list(lapply(x, balance, charge = charge, tol = tol)))
  }

  x <- as_reaction(x)
  # check balance
  if(is_balanced(x, charge = charge)) return(x)

  # solve the null space of the (molecule) matrix
  xmat <- as.matrix(x$mol)
  null_space <- MASS::Null(xmat)

  # check dimensions (need dimensions of length(x$mol), 1L)
  if(identical(dim(null_space), c(length(x$mol), 1L))) {
    ns <- null_space[, 1, drop = TRUE]
    # check for all zero null space
    if(any(abs(ns) < tol)) stop("Could not balance reaction: null space contains zero values ",
                                "for mols ", paste(x$mol[abs(ns) < tol], collapse = ", "))
    # scale such that min(ns) == 1, first value is positive
    ns <- ns / min(abs(ns))
    x$coefficient <- ns * sign(ns[1])
    # check balance
    if(!is_balanced(x, charge = FALSE)) {
      stop("Could not balance reaction: coefficients did not result in ",
           "an element-balanced reaction: ",
           paste(x$coefficient, as.character(x$mol), collapse = " // "))
    }
    if(charge && !is_balanced(x, charge = TRUE)) {
      # add electron to the correct side
      total_charge <- charge(x)
      x$mol <- c(x$mol, as_mol(electron_))
      x$coefficient <- c(x$coefficient, total_charge)
      # combines all electrons into one component
      simplify_reaction(x)
    }
    # return reaction
    x
  } else {
    stop("Could not balance reaction: null space has incorrect dimensions: ",
         paste(dim(null_space), collapse = ", "))
  }
}

