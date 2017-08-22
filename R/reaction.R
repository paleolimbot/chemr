
#' Create reaction objects
#'
#' @param lhs The lefthand side
#' @param counts_lhs The counts for molecules on the lefthand side
#' @param rhs The righthand side
#' @param counts_rhs The counts for molecules on the righthand side
#' @param validate Flag to validate molecules in the reaction
#' @param x An object to convert to a reaction object
#' @param coefficients Coefficients corresponding to as_reaction input
#' @param ... Passed to/from methods
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' reaction(~H2O + `H+`, ~`H3O+`)
#' reaction(c("H2O", "H+"), "H3O+")
#'
reaction <- function(lhs, rhs, counts_lhs = 1, counts_rhs = 1, validate = TRUE) {
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
  r <- new_reaction(list(mol = c(lhs, rhs),
                         coefficients = c(counts_lhs, counts_rhs)))
  # need to validate reaction counts regardless of validate arg
  validate_reaction(r)
  # return r
  r
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
as_reaction.mol <- function(x, coefficients, validate = TRUE, ...) {
  if(length(coefficients) != length(x)) {
    stop("length(coefficients) is not equal to length(x)")
  }
  r <- new_reaction(list(mol = x, coefficients = coefficients))
  if(validate) validate_reaction(r)
  r
}

#' @rdname reaction
#' @export
as_reaction.formula <- function(x, validate = TRUE, ...) {

  counts_lhs <- count_names(x[[2]])
  counts_rhs <- count_names(x[[3]])

  mols_lhs <- as_mol(names(counts_lhs), validate = validate, ...)
  mols_rhs <- as_mol(names(counts_rhs), validate = validate, ...)

  reaction(lhs = mols_lhs, rhs = mols_rhs,
           counts_lhs =  stats::setNames(counts_lhs, NULL),
           counts_rhs = stats::setNames(counts_rhs, NULL),
           validate = FALSE)
}

#' @rdname reaction
#' @export
as_reaction.character <- function(x, validate = TRUE, ...) {
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
           validate = FALSE)

}

parse_side <- function(side, validate = TRUE) {
  components <- stringr::str_split(side, "\\s+\\+\\s+")[[1]]
  # component is a count plus a molecule id
  component_match <- stringr::str_match(components, "^\\s*([0-9.]*)(.*?)\\s*$")
  bad_components <- components[is.na(component_match[, 1, drop = TRUE])]
  if(length(bad_components) > 0) stop("Bad components in reaction: ",
                                      paste(bad_components, collapse = ", "))
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
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' r <- new_reaction(list(mol = as_mol(c("H2O", "H+", "H3O+")),
#'                        coefficients = c(1, 1, -1)))
#' validate_reaction(r)
#'
new_reaction <- function(x) {
  # check base type
  if(!is.list(x)) stop("x must be a list")
  # return structure
  structure(x, class = "reaction")
}

#' @rdname new_reaction
#' @export
validate_reaction <- function(x) {
  # check type
  if(!is.list(x)) stop("x must be a list")
  if(!inherits(x, "reaction")) stop("x must inherit 'reaction'")

  # check required names
  if(!("mol" %in% names(x))) stop("Required component 'mol' missing from x")
  if(!("coefficients" %in% names(x))) stop("Required component 'coefficients' missing from x")
  # check types
  if(!is_mol(x$mol)) stop("x$mol is not a mol vector")
  if(!is.numeric(x$coefficients)) stop("x$coefficients is not numeric")
  # check lengths
  if(length(x$mol) != length(x$coefficients)) stop("length(x$mol) != length(x$coefficients)")

  # return x, invisibly
  invisible(x)
}

#' @rdname new_reaction
#' @export
is_reaction <- function(x) {
  inherits(x, "reaction")
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
  new_reaction(list(mol = x$mol[i], coefficients = x$coefficients[i]))
}

#' @export
#' @rdname reactionsubset
lhs <- function(x) {
  x[x$coefficients >= 0]
}

#' @export
#' @rdname reactionsubset
rhs <- function(x) {
  x[x$coefficients < 0]
}

#' Coerce reactions to a character vector
#'
#' @param x A reaction object
#' @param ... Ignored
#'
#' @return A character vector
#' @export
#'
#' @examples
#' r <- as_reaction(2*H2 + O2 ~ 2*H2O)
#' as.character(r)
#' print(r)
#'
as.character.reaction <- function(x, ...) {
  sides <- vapply(list(lhs(x), rhs(x)), function(r) {
    coeffs <- abs(r$coefficients)
    coeffs <- ifelse(coeffs == 1, "", as.character(coeffs))
    paste0(coeffs, as.character(r$mol), collapse = " + ")
  }, character(1))
  paste(sides, collapse = " = ")
}

#' @rdname as.character.reaction
#' @export
print.reaction <- function(x, ...) {
  cat("<reaction>", as.character(x))
  invisible(x)
}
