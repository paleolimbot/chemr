
#' Create reaction objects
#'
#' @param lhs The lefthand side
#' @param counts_lhs The counts for molecules on the lefthand side
#' @param rhs The righthand side
#' @param counts_rhs The counts for molecules on the righthand side
#' @param validate Flag to validate molecules in the reaction
#' @param x An object to convert to a reaction object
#' @param ... Passed to/from methods
#'
#' @return A reaction object
#' @export
#'
#' @examples
#' reaction(~H2O + `H+`, ~`H3O+`)
#' reaction(c("H2O", "H+"), "H3O+")
#'
reaction <- function(lhs, rhs, counts_lhs = rep(1, length(lhs)),
                     counts_rhs = rep(1, length(rhs)), validate = TRUE) {
  lhs <- as_mol(lhs, validate = validate)
  rhs <- as_mol(rhs, validate = validate)
  r <- new_reaction(list(lhs = lhs,
                         rhs = rhs,
                         counts_lhs = counts_lhs,
                         counts_rhs = counts_rhs))
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
#' r <- new_reaction(list(lhs = as_mol(c("H2O", "H+")), rhs = as_mol("H3O+"),
#'                        counts_lhs = c(1, 1), counts_rhs = 1))
#' validate_reaction(r)
#'
new_reaction <- function(x) {
  if(!is.list(x)) stop("x must be a list")
  structure(x, class = "reaction")
}

#' @rdname new_reaction
#' @export
validate_reaction <- function(x) {
  # check type
  if(!is.list(x)) stop("x must be a list")
  if(!inherits(x, "reaction")) stop("x must inherit 'reaction'")
  # check lhs, rhs
  if(!("lhs" %in% names(x))) stop("x missing element 'lhs'")
  if(!("rhs" %in% names(x))) stop("x missing element 'rhs'")
  if(!is_mol(x$lhs)) stop("x$lhs does not inherit 'mol'")
  if(!is_mol(x$rhs)) stop("x$rhs does not inherit 'mol'")
  # check counts types
  if(!("counts_lhs" %in% names(x))) stop("x missing element 'lhs'")
  if(!("counts_rhs" %in% names(x))) stop("x missing element 'rhs'")
  if(!is.numeric(x$counts_lhs)) stop("x$counts_lhs is not numeric")
  if(!is.numeric(x$counts_rhs)) stop("x$counts_rhs is not numeric")
  # check counts lengths
  if(length(x$counts_lhs) != length(x$lhs)) stop("x$counts_lhs must have same length as x$lhs")
  if(length(x$counts_rhs) != length(x$rhs)) stop("x$counts_rhs must have same length as x$rhs")

  # return x, invisibly
  invisible(x)
}

#' @rdname new_reaction
#' @export
is_reaction <- function(x) {
  inherits(x, "reaction")
}
