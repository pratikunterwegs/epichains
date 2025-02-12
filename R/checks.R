#' Check if offspring argument is specified as a character string
#'
#' @param offspring_dist Offspring distribution: a character string
#' corresponding to the R distribution function (e.g., "pois" for Poisson,
#' where \code{\link{rpois}} is the R function to generate Poisson random
#' numbers).
#' @keywords internal
check_offspring_valid <- function(offspring_dist) {
  if (!checkmate::test_string(offspring_dist)) {
    stop(sprintf(
      "%s %s",
      "'offspring_dist' must be specified as a character string.",
      "Did you forget to enclose it in quotes?"
    ))
  }
}


#' Check if constructed random number generator for offspring exists
#'
#' @param roffspring_name Constructed random offspring sampler: a character
#' string corresponding to the R distribution function (e.g., "rpois" for
#' Poisson.
#' @keywords internal
check_offspring_func_valid <- function(roffspring_name) {
  if (!(exists(roffspring_name)) ||
        !checkmate::test_function(get(roffspring_name))) {
    stop("Function ", roffspring_name, " does not exist.")
  }
}


#' Check if the serials_dist argument is specified as a function
#'
#' @param serials_dist The serial interval distribution function; the name of a
#' user-defined named or anonymous function with only one argument `n`,
#' representing the number of serial intervals to generate.
#'
#' @keywords internal
check_serial_valid <- function(serials_dist) {
  if (!checkmate::test_function(serials_dist)) {
    stop(sprintf(
      "%s %s",
      "The `serials_dist` argument must be a function",
      "(see details in ?sim_chain_tree)."
    ))
  }
}


#' Check that nchains is greater than 0 and not infinity
#'
#' @param nchains Number of chains to simulate.
#'
#' @keywords internal
check_nchains_valid <- function(nchains) {
  if (!checkmate::test_count(nchains, positive = TRUE)) {
    stop("`nchains` must be > 0 but less than `Inf`")
  }
}
