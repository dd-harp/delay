# -------------------------------------------------------------------------------- #
#
#   Discretize continuous-time output
#   Sean Wu (slwu89@berkeley.edu)
#   January 2020
#
# -------------------------------------------------------------------------------- #

#' Discretize Output
#'
#' Modified from \code{smfsb} package
#'
#' @param out a matrix
#' @param out dt the discretization lattice
#'
#' @useDynLib delay discretize_C
#' @export
discretize <- function(out, dt=1){
  storage.mode(out) <- "double"
  .Call(discretize_C,as.matrix(out),as.numeric(dt))
}
