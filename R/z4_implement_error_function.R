#' Generate Management Errors
#'
#' Generate a beta random deviate around a target exploitation rate
#' with variability given by \code{SUM}, the sum of the two beta shape parameters.
#'
#' @param U_target Numeric vector of length 1: the mean target exploitation rate
#' @param SUM Numeric vector of length 1: the sum of the two shape parameters of the beta distribution.
#'   The larger this value is, the less implementation error there is.

implement_error = function(U_target, SUM) {
  Ua = U_target * SUM
  Ub = SUM - Ua

  rbeta(1, Ua, Ub)
}
