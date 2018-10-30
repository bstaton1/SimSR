#' Generate Management Errors
#'
#' Generate a beta random deviate around a target exploitation rate
#' with variability given by \code{SUM}, the sum of the two beta shape parameters.
#'
#' @param U_target Numeric vector: the mean target exploitation rate
#' @param SUM Numeric vector: the sum of the two shape parameters of the beta distribution.
#'   The larger this value is, the less implementation error there is.

implement_error = function(U_target, SUM) {
  Ua = U_target * SUM
  Ub = SUM - Ua

  rbeta(1, Ua, Ub)
}


# U_target = runif(1000, 0.01, 0.99)
# U_real = sapply(U_target, function(x) implement_error(U_target = x, SUM = 10))
#
# windows()
# plot(U_target ~ U_real)
