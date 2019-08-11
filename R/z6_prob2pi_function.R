#' Maturity Conversion Function 1
#'
#' Converts a vector from the marginal probability of maturity at each age
#' to a vector containing the conditional probabilities of maturity at ages
#' \code{a:(a_max-1)} given not yet matured.

pi2prob = function(pi) {
  A = length(pi)
  prob = numeric(A-1)
  prob[1] = pi[1]

  for (a in 2:(A-1)) {
    prob[a] = pi[a]/(1 - sum(pi[1:(a-1)]))
  }
  prob
}

#' Maturity Conversion Function 2
#'
#' Reverses the action of \code{\link{pi2prob}}

prob2pi = function(prob) {
  A = length(prob) + 1
  pi = numeric(A)
  pi[1] = prob[1]
  for (a in 2:(A-1)) { pi[a] = prob[a] * (1 - sum(pi[1:(a-1)]))}
  pi[A] = 1 - sum(pi[1:(A - 1)])

  pi
}

#' Maturity Conversion Function 3
#'
#' Converts from \code{pi} to a set of logit-scale coefficients that is used to obtain
#' the \code{prob} vector

pi2beta = function(pi) {
  prob = pi2prob(pi)
  lprob = StatonMisc::logit(prob)
  A = length(pi)

  b = numeric(A-1)
  b[1] = lprob[1]
  for (a in 2:(A-1)) b[a] = lprob[a] - b[1]

  # m = matrix(0, A-1, A-1)
  # m[,1] = 1
  # for (a in 2:(A-1)) m[a,a] = 1
  # prob2pi(StatonMisc::expit(m %*% b))

  b
}
