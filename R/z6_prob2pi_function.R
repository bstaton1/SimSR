#' inverse logit
#'
#' @param x the value on the logit scale to convert to the natural scale.
expit = function(x) {
  exp(x)/(1 + exp(x))
}

#' Convert Beta Deviates to Dirichlet Deviates
#'

prob2pi = function(prob) {
  A = length(prob) + 1
  pi = numeric(A)
  pi[1] = prob[1]
  for (a in 2:(A-1)) { pi[a] = prob[a] * (1 - sum(pi[1:(a-1)]))}
  pi[A] = 1 - sum(pi[1:(A - 1)])

  pi
}
