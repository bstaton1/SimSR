
#' Convert Beta Deviates to Dirichlet Deviates
#'

pi2prob = function(pi) {
  A = length(pi)
  prob = numeric(A-1)
  prob[1] = pi[1]

  for (a in 2:(A-1)) {
    prob[a] = pi[a]/(1 - sum(pi[1:(a-1)]))
  }
  prob
}

prob2pi = function(prob) {
  A = length(prob) + 1
  pi = numeric(A)
  pi[1] = prob[1]
  for (a in 2:(A-1)) { pi[a] = prob[a] * (1 - sum(pi[1:(a-1)]))}
  pi[A] = 1 - sum(pi[1:(A - 1)])

  pi
}

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
