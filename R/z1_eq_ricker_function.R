#' Calculate Equilibrium Quantities from a Ricker Function
#'
#' Given known (or estimated) values of alpha, beta,
#'   U_msy, and S_msy, for a single stock,
#'   calculate Seq and Ceq at each of a set of exploitation rates
#'
#' @param alpha a numeric vector of length 1.
#' @param beta a numeric vector of length 1.
#' @param U_msy a numeric vector of length 1.
#' @param S_msy a numeric vector of length 1.
#' @param U_range a numeric vector storing the exploitation
#'   rates to calculate the equilibrium quantities at.

eq_ricker = function(alpha, beta, U_msy, S_msy, U_range) {

  # calculate Seq
  Seq = 1/beta * (log(alpha) + log(1 - U_range))
  Seq[Seq < 0] = 0

  # calculate Ceq
  Ceq = Seq * U_range/(1 - U_range)
  Ceq[is.na(Ceq)] = 0

  # determine if U would overfish it or push it to extinction
  overfished = ifelse(U_range > U_msy, 1, 0)
  extinct = ifelse(Seq == 0, 1, 0)

  # output
  return(list(S = Seq, C = Ceq, overfished = overfished, extinct = extinct))
}
