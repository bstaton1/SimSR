#' Initialize Simulation
#'
#' Handles setting up the dimensional variables
#' and drawing random quantities for driving the
#' operating model and observation.
#'
#' @param nt Numeric vector of length 1: the number of
#'   calendar years with observations to simulate
#' @param a_min Numeric vector of length 1: youngest age at maturity
#' @param a_max Numeric vector of length 1: oldest age at maturity
#' @param U_msy Numeric vector storing U_msy for each stock to simulate
#' @param S_msy Same as \code{S_msy}, except S_msy. Must be the same
#'   length as \code{U_msy}
#' @param phi Numeric vector of length 1: the autocorrelation coefficient for recruitment residuals.
#'   Set to zero to turn off the AR(1) process
#' @param pi_grand Average maturity schedule across all stocks and years. Must have the same length as
#'   suggested by the number of ages of maturity calculated from supplied arguments \code{a_min} and \code{a_max}.
#'   Will be scaled to sum to 1 if does not already.
#' @param U_SUM Numeric vector of length 1: beta sample size of implementation error
#' @param min_S_cv Numeric vector of length 1: smallest observation CV for any substocks's escapement.
#' @param max_S_cv Numeric vector of length 1: largest observation CV for any substocks's escapement.
#' @param min_C_cv Numeric vector of length 1: smallest observation CV for any year's total harvest.
#' @param max_C_cv Numeric vector of length 1: largest observation CV for any year's total harvest.
#' @param x_ESS Numeric vector of length 1: the effective sample size for scale sampling on stocks that have age data.
#' @return List object with elements needed to drive many of the other functions in this package
#'   via the \code{params} argument.

#' @export

init_sim = function(nt = 42, a_min = 4, a_max = 7, U_msy, S_msy, U_SUM = 100, phi = 0.3, pi_grand = c(0.2, 0.4, 0.37, 0.03),
                    min_S_cv = 0.1, max_S_cv = 0.2, min_C_cv = 0.1, max_C_cv = 0.2, x_ESS = 100) {

  if (length(U_msy) != length(S_msy)) {
    stop ("U_msy and S_msy must have the same length")
  }

  # dimensions
  ns = length(U_msy)
  na = a_max - a_min + 1
  ages = a_min:a_max
  ny = nt + na - 1

  if (length(pi_grand) != na) {
    stop ("pi_grand supplied with a different number of elements than suggested by a_min and a_max, not allowed")
  }
  # ensure pi_grand sums to 1
  pi_grand = pi_grand/sum(pi_grand)

  # obtain leading parameters on other scale
  alpha = exp(U_msy)/(1 - U_msy)
  beta = U_msy/S_msy
  log_alpha = log(alpha)

  # fishery parameters
  v = rep(1, ns)

  # create covariance matrix
  Rwish = matrix(rep(-0.01, ns^2), ns, ns)
  diag(Rwish) = rep(0.22, ns)
  dfwish = 35
  Sigma = solve(rWishart(1, df = dfwish, Sigma = Rwish)[,,1])
  sigma = sqrt(diag(Sigma))

  # obtain correlation matrix
  rho_mat = matrix(NA, ns, ns)
  for (i in 1:ns) {
    for (j in 1:ns) {
      rho_mat[i,j] = Sigma[i,j]/(sigma[i] * sigma[j])
    }
  }

  # OBSERVATION VARIABILITY
  cv_S_ts_obs = matrix(runif(ns, min_S_cv, max_S_cv), nt, ns, byrow = T)
  cv_C_t_obs = rep(runif(1, min_C_cv, max_C_cv), nt)

  out = list(
    ns = ns,
    nt = nt,
    a_min = a_min,
    a_max = a_max,
    na = na,
    ages = ages,
    ny = ny,
    U_msy = U_msy,
    S_msy = S_msy,
    alpha = alpha,
    beta = beta,
    log_alpha = log_alpha,
    phi = phi,
    U_SUM = U_SUM,
    v = v,
    sigma = sigma,
    Sigma = Sigma,
    rho_mat = rho_mat,
    cv_S_ts_obs = cv_S_ts_obs,
    cv_C_t_obs = cv_C_t_obs,
    pi_grand = pi_grand,
    x_ESS = x_ESS
  )

  return(out)
}
