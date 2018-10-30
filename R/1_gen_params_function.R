#' Initialize Simulation
#'
#' Handles setting up the dimensional variables
#' and drawing random quantities for driving the
#' operating model and observation
#'
#' @param nt Numeric vector of length 1: the number of
#'   calendar years with observations to simulate
#' @param a_min Numeric vector of length 1: youngest age at maturity
#' @param a_max Numeric vector of length 1: oldest age at maturity
#' @param U_msy Numeric vector storing U_msy for each stock to simulate
#' @param S_msy Same as \code{S_msy}, except S_msy. Must be the same
#'   length as \code{U_msy}
#' @param rho Numeric vector of length 1: the average correlation between
#'   recruitment residuals among substocks
#' @param min_sigR Numeric vector of length 1: the smallest lognormal SD that
#'   can be drawn for a substock's recruitment residuals
#' @param max_sigR Similar to \code{min_sigR}
#'
#' @export

init_sim = function(nt = 42, a_min = 4, a_max = 7, U_msy, S_msy, rho = 0.5, min_sigR = 0.4, max_sigR = 0.6) {

  if (length(U_msy) != length(S_msy)) {
    stop ("U_msy and S_msy must have the same length")
  }

  # dimensions
  ns = length(U_msy)
  na = a_max - a_min + 1
  ages = a_min:a_max
  ny = nt + na - 1

  # obtain leading parameters on other scale
  alpha = exp(U_msy)/(1 - U_msy)
  beta = U_msy/S_msy
  log_alpha = log(alpha)
  phi = 0.7

  # fishery parameters
  max_p_overfished = 0.1
  U_SUM = 50
  v = rep(1, ns)

  # create covariance matrix
  sigma = runif(ns, min_sigR, max_sigR)
  Sigma = matrix(1, ns, ns)
  rho_mat = matrix(NA, ns, ns)
  for (i in 1:ns) {
    for (j in 1:ns) {
      Sigma[i,j] = sigma[i] * sigma[j] * rho
    }
  }
  diag(Sigma) = sigma^2

  # obtain correlation matrix
  for (i in 1:ns) {
    for (j in 1:ns) {
      rho_mat[i,j] = Sigma[i,j]/(sigma[i] * sigma[j])
    }
  }

  # OBSERVATION VARIABILITY
  cv_S_ts_obs = matrix(runif(ns, 0.1, 0.2), nt, ns, byrow = T)
  cv_C_t_obs = rep(runif(1, 0.1, 0.2), nt)

  x_ESS = 100

  out = list(
    ns = ns,
    nt = nt,
    a_min = a_min,
    a_max = a_max,
    na = na,
    ages = ages,
    ny = ny,
    pi = pi,
    U_msy = U_msy,
    S_msy = S_msy,
    alpha = alpha,
    beta = beta,
    log_alpha = log_alpha,
    phi = phi,
    max_p_overfished = max_p_overfished,
    U_SUM = U_SUM,
    v = v,
    sigma = sigma,
    rho = rho,
    Sigma = Sigma,
    rho_mat = rho_mat,
    cv_S_ts_obs = cv_S_ts_obs,
    cv_C_t_obs = cv_C_t_obs,
    x_ESS = x_ESS
  )

  return(out)
}
