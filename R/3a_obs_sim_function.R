#' Generate Observed States
#'
#' Based on the true states obtained using \code{ricker_sim}, introduce observation errors
#'
#' @param params A list created using \code{init_sim()}.
#' @param true A list created using \code{ricker_sim()}.
#'
#' @return A list containing the observed states,
#'   which are not necessaryily in a format ready for model-fitting
#'
#' @export

obs_sim = function(params, true) {

  output = with(append(params, true), {

    # generate lognormal observation error sd
    sig_S_ts_obs = sqrt(log(cv_S_ts_obs^2 + 1))
    sig_C_t_obs = sqrt(log(cv_C_t_obs^2 + 1))

    # observe spawners
    S_ts_obs = matrix(NA, nt, ns)
    for (s in 1:ns) {
      for (t in 1:nt) {
        S_ts_obs[t,s] = rlnorm(1, log(S_ts[t,s]), sig_S_ts_obs[t,s])
      }
    }

    # observe total harvest
    C_tot_t_obs = rlnorm(nt, log(C_tot_t), sig_C_t_obs)

    # observe age composition
    x_tas_obs = array(NA, dim = c(nt, na, ns))
    for (s in 1:ns) {
      for (t in 1:nt) {
        x_tas_obs[t,,s] = t(rmultinom(n = 1, size = x_ESS, prob = q_tas[t,,s]))
      }
    }

    N_tot_t_obs = C_tot_t_obs + rowSums(S_ts_obs)
    U_t_obs = C_tot_t_obs/N_tot_t_obs

    # bundle output
    list(
      C_tot_t_obs = C_tot_t_obs,
      S_ts_obs = S_ts_obs,
      x_tas_obs = x_tas_obs,
      sig_S_ts_obs = sig_S_ts_obs,
      sig_C_t_obs = sig_C_t_obs,
      U_t_obs = U_t_obs
      )
  })

  # return output
  return(output)
}




