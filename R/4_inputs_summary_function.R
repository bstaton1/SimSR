#' Summarize important input quantities
#'
#' Reformats key information from the
#'   three primary input lists for writing to output files
#'
#' @param params a list object created with \code{init_sim()}
#' @param obs a list object created with the observation generation
#'   functions workflow.
#' @param true a list object created by \code{ricker_sim()}
#' @param seed An identifier for this set of inputs
#'
#' @export

input_summary = function(params, obs, true, seed) {

  output = with(append(append(params, obs), true), {
    # get mean rho
    rho_mat2 = rho_mat
    diag(rho_mat2) = NA
    mean_rho = mean(rho_mat2, na.rm = T)
    mean_sigma_R = mean(sigma)

    mrp = gen_mgmt(params)$mgmt

    # names of the things being summarized
    leading_p = c(rep("alpha", ns), rep("beta", ns), rep("sigma_R", ns), rep("U_msy", ns), rep("S_msy", ns), rep("pi", na))
    derived_p = c("mean_rho", "mean_sigma_R")
    mrp_p = names(mrp)
    n_obs_p = c(rep("n_S_obs", ns), rep("n_R_obs", ns), rep("n_SR_obs", ns), rep("age_stock", ns))
    state_p = c(rep("R", ns * ny), rep("S", ns * nt), rep("H", nt), rep("U", nt))
    obs_state_p = c(rep("R_obs", ns * ny), rep("S_obs", ns * nt), rep("H_obs", nt), rep("U_obs", nt))

    # the stock attribute for each thing being summarized
    leading_s = c(rep(1:ns, 5), rep(NA, na))
    derived_s = c(rep(NA, 2))
    mrp_s = rep(NA, length(mrp_p))
    n_obs_s = rep(1:ns, 4)
    state_s = c(rep(1:ns, each = ny), rep(1:ns, each = nt), rep(NA, nt * 2))
    obs_state_s = c(rep(1:ns, each = ny), rep(1:ns, each = nt), rep(NA, nt * 2))

    # the year attribute for each thing being summarized
    leading_y = rep(NA, length(leading_s))
    derived_y = rep(NA, length(derived_s))
    mrp_y = rep(NA, length(mrp_p))
    n_obs_y = rep(NA, length(n_obs_s))
    state_y = c(rep(1:ny, ns), rep(1:nt, ns), rep(1:nt, 2))
    obs_state_y = c(rep(1:ny, ns), rep(1:nt, ns), rep(1:nt, 2))

    # value of the parameter/quantity
    leading_v = c(alpha, beta, sigma, U_msy, S_msy, pi_grand)
    derived_v = c(mean_rho, mean_sigma_R)
    mrp_v = unname(mrp)
    n_obs_v = c(n_S_obs, n_R_obs, n_SR_obs, ifelse(1:ns %in% age_comp_stocks, 1, 0))
    state_v = c(as.numeric(R_ys), as.numeric(S_ts), C_tot_t, U_real)
    obs_state_v = c(as.numeric(R_ys_obs), as.numeric(S_ts_obs), C_tot_t_obs, U_t_obs)

    # combine to data.frames
    leading_df = data.frame(seed = seed, stock = leading_s, year = leading_y, param = leading_p, value = leading_v)
    derived_df = data.frame(seed = seed, stock = derived_s, year = derived_y, param = derived_p, value = derived_v)
    mrp_df = data.frame(seed = seed, stock = mrp_s, year = mrp_y, param = mrp_p, value = mrp_v)
    n_obs_df = data.frame(seed = seed, stock = n_obs_s,year = n_obs_y, param = n_obs_p, value = n_obs_v)
    state_df = data.frame(seed = seed, stock = state_s, year = state_y, param = state_p, value = state_v)
    obs_state_df = data.frame(seed = seed, stock = obs_state_s, year = obs_state_y, param = obs_state_p, value = obs_state_v)

   rbind(
     leading_df,
     derived_df,
     mrp_df,
     n_obs_df,
     state_df,
     obs_state_df
   )
  })

  # return output
  return(output)
}
