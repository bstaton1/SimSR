#' Summarize true parameters
#'
#' Reformats key output from the \code{params} and \code{obs}
#' lists for writing to output files.
#'
#' @param params a list object created with \code{init_sim()}
#' @param obs a list object created with the observation generation
#'   functions workflow.
#' @param seed An identifier for this set of parameters
#'
#' @export

params_summary = function(params, obs, seed) {

  output = with(append(params, obs), {
    # get mean rho
    rho_mat2 = rho_mat
    diag(rho_mat2) = NA
    mean_rho = mean(rho_mat2, na.rm = T)

    # names of the parameters
    p = c(rep("alpha", ns), rep("beta", ns), rep("sigma_R", ns), rep("U_msy", ns),
          rep("S_msy", ns), "mean_rho", "mean_sigma_R", "S_obj", "U_obj", "U_MSY", "S_MSY",
          rep("n_S_obs", ns), rep("n_R_obs", ns), rep("n_SR_obs", ns), rep("age_stock", ns))

    # stock identifier
    s = c(rep(1:ns, 5), rep(NA, 6), rep(1:ns, 4))

    # value of the parameter/quantity
    v = c(alpha, beta, sigma, U_msy, S_msy, mean_rho, mean(sigma))

    # get drainage-wide parameter estimates
    dw = unname(gen_mgmt(params = params)$mgmt[c("S_obj", "U_obj", "U_MSY", "S_MSY")])
    v = c(v, dw)

    # include how many years were observed of different types
    v = c(v, n_S_obs, n_R_obs, n_SR_obs, ifelse(1:ns %in% age_comp_stocks, 1, 0))

    data.frame(seed = seed, stock = s, param = p, value = v)

  })

  # return output
  return(output)
}
