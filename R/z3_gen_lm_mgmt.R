#' Obtain estimates of U_msy and S_msy given alpha and beta
#'
#' Converts alpha and beta into Smsy and Umsy.
#'
#' @param @alpha a numeric vector representing the alpha parameter from a SRA.
#'   Can be length > 1.
#' @param @beta a numeric vector representing the beta parameter from a SRA.
#'   Can be length > 1.
#'
#' @export

gen_lm_mgmt = function(alpha, beta) {
  log_alpha = log(alpha)

  U_msy = log_alpha * (0.5 - (0.65 * log_alpha ^1.27)/(8.7 + log_alpha^1.27))
  U_msy[U_msy == "NaN"] = 0
  U_msy[U_msy < 0] = 0
  S_msy = U_msy/beta

  return(list(U_msy = U_msy, S_msy = S_msy))
}
