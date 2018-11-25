#' Generate Quantities used in Management
#'
#' Based on vectors of substock-specific parameters, obtain the biological reference points
#'   for the aggregate stock if all substocks were fished at the same rate.
#'
#' @param params A list created using \code{init_sim()}.
#' @param U_range A sequence at which to calculate the different equilibrium states
#'
#' @export

gen_mgmt = function(params, U_range = seq(0, 1, 0.01)) {

  output = with(params, {

    # key parameters for each substock
    sub_params = cbind(alpha, beta, U_msy, S_msy)

    # determine the equilibrium quantities for each substock at various exploitation rates (U_range)
    Seq_s = apply(sub_params, 1, function(x) eq_ricker(alpha = x["alpha"], beta = x["beta"], U_msy = x["U_msy"], S_msy = x["S_msy"], U_range = U_range)$S)
    Ceq_s = apply(sub_params, 1, function(x) eq_ricker(alpha = x["alpha"], beta = x["beta"], U_msy = x["U_msy"], S_msy = x["S_msy"], U_range = U_range)$C)
    overfished_s = apply(sub_params, 1, function(x) eq_ricker(alpha = x["alpha"], beta = x["beta"], U_msy = x["U_msy"], S_msy = x["S_msy"], U_range = U_range)$overfished)
    extinct_s = apply(sub_params, 1, function(x) eq_ricker(alpha = x["alpha"], beta = x["beta"], U_msy = x["U_msy"], S_msy = x["S_msy"], U_range = U_range)$extinct)

    # sum across substocks
    Seq = rowSums(Seq_s)
    Ceq = rowSums(Ceq_s)
    overfished = rowSums(overfished_s)/ns
    extinct = rowSums(extinct_s)/ns

    # system-wide BRPs
    S_MSY = Seq[which.max(Ceq)]
    U_MSY = U_range[which.max(Ceq)]

    # system-wide MRPs
    if (all(overfished > max_p_overfished)) {
      S_obj = NA
      U_obj = NA
    } else {
      S_obj = min(Seq[which(overfished < max_p_overfished)])
      U_obj = max(U_range[which(overfished < max_p_overfished)])
    }

    # if(S_obj %in% c("Inf", "-Inf")) S_obj = NA
    # if(U_obj %in% c("Inf", "-Inf")) U_obj = NA
    #
    list(
      mgmt = c(S_obj = S_obj, U_obj = U_obj, S_MSY = S_MSY, U_MSY = U_MSY),
      Seq = Seq,
      Ceq = Ceq,
      Req_s = log(alpha)/beta,
      overfished = overfished,
      extinct = extinct
    )
  })

  return(output)
}
