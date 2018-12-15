#' Generate Quantities used in Management
#'
#' Based on vectors of substock-specific parameters, obtain the biological reference points
#'   for the aggregate stock if all substocks were fished at the same rate. Calculates
#'   a set of quantities called Sstar_p and Ustar_p. Sstar_p is the smallest aggregate escapement
#'   where you can still have no more than p*100% of the substocks not overfished.
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
    if (all(overfished > 0.1)) {
      Sstar_0.1 = NA
      Ustar_0.1 = NA
    } else {
      Sstar_0.1 = min(Seq[which(overfished <= 0.1)])
      Ustar_0.1 = max(U_range[which(overfished <= 0.1)])
    }
    if (all(overfished > 0.3)) {
      Sstar0.3 = NA
      Ustar_0.3 = NA
    } else {
      Sstar_0.3 = min(Seq[which(overfished <= 0.3)])
      Ustar_0.3 = max(U_range[which(overfished <= 0.3)])
    }
    if (all(overfished > 0.5)) {
      Sstar_0.5 = NA
      Ustar_0.5 = NA
    } else {
      Sstar_0.5 = min(Seq[which(overfished <= 0.5)])
      Ustar_0.5 = max(U_range[which(overfished <= 0.5)])
    }

    list(
      mgmt = c(
        Sstar_0.1 = Sstar_0.1, Sstar_0.3 = Sstar_0.3, Sstar_0.5 = Sstar_0.5,
        Ustar_0.1 = Ustar_0.1, Ustar_0.3 = Ustar_0.3, Ustar_0.5 = Ustar_0.5,
        S_MSY = S_MSY, U_MSY = U_MSY),
      Seq = Seq,
      Ceq = Ceq,
      Req_s = log(alpha)/beta,
      overfished = overfished,
      extinct = extinct
    )
  })

  return(output)
}
