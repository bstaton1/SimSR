#' Obtain Brood Year Recruits
#'
#' Reconstruct the brood table to get the number of fish that were produced by a single
#' spawning event. For indexing, recruit abundance in brood year \code{a_max+1}
#' was produced by spawner abundance in calendar year 1.
#'
#'
#' @param params A list created using \code{\link{init_sim}}.
#' @param obs A list created using \code{\link{obs_sim}}.
#'
#' @details Currently, a brood year recruitment event becomes an NA if any
#'   one of the return years was missing in the escapement counts.
#'   For stocks without age composition sampling, the average
#'   age composition observed for monitored stocks is used to allocate
#'   abundance-at-age-and-year to brood year recruits. The assumption is made that
#'   all substocks received the same exploitation rate.
#'
#' @return An updated list which contains all of the contents of the \code{obs} list
#'   supplied as an argument
#'
#' @export

gen_Rys_obs = function(params, obs) {

  output = with(append(params, obs), {
    # containers
    H_ts_obs = matrix(NA, nt, ns)
    R_ys_obs = matrix(NA, ny, ns)
    S_tas_obs = array(NA, dim = c(nt, na, ns))
    H_tas_obs = array(NA, dim = c(nt, na, ns))
    N_tas_obs = array(NA, dim = c(nt, na, ns))

    # obtain age proportions for stocks that have age data
    q_tas_obs = array(NA, dim = c(nt, na, length(age_comp_stocks)))
    for (s in 1:length(age_comp_stocks)) {
      q_tas_obs[,,s] = t(apply(x_tas_obs[,,s], 1, function(x) x/sum(x)))
    }

    # obtain the average age proportion across all stocks that have data
    q_ta_ave = matrix(NA, nt, na)
    for (t in 1:nt) {
      for (a in 1:na) {
        q_ta_ave[t,a] = mean(q_tas_obs[t,a,], na.rm = T)
      }
    }

    for (s in 1:ns) {

      # determine appropriate age comps to use
      if (s %in% age_comp_stocks) {
        q_ta_use = q_tas_obs[,,which(age_comp_stocks == s)]
      } else {
        q_ta_use = q_ta_ave
      }

      # generate calendar year quantities for each stock
      H_ts_obs[,s] = (S_ts_obs[,s] * U_t_obs)/(1 - U_t_obs)
      S_tas_obs[,,s] = apply(q_ta_use, 2, function(x) x * S_ts_obs[,s])
      H_tas_obs[,,s] = apply(q_ta_use, 2, function(x) x * H_ts_obs[,s])
      N_tas_obs[,,s] = S_tas_obs[,,s] + H_tas_obs[,,s]

      # generate brood year recruits for each stock
      for (y in 1:ny) {
        if (y <= (nt - na + 1)) {
          brd.yr.runs = diag(N_tas_obs[y:(y+na-1),,s])
          R_ys_obs[y+na-1,s] = sum(brd.yr.runs, na.rm = all(!is.na(brd.yr.runs)))
        } else {
          next()
        }
      }
    }

    R_ys_obs[R_ys_obs == "NaN"] = NA

    # brood year indices
    S_ind = 1:(nt - a_max)
    R_ind = (a_max + 1):(ny - na + 1)

    list(
      R_ys_obs = R_ys_obs,
      n_S_obs = apply(S_ts_obs, 2, function(x) sum(!is.na(x))),
      n_R_obs = apply(R_ys_obs, 2, function(x) sum(!is.na(x))),
      n_SR_obs = sapply(1:ns, function(s) {
        sum(!is.na(S_ts_obs[S_ind,s]) & !is.na(R_ys_obs[R_ind,s]))
      })
    )
  })

  # add it to the observed data set
  obs = append(obs, output)

  # return output
  return(obs)
}
