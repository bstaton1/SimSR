#' Generate True States
#'
#' Simulates Ricker dynamics based on the settings supplied via the input argument.
#'
#' @param params A list created using \code{init_sim()}.
#'
#' @return A list with the true states simulated.
#'
#' @export

ricker_sim = function(params) {

  output = with(params, {

    # obtain management reference points
    mgmt = gen_mgmt(params = params)
    Req_s = mgmt$Req_s
    U_obj = unname(mgmt$mgmt["U_obj"])

    # year zero resid
    log_resid_0 = 0

    # containers (brood years)
    log_R_ys_mean1 = matrix(NA, ny, ns)
    log_R_ys_mean2 = matrix(NA, ny, ns)
    R_ys = matrix(NA, ny, ns)
    log_R_ys = matrix(NA, ny, ns)
    log_resid_ys = matrix(NA, ny, ns)
    p_yas = array(NA, dim = c(ny, na, ns))

    # containers (calendar years)
    S_ts = matrix(NA, nt, ns)
    N_tas = array(NA, dim = c(nt, na, ns))
    N_ts = matrix(NA, nt, ns)
    C_ts = matrix(NA, nt, ns)
    N_tot_t = rep(NA, nt)
    U_real = rep(NA, nt)

    # create time-varying and substock-specific maturity schedules
    B_grand = c(-1.4, 1.4, 4)
    m = rbind(
      c(1,0,0),
      c(1,1,0),
      c(1,0,1)
    )
    stock_effect = rnorm(ns, 0, 0.2)
    Sigma_mat = matrix(0.1 * 0.1 * 0.90, ns, ns)
    diag(Sigma_mat) = rep(0.1^2, ns)
    stock_year_effect = mvtnorm::rmvnorm(ny, mean = stock_effect, Sigma_mat)
    for (y in 1:ny) {
      for (s in 1:ns) {
        p_yas[y,,s] = prob2pi(expit(m %*% (B_grand + c(stock_year_effect[y,s],0,0))))
      }
    }
    pi = prob2pi(expit(m %*% B_grand))

    # first brood year recruits
    log_R_ys_mean1[1,] = log(Req_s)
    log_R_ys_mean2[1,] = log_R_ys_mean1[1,] + phi * log_resid_0
    log_R_ys[1,] = mvtnorm::rmvnorm(1, log_R_ys_mean2[1,], Sigma)
    R_ys[1,] = exp(log_R_ys[1,])
    log_resid_ys[1,] = log_R_ys[1,] - log_R_ys_mean1[1,]

    # 2:a_max brood year recruits
    for (y in 2:a_max) {
      for (s in 1:ns) {
        log_R_ys_mean1[y,s] = log(Req_s[s])
        log_R_ys_mean2[y,s] = log_R_ys_mean1[y,s] + phi * log_resid_ys[y-1,s]
      }
      log_R_ys[y,] = mvtnorm::rmvnorm(1, log_R_ys_mean2[y,], Sigma)
      R_ys[y,] = exp(log_R_ys[y,])
      log_resid_ys[y,] = log_R_ys[y,] - log_R_ys_mean1[y,]
    }

    # fill N_tas with these recruits
    for (s in 1:ns) {
      for (t in 1:a_max) {
        for (a in 1:na) {
          # N_tas[t,a,s] = R_ys[t+na-a,s] * pi[a]
          N_tas[t,a,s] = R_ys[t+na-a,s] * p_yas[t+na-a,a,s]
        }
      }
    }

    # initialize N_ts, S_ts, C_ts
    for (t in 1:a_max) {
      for (s in 1:ns) {
        N_ts[t,s] = sum(N_tas[t,1:na,s])
      }

      N_tot_t[t] = sum(N_ts[t,])
      U_real[t] = implement_error(U_target = U_obj, SUM = U_SUM)
      for (s in 1:ns) {
        C_ts[t,s] = N_ts[t,s] * (U_real[t] * v[s])
        S_ts[t,s] = N_ts[t,s] * (1 - U_real[t] * v[s])
      }
    }

    # carry out the rest of the time series
    for (y in (a_max+1):ny) {
      # create brood year recruits
      for (s in 1:ns) {
        log_R_ys_mean1[y,s] = log(S_ts[y-a_max,s] * exp(log_alpha[s] - beta[s] * S_ts[y-a_max,s]))
        log_R_ys_mean2[y,s] = log_R_ys_mean1[y,s] + phi * log_resid_ys[y-1,s]
      }
      log_R_ys[y,] = mvtnorm::rmvnorm(1, log_R_ys_mean2[y,], Sigma)
      R_ys[y,] = exp(log_R_ys[y,])
      log_resid_ys[y,] = log_R_ys[y,] - log_R_ys_mean1[y,]

      for (s in 1:ns) {
        # place these recruits in calendar year at age
        for (a in 1:na) {
          if (y-a_min+a > nt) {
            next()
          } else {
            N_tas[y-a_min+a,a,s] = R_ys[y,s] * p_yas[t+na-a,a,s]
          }
        }

        # calendar year processes at the first full calendar year
        t = y-a_min+1
        N_ts[t,s] = sum(N_tas[t,1:na,s])
      }

      N_tot_t[t] = sum(N_ts[t,])
      U_real[t] = implement_error(U_target = U_obj, SUM = U_SUM)
      for (s in 1:ns) {
        S_ts[t,s] = N_ts[t,s] * (1 - U_real[t] * v[s])
        C_ts[t,s] = N_ts[t,s] * (U_real[t] * v[s])
      }
    }

    # totals across stocks
    C_tot_t = rowSums(C_ts)
    N_tot_ta = apply(N_tas, 2, rowSums)

    q_ta = apply(N_tot_ta, 2, function(x) x/N_tot_t)

    w_y = rowMeans(log_resid_ys)

    q_tas = array(NA, dim = dim(N_tas))
    for (s in 1:ns) {
      q_tas[,,s] = t(apply(N_tas[,,s], 1, function(x) x/sum(x)))
    }

    # package output
    list(
      R_ys = R_ys,
      log_resid_ys = log_resid_ys,
      U_real = U_real,
      N_ts = N_ts,
      S_ts = S_ts,
      C_tot_t = C_tot_t,
      w_y = w_y,
      q_ta = q_ta,
      q_tas = q_tas,
      p_yas = p_yas
      )
  })

  # return output
  return(output)
}

