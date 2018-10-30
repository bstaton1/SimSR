#' Sample Observed States
#'
#' Based on the observed states obtained using \code{obs_sim},
#' introduce a sampling frequency.
#'
#' @param params A list created using \code{init_sim()}.
#' @param obs A list created using \code{obs_sim()}.
#'
#' @return A list containing the filtered observed states.
#'   Any year not sampled will be an NA.
#' @examples
#' x = replicate(100, sum(sample_TS(params, minSRobs = 5, size = 7)))
#' hist(x)
#' count_obs_sr_pairs(sample_TS(params, minSRobs = 5))
#'
#' @export

obs_filter = function(params, obs, mimic = NULL) {

  output = with(append(params, obs), {

    ### ESCAPEMENT OBSERVATION FILTERING ###
    S_ts_obs_filtered = S_ts_obs
    # S_ts_obs_filtered = matrix(0, nt, ns)

    # matrix of observed years for each stock
    NA_yrs = replicate(ns, !as.logical(sample_TS(params, mimic = mimic, minSRobs = 5)))
    # S_ts_obs_filtered[NA_yrs] = NA

    for (s in 1:ns) {
      S_ts_obs_filtered[NA_yrs[,s],s] = NA
    }

    ### AGE COMP DATA FILTERING ###
    age_comp_stocks = sample(1:ns, ceiling(ns/2))
    x_tas_obs_filtered = array(NA, dim = c(nt, na, length(age_comp_stocks)))
    for (s in 1:length(age_comp_stocks)) {
      x_tas_obs_filtered[!NA_yrs[,age_comp_stocks[s]],,s] =
        x_tas_obs[!NA_yrs[,age_comp_stocks[s]],,age_comp_stocks[s]]
    }
    # x_tas_obs_filtered[is.na(x_tas_obs_filtered)] = 0

    # return output
    list(
      S_ts_obs_filtered = S_ts_obs_filtered,
      x_tas_obs_filtered = x_tas_obs_filtered,
      age_comp_stocks = age_comp_stocks
    )

  })

  # change out the full observed time series with the filtered time series
  obs$S_ts_obs = output$S_ts_obs_filtered
  obs$x_tas_obs = output$x_tas_obs_filtered
  obs$age_comp_stocks = output$age_comp_stocks

  return(obs)
}

count_obs_sr_pairs = function(x) {
  N_ta = matrix(1, length(x), ncol = 4)
  N_ta[x == 0,] = NA

  nt = length(x)
  na = 4
  a_max = 7
  ny = nt + na - 1
  R_y = rep(NA, ny)
  for (y in 1:ny) {
    if (y <= (nt - 4)) {
      brd.yr.runs = diag(N_ta[y:(y+na),])
      R_y[y+na-1] = sum(brd.yr.runs, na.rm = all(!is.na(brd.yr.runs)))
    } else {
      next()
    }
  }

  S_ind = 1:(nt - a_max - 1)
  R_ind = (a_max + 1):(ny - na)

  sum(!is.na(x[S_ind]) & !is.na(R_y[R_ind]))
}

sample_TS = function(params, mimic = NULL, minSRobs = 3, size = 7) {

  with(params, {
    if (is.null(mimic)) {
      mimic = data.frame(
        stock = rep(1:ns, each = nt),
        obs = rbinom(nt * ns, size = 1, prob = 0.6)
      )
    } else {
      if (any(colnames(mimic) %!in% c("stock", "obs"))) {
        stop("mimic must have columns 'stock' and 'obs'")
      }
    }

    mimic$stratum = rep(rep(1:ceiling(nt/size), each = size), ns)[1:nrow(mimic)]
    counts = with(mimic, tapply(obs, list(stratum, stock), sum))

    counts = reshape2::melt(counts)
    colnames(counts) = c("stratum", "stock", "count")
    counts = cbind(counts, nocount = size - counts$count)


    fit = glm(cbind(count, nocount) ~ stratum,
              data = counts, family = binomial)

    # funciton to generate the years sampled in a strata for a stock
    sample_strata = function(fit, stratum) {
      # obtain probability any year in this strata was sampled
      p = predict(fit,
                  newdata = data.frame(stratum = stratum),
                  type = "response")

      # generate the vector
      rbinom(n = size, size = 1, prob = p)
    }

    # create the time series sampled
    # if the number of SR pair years is less than a threshold, try again
    x = rep(0, nt)
    while(count_obs_sr_pairs(x) < minSRobs) {
      # cat("\r", i)
      x = as.numeric(sapply(1:(nt/size), function(s) sample_strata(fit, s)))
      # i = i + 1
    }

    x
  })
}
