#' Sample Observed States
#'
#' Based on the observed states obtained using \code{obs_sim()},
#' introduce a sampling frequency (i.e., insert missing years of data)
#'
#' @param params A list created using \code{init_sim()}.
#' @param obs A list created using \code{obs_sim()}.
#' @param mimic A data frame with columns \code{"stock"} (a stock identifier)
#'   and \code{"obs"} (a binary indicator for whether that stock was monitored that year).
#'   This function attempts to mimic this pattern of observation frequency. If \code{NULL} (default),
#'   then each year and stock will have a 60\% of being monitored.
#' @param minSRobs Numeric vector of length 1: the minimum number of fully observed
#'   spawner and recruit brood year pairs allowed for any given stock. Defaults to 3.
#' @param p_age Numeric vector of length 1: the fraction of the stocks to be randomly-assigned
#'   to have age composition monitored in the same years they have escapement monitored.
#'   Defaults to 0.5, and the number of stocks is rounded up if a fractional number is specified.
#'
#' @return A list containing the filtered observed states.
#'   Any year not sampled will be an NA.
#' @export

obs_filter = function(params, obs, mimic = NULL, minSRobs = 3, p_age = 0.5) {

  output = with(append(params, obs), {

    ### ESCAPEMENT OBSERVATION FILTERING ###
    S_ts_obs_filtered = S_ts_obs

    # matrix of observed years for each stock
    NA_yrs = replicate(ns, !as.logical(sample_TS(params, mimic = mimic, minSRobs = minSRobs)))

    for (s in 1:ns) {
      S_ts_obs_filtered[NA_yrs[,s],s] = NA
    }

    ### AGE COMP DATA FILTERING ###
    age_comp_stocks = sort(sample(1:ns, ceiling(ns * p_age)))
    x_tas_obs_filtered = array(NA, dim = c(nt, na, length(age_comp_stocks)))
    for (s in 1:length(age_comp_stocks)) {
      x_tas_obs_filtered[!NA_yrs[,age_comp_stocks[s]],,s] =
        x_tas_obs[!NA_yrs[,age_comp_stocks[s]],,age_comp_stocks[s]]
    }

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

sample_TS = function(params, mimic = NULL, minSRobs = 3, size = 7) {
  require(StatonMisc)
  with(params, {

    # determine the frequency to mimic
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


    fit = suppressWarnings(
      glm(cbind(count, nocount) ~ stratum,
          data = counts, family = binomial)
    )

    # function to generate the years sampled in a strata for a stock
    sample_strata = function(fit, stratum) {
      # obtain probability any year in this strata was sampled
      p = predict(
        fit,
        newdata = data.frame(stratum = stratum),
        type = "response")

      # generate the vector
      rbinom(n = size, size = 1, prob = p)
    }

    # create the time series sampled
    # if the number of SR pair years is less than a threshold, try again
    x = rep(0, nt)
    while(count_obs_sr_pairs(x, a_max, na) < minSRobs) {
      # cat("\r", i)
      x = as.numeric(sapply(1:(nt/size), function(s) sample_strata(fit, s)))
      # i = i + 1
    }

    x
  })
}

