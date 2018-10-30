#' Create a plot of the true states
#'
#' @param params A list created with \code{init_sim()}
#' @param true A list created with \code{ricker_sim()}
#' @param new_window logical. Do you want the plots to open in a new window?
#' @param type character vector. Current accepted types are \code{"SRA"} and \code{"spawners"}
#'
#' @export

true_plot = function(params, true, new_window = T, type = "SRA") {

  if (new_window) windows()

  with(append(params, true), {

    if (type == "SRA") {
      S_ind = 1:(nt - a_min)
      R_ind = (a_max + 1):ny
      par(mfrow = c(3,5), mar = c(2,2,2,2))
      for (s in 1:ns) {
        # plot(R_ts[,s], type = "l")
        plot(R_ys[R_ind,s] ~ S_ts[S_ind,s], pch = 16, ylim = c(0, max(R_ys[R_ind,s])), xlim = c(0, max(S_ts[S_ind,s])))
        curve(x * exp(log_alpha[s] - beta[s] * x), from = 0, to = max(S_ts[S_ind,s]), add = T)
        lines(R_ys[R_ind,s] ~ S_ts[S_ind,s], col = "grey")
        abline(c(0,1), col = "grey", lty = 2)
      }
    } else {
      if (type == "spawners") {
        par(mfrow = c(3,5), mar = c(2,2,2,2))
        for (s in 1:ns) {
          plot(S_ts[,s], type = "o", pch = 16)
        }
      } else {
        matplot(log_resid_ys, type = "l", col = "grey", lty = 1, main = "log_resids")
        lines(w_y, lwd = 2)
        abline(h = 0, lty = 2)
      }
    }
  })
}

