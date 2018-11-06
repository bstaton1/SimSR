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
      par(mfrow = c(3,5), mar = c(2,1.5,1,1), oma = c(0, 0.5, 0, 0))
      for (s in 1:ns) {
        # plot(R_ts[,s], type = "l")
        plot(R_ys[R_ind,s] ~ S_ts[S_ind,s], pch = 16, ylim = c(0, max(R_ys[R_ind,s])), xlim = c(0, max(log(alpha[s])/beta[s], S_ts[S_ind,s])))
        curve(x * exp(log_alpha[s] - beta[s] * x), from = 0, to = max(log(alpha[s])/beta[s], S_ts[S_ind,s]), add = T)
        usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
        text(x = usr[2], y = usr[4] - 0.1 * ydiff, paste("U_msy =", round(U_msy[s], 2)), pos = 2, col = "grey50", font = 2)
        text(x = usr[2], y = usr[4] - 0.2 * ydiff, paste("S_msy =", round(S_msy[s], -2)), pos = 2, col = "grey50", font = 2)

        # lines(R_ys[R_ind,s] ~ S_ts[S_ind,s], col = "grey")
        abline(c(0,1), col = "grey", lty = 2)
      }
      plot(1, 1, type = "n", ann = F, axes = F,
           xlim = c(0,1), ylim = c(0,1))
      text(x = 0.05, y = 0.8, pos = 4, paste("Mean U =", round(mean(U_real),2)))
      text(x = 0.05, y = 0.7, pos = 4, paste("Min U =", round(min(U_real), 2)))
      text(x = 0.05, y = 0.6, pos = 4, paste("Max U =", round(max(U_real), 2)))
      box()
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

