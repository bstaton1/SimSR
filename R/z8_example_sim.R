#' An example simulation wrapper
#'
#' @param U_SUM passed to \code{init_sim()}
#' @param min_sigR passed to \code{init_sim()}
#' @param max_sigR passed to \code{init_sim()}
#' @param max_p_overfished passed to \code{init_sim()}
#' @param do_plot logical. Do you want to make an SRA plot?
#' @param save logical. Do you want to save a file containing the plot?
#'   Only takes effect if \code{do_plot = TRUE}
#'
#' @examples
#' example_sim(
#'   U_SUM = 100,
#'   min_sigR = 0.3,
#'   max_sigR = 0.5,
#'   max_p_overfished = 0.3,
#'   do_plot = F,
#'   save = F
#'   )
#'
#' @export

example_sim = function(U_SUM = 100, min_sigR = 0.3, max_sigR = 0.5,
                       max_p_overfished = 0.3, do_plot = F, save = F) {
  Umsy = c(
    0.6400932, 0.4595393, 0.5518526, 0.3730096,
    0.6774239, 0.4963993, 0.4309589, 0.4991643,
    0.7116140, 0.6075838, 0.4382813, 0.5283336
  )

  Smsy = c(
    7974,  3252, 10832,  1064,
    2576,  6330,  7331,   623,
    2495,  1702,   274,  1106
  )

  # settings for operating model
  params = init_sim(
    U_msy = Umsy,
    S_msy = Smsy,
    U_SUM = U_SUM,
    min_sigR = min_sigR,
    max_sigR = max_sigR,
    max_p_overfished = max_p_overfished
  )

  # create true states
  true = ricker_sim(params)

  # observed states
  obs = obs_sim(params, true)
  obs = obs_filter(params, obs)

  ppi = 600
  if (save & do_plot) png("SRA.png", h = 5 * ppi, w = 8 * ppi, res = ppi)
  if (do_plot) {
    true_plot(params = params,
              true = true, new_window = F,
              type = "SRA")
    matplot(true$log_resid_ys, type = "l", col = "grey")
    lines(true$w_y, lwd = 2)
  }
  if (save & do_plot) dev.off()

  list(
    params = params,
    true = true,
    obs = gen_Rys_obs(params, obs)
  )
}

example_sim(
  U_SUM = 100,
  min_sigR = 0.3,
  max_sigR = 0.5,
  max_p_overfished = 0.3,
  do_plot = T,
  save = T
)
