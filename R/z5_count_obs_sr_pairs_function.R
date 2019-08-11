#' Count the number of observed SR pairs
#'
#' @param x Numeric vector consisting of 0 (not-observed) and 1 (observed) values
#'   with length equal to the number of possibly observable years
#' @param na Numeric vector of length 1: the ages at which fish can return
#' @param a_max Numeric vector of length 1: the maximum age at which fish can return
#' @export

count_obs_sr_pairs = function(x, a_max, na) {
  nt = length(x)
  N_ta = matrix(1, nt, na)
  N_ta[x == 0,] = 0
  ny = nt + na - 1
  R_y = rep(NA, ny)
  for (y in 1:ny) {
    if (y <= (nt - na + 1)) {
      brd.yr.runs = diag(N_ta[y:(y+na-1),])
      R_y[y+na-1] = sum(brd.yr.runs)
    } else {
      next()
    }
  }

  S_ind = 1:(nt - a_max)
  R_ind = (a_max + 1):(ny - na + 1)

  sum((x[S_ind] == 1) & (R_y[R_ind] == na))
}
