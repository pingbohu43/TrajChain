#' Kernel estimator of the ratio of the shape function at two ages
#'
#' Computes the kernel ratio estimator of equation (4) of Hu et al.,
#' \deqn{\hat r(t, t - v) = \frac{\sum_{i=1}^n \int K_h(t-u)\, Y_i(u)\, dN_i(u)}
#'   {\sum_{i=1}^n \int K_h(t-u)\, Y_i(u - v)\, dN_i(u)},}{
#'   r(t, t - v) = sum_i int K_h(t - u) Y_i(u) dN_i(u) / sum_i int K_h(t - u) Y_i(u - v) dN_i(u),}
#' which estimates \eqn{r(t, t - v) = f(t)/f(t - v)}, the relative change of the
#' shape function between two ages \eqn{v} apart. The counting process
#' \eqn{N_i} jumps at the follow-up visits of subject \eqn{i} for which both the
#' current measurement and the measurement at the preceding visit (scheduled
#' \eqn{v} earlier) are observed.
#'
#' @param Y Numeric matrix of marker values with one row per subject and one
#'   column per scheduled visit; column 1 is the baseline measurement at study
#'   entry. Visits that were not observed (after death, dropout or end of
#'   follow-up) are `NA`.
#' @param entry_age Numeric vector of ages at study entry (\eqn{A_i}), one per
#'   row of `Y`.
#' @param visit_gaps Scheduled time between consecutive visits: a single number
#'   for equally spaced visits, or a vector of length `ncol(Y) - 1` giving
#'   \eqn{v_1, \dots, v_k}. Visit \eqn{l} of subject \eqn{i} takes place at age
#'   \eqn{A_i + v_1 + \dots + v_l}.
#' @param t Ages at which to evaluate the ratio.
#' @param gap The age difference \eqn{v} of the ratio. It must equal one of the
#'   gaps in `visit_gaps`; it can be omitted when all gaps are equal.
#' @param bandwidth Bandwidth \eqn{h} on the age scale, or `"rot"` for the
#'   rule-of-thumb \eqn{\hat\sigma_A n^{-1/6}}{sd(A) n^(-1/6)}.
#' @param kernel Kernel; see [kernel_function()].
#'
#' @return A data frame with one row per element of `t` and columns
#'   \describe{
#'     \item{`t`}{evaluation age;}
#'     \item{`ratio`}{\eqn{\hat r(t, t - v)} (`NA` when no observed pair falls
#'       in the kernel window);}
#'     \item{`numerator`, `denominator`}{the two kernel sums;}
#'     \item{`n_pairs`}{number of observed visit pairs whose current visit lies
#'       within one bandwidth of `t`.}
#'   }
#'
#' @template ref
#' @seealso [trajchain()], which chains these ratios into the shape function.
#'
#' @examples
#' sim <- simulate_cohort(n = 500, seed = 1)
#' ratio_estimate(sim$Y, sim$entry_age, sim$visit_gaps, t = c(0.3, 0.5, 0.7))
#'
#' # true ratios f(t) / f(t - v)
#' v <- sim$visit_gaps[1]
#' shape_monotone(c(0.3, 0.5, 0.7)) / shape_monotone(c(0.3, 0.5, 0.7) - v)
#'
#' @export
ratio_estimate <- function(Y, entry_age, visit_gaps, t, gap = NULL,
                           bandwidth = "rot", kernel = "epanechnikov") {
  Y <- .check_Y(Y)
  A <- .check_entry_age(entry_age, nrow(Y))
  gaps <- .check_gaps(visit_gaps, ncol(Y) - 1L)
  if (!is.numeric(t) || !length(t) || any(!is.finite(t))) {
    .stopf("`t` must be a numeric vector of finite ages.")
  }
  kern <- kernel_function(kernel)
  if (is.null(gap)) {
    if (any(abs(gaps - gaps[1L]) > 1e-8 * gaps[1L])) {
      .stopf("The visits are unequally spaced; specify which `gap` the ratio refers to.")
    }
    gap <- gaps[1L]
  }
  .check_number(gap, "gap", positive = TRUE)
  cols <- which(abs(gaps - gap) <= 1e-8 * max(gap, gaps)) + 1L
  if (!length(cols)) {
    .stopf("`gap` = %s does not match any of the visit gaps (%s).", .format_num(gap),
           paste(.format_num(unique(gaps)), collapse = ", "))
  }
  if (is.character(bandwidth)) {
    match.arg(bandwidth, "rot")
    h <- bw_rot(A, -1 / 6)
  } else {
    h <- .check_number(bandwidth, "bandwidth", positive = TRUE)
  }
  equal <- all(abs(gaps - gaps[1L]) <= 1e-12 * gaps[1L])
  offsets <- if (equal) gaps[1L] * (0:length(gaps)) else c(0, cumsum(gaps))
  obs <- !is.na(Y)
  u <- ycur <- ylag <- numeric(0)
  for (cc in cols) {
    ok <- obs[, cc] & obs[, cc - 1L]
    u <- c(u, A[ok] + offsets[cc])
    ycur <- c(ycur, Y[ok, cc])
    ylag <- c(ylag, Y[ok, cc - 1L])
  }
  S <- .ksum(t, u, cbind(ycur, ylag), h, kern)
  ratio <- S[, 1L] / S[, 2L]
  ratio[S[, 2L] == 0] <- NA_real_
  data.frame(t = t, ratio = ratio, numerator = S[, 1L], denominator = S[, 2L],
             n_pairs = .count_window(t, u, h))
}
