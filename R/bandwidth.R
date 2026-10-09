#' Rule-of-thumb bandwidths based on the spread of entry ages
#'
#' Computes \eqn{h = c\,\hat\sigma_A\, n^{\alpha}}{h = c * sd(A) * n^alpha},
#' where \eqn{\hat\sigma_A}{sd(A)} is the sample standard deviation of the
#' entry ages and \eqn{n} is the number of subjects. This is the form of all
#' bandwidths used in the paper.
#'
#' @param entry_age Numeric vector of entry ages \eqn{A_i} (one per subject).
#'   Non-finite values are ignored.
#' @param exponent The exponent \eqn{\alpha}. The default \eqn{-1/6} gives the
#'   "small" bandwidth for the ratio estimator (\eqn{mh \to 0}); the paper's
#'   "large" bandwidth uses \eqn{-0.14} and the mean-trajectory bandwidth uses
#'   \eqn{-1/5}.
#' @param constant The multiplicative constant \eqn{c} (\eqn{\kappa} in the
#'   paper). Use `2.34` together with `exponent = -1/5` for the bandwidth
#'   \eqn{h'} of the mean-trajectory estimators (the rule-of-thumb constant of
#'   the Epanechnikov kernel), and a value selected by [cv_bandwidth()] for the
#'   ratio estimator.
#'
#' @return A single positive number, on the same scale as `entry_age`.
#'
#' @details [trajchain()] uses, by default, the ratio bandwidth
#'   \eqn{h = \hat\sigma_A n^{-1/6}}{h = sd(A) n^(-1/6)} and the
#'   mean-trajectory bandwidth \eqn{h' = 2.34\,\hat\sigma_A n^{-1/5}}{h' = 2.34 sd(A) n^(-1/5)},
#'   as in the authors' simulation and data-analysis code.
#'
#' @template ref
#' @seealso [cv_bandwidth()] to choose the constant by cross-validation.
#'
#' @examples
#' set.seed(1)
#' age <- runif(200, 25, 80)
#' bw_rot(age)                                   # ratio estimator
#' bw_rot(age, exponent = -1/5, constant = 2.34) # mean trajectory
#'
#' @export
bw_rot <- function(entry_age, exponent = -1 / 6, constant = 1) {
  if (!is.numeric(entry_age)) .stopf("`entry_age` must be numeric.")
  .check_number(exponent, "exponent")
  .check_number(constant, "constant", positive = TRUE)
  A <- as.numeric(entry_age)
  A <- A[is.finite(A)]
  n <- length(A)
  if (n < 2L) .stopf("At least two finite entry ages are needed to compute a bandwidth.")
  s <- stats::sd(A)
  if (s <= 0) .stopf("The entry ages have zero spread; a rule-of-thumb bandwidth is undefined.")
  constant * (s * n^exponent)
}

# Resolve the `bandwidth` argument of trajchain() for the ratio estimator.
# Returns list(h = <age units>, constant = <c or NA>, cv = <cv object or NULL>).
.resolve_ratio_bandwidth <- function(bandwidth, design, kern, cv_args = list()) {
  A <- design$entry_age
  if (inherits(bandwidth, "trajchain_cv")) {
    if (!identical(bandwidth$type, "bandwidth")) {
      .stopf("`bandwidth` must be a cross-validation object returned by cv_bandwidth().")
    }
    return(list(h = bandwidth$constant * bw_rot(A, -1 / 6), constant = bandwidth$constant,
                cv = bandwidth))
  }
  if (is.character(bandwidth)) {
    bandwidth <- match.arg(bandwidth, c("rot", "cv"))
    if (bandwidth == "rot") {
      return(list(h = bw_rot(A, -1 / 6), constant = 1, cv = NULL))
    }
    cv <- do.call(.cv_bandwidth_core, c(list(design = design, kern = kern), cv_args))
    return(list(h = cv$constant * bw_rot(A, -1 / 6), constant = cv$constant, cv = cv))
  }
  h <- .check_number(bandwidth, "bandwidth", positive = TRUE)
  list(h = h, constant = NA_real_, cv = NULL)
}

.resolve_mean_bandwidth <- function(bandwidth_mu, design) {
  if (is.character(bandwidth_mu)) {
    match.arg(bandwidth_mu, "rot")
    return(2.34 * bw_rot(design$entry_age, -1 / 5))
  }
  .check_number(bandwidth_mu, "bandwidth_mu", positive = TRUE)
}
