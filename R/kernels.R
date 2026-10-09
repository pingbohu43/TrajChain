#' Kernel functions with compact support
#'
#' Returns one of the kernel functions supported by the package. All kernels
#' are symmetric probability densities supported on \eqn{[-1, 1]}, as required
#' by the theory in Hu et al. The default, the Epanechnikov kernel
#' \eqn{K(x) = (3/4)(1 - x^2) I(|x| \le 1)}{K(x) = 0.75 (1 - x^2) I(|x| <= 1)},
#' is the kernel used throughout the paper.
#'
#' @param kernel Either the name of a kernel (`"epanechnikov"`, `"biweight"`,
#'   `"triweight"`, `"triangular"` or `"uniform"`), or a user-supplied
#'   vectorized function. A user-supplied function must return a non-negative
#'   value for every element of its argument, preserve the dimensions of a
#'   matrix argument and vanish outside \eqn{[-1, 1]}.
#'
#' @return A function of one argument `u` (a numeric vector or matrix) that
#'   returns \eqn{K(u)} with the same dimensions as `u`.
#'
#' @details Kernel weights enter the estimators only through ratios, so the
#'   normalizing factor \eqn{1/h} of \eqn{K_h(\cdot) = K(\cdot/h)/h} cancels
#'   and is omitted. The rule-of-thumb constants used by [bw_rot()] (for
#'   example \eqn{2.34} for the mean-trajectory bandwidth) were calibrated for
#'   the Epanechnikov kernel.
#'
#' @examples
#' K <- kernel_function("epanechnikov")
#' K(c(-1.5, -0.5, 0, 0.5, 1.5))
#' integrate(K, -1, 1)$value
#'
#' @export
kernel_function <- function(kernel = "epanechnikov") {
  if (is.function(kernel)) return(.validate_kernel(kernel))
  if (!is.character(kernel) || length(kernel) != 1L) {
    .stopf("`kernel` must be a kernel name or a function.")
  }
  kernel <- match.arg(kernel, c("epanechnikov", "biweight", "triweight",
                                "triangular", "uniform"))
  switch(kernel,
    epanechnikov = function(u) {
      out <- 0.75 * (1 - u^2)
      out[abs(u) > 1] <- 0
      out
    },
    biweight = function(u) {
      out <- (15 / 16) * (1 - u^2)^2
      out[abs(u) > 1] <- 0
      out
    },
    triweight = function(u) {
      out <- (35 / 32) * (1 - u^2)^3
      out[abs(u) > 1] <- 0
      out
    },
    triangular = function(u) {
      out <- 1 - abs(u)
      out[abs(u) > 1] <- 0
      out
    },
    uniform = function(u) {
      out <- 0.5 + 0 * u
      out[abs(u) > 1] <- 0
      out
    }
  )
}

.validate_kernel <- function(kern) {
  probe <- matrix(c(-2, -1.5, -0.5, 0, 0.5, 1.5, 2, 0.25), nrow = 2)
  val <- tryCatch(kern(probe), error = function(e) NULL)
  if (is.null(val) || !is.numeric(val) || length(val) != length(probe)) {
    .stopf("A user-supplied `kernel` must be a vectorized function returning one value per input.")
  }
  if (!identical(dim(val), dim(probe))) {
    .stopf("A user-supplied `kernel` must preserve the dimensions of a matrix argument.")
  }
  if (any(!is.finite(val)) || any(val < 0)) {
    .stopf("A user-supplied `kernel` must return finite, non-negative values.")
  }
  if (any(val[abs(probe) > 1] != 0)) {
    .stopf("A user-supplied `kernel` must vanish outside [-1, 1] (compact support is required).")
  }
  kern
}

.kernel_label <- function(kernel) {
  if (is.function(kernel)) "user-supplied" else match.arg(kernel, c(
    "epanechnikov", "biweight", "triweight", "triangular", "uniform"))
}
