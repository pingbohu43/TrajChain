# Internal utilities: error helpers, argument checks, a local RNG seed,
# kernel sums and a numerical greatest common divisor.

.stopf <- function(...) stop(sprintf(...), call. = FALSE)
.warnf <- function(...) warning(sprintf(...), call. = FALSE)

.is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.check_number <- function(x, name, positive = FALSE, nonnegative = FALSE) {
  if (!.is_number(x)) .stopf("`%s` must be a single finite number.", name)
  if (positive && x <= 0) .stopf("`%s` must be positive.", name)
  if (nonnegative && x < 0) .stopf("`%s` must be non-negative.", name)
  invisible(as.numeric(x))
}

.check_count <- function(x, name, min = 1L) {
  if (!.is_number(x) || abs(x - round(x)) > 1e-8 || x < min) {
    .stopf("`%s` must be a whole number >= %d.", name, as.integer(min))
  }
  invisible(as.integer(round(x)))
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .stopf("`%s` must be TRUE or FALSE.", name)
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Local RNG seed. When `seed` is not NULL the seed is set and the caller's
# random number stream is restored when the calling function exits, so that
# functions of this package never change the user's global RNG state.
# ---------------------------------------------------------------------------
.get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

.check_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!.is_number(seed) || abs(seed - round(seed)) > 0) {
    .stopf("`seed` must be NULL or a single whole number.")
  }
  as.integer(seed)
}

# ---------------------------------------------------------------------------
# Kernel-weighted sums.
#   .ksum(s, u, W, h, kern)[j, ] = sum_i K((s_j - u_i) / h) * W[i, ]
# Computed in blocks of evaluation points so that memory stays bounded.
# ---------------------------------------------------------------------------
.ksum <- function(s, u, W, h, kern, block = 2e6) {
  W <- as.matrix(W)
  ns <- length(s)
  out <- matrix(0, nrow = ns, ncol = ncol(W))
  nu <- length(u)
  if (ns == 0L || nu == 0L) return(out)
  step <- max(1L, as.integer(block %/% nu))
  for (first in seq.int(1L, ns, by = step)) {
    idx <- first:min(ns, first + step - 1L)
    K <- kern(outer(s[idx], u, "-") / h)
    out[idx, ] <- K %*% W
  }
  out
}

# Number of points u_i with s_j - h <= u_i <= s_j + h, for each s_j.
.count_window <- function(s, u, h) {
  if (length(u) == 0L) return(numeric(length(s)))
  u <- sort(u)
  hi <- findInterval(s + h, u)                         # #{u <= s + h}
  lo <- findInterval(s - h, u, left.open = TRUE)       # #{u <  s - h}
  as.numeric(hi - lo)
}

# ---------------------------------------------------------------------------
# Greatest common divisor of positive real numbers (e.g. visit gaps measured
# in years), computed by the Euclidean algorithm with a relative tolerance.
# ---------------------------------------------------------------------------
.gcd_numeric <- function(x, rel_tol = 1e-8) {
  x <- sort(unique(as.numeric(x)))
  tol <- rel_tol * max(x)
  gcd2 <- function(a, b) {
    while (b > tol) {
      r <- a %% b
      if (r > b - tol) r <- 0
      a <- b
      b <- r
    }
    a
  }
  g <- x[1L]
  for (xi in x[-1L]) g <- gcd2(xi, g)
  g
}

# Pointwise sd / quantile that tolerate non-finite values.
.safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

.safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else as.numeric(stats::quantile(x, probs = p, names = FALSE))
}

.format_num <- function(x, digits = 4) as.character(signif(x, digits))
