#' Shape functions used in the simulation studies
#'
#' The two shape functions of Section 6 of Hu et al., both probability
#' densities on \eqn{[0, 1]}:
#' * `shape_monotone()`: the monotonically increasing function
#'   \eqn{f(t) = 8/\{11 + 11\exp(5 - 10t)\} + 7/11};
#' * `shape_ushape()`: the U-shaped function
#'   \eqn{f(t) = 4/5 + (60/37)(t - 3/10)^2}.
#'
#' @param t Numeric vector of ages on the unit scale.
#' @return Numeric vector of \eqn{f(t)}.
#' @template ref
#' @examples
#' curve(shape_monotone(x), 0, 1, ylim = c(0.5, 1.6), ylab = "f(t)", xlab = "t")
#' curve(shape_ushape(x), 0, 1, add = TRUE, lty = 2)
#' integrate(shape_monotone, 0, 1)$value
#' integrate(shape_ushape, 0, 1)$value
#' @name shape_functions
NULL

#' @rdname shape_functions
#' @export
shape_monotone <- function(t) 8 / (11 + 11 * exp(5 - 10 * t)) + 7 / 11

#' @rdname shape_functions
#' @export
shape_ushape <- function(t) 4 / 5 + (60 / 37) * (t - 3 / 10)^2

.lognormal_par <- function(mean, var) {
  if (var <= 1e-12) return(list(meanlog = log(mean), sdlog = 0))
  s2 <- log1p(var / mean^2)
  list(meanlog = log(mean) - 0.5 * s2, sdlog = sqrt(s2))
}

.resolve_shape <- function(shape) {
  if (is.function(shape)) return(shape)
  shape <- match.arg(shape, c("monotone", "ushape"))
  if (shape == "monotone") shape_monotone else shape_ushape
}

.check_par <- function(x, name, names_needed) {
  if (!is.numeric(x) || any(!is.finite(x)) || is.null(names(x)) ||
      !all(names_needed %in% names(x))) {
    .stopf("`%s` must be a named numeric vector with elements %s.", name,
           paste(sprintf("'%s'", names_needed), collapse = ", "))
  }
  x
}

#' Simulate longitudinal cohort data with death and delayed entry
#'
#' Generates data from the simulation designs of Section 6 of Hu et al. on the
#' unit age scale \eqn{[0, 1]}. In the underlying population, the latent
#' variable \eqn{Z^0} is log-normal, the age at death \eqn{T^0} given
#' \eqn{Z^0 = z} is Weibull with shape \eqn{t_1} and scale
#' \eqn{t_2 z^{-1/t_1}}, and the marker is \eqn{Y^0(t) = Z^0 f(t) e_t} with
#' independent log-normal noise \eqn{e_t} of mean one. Entry ages \eqn{A^0}
#' follow either
#' * random truncation: \eqn{A^0 \sim U(a_1, a_2)}, independent of
#'   \eqn{(Z^0, T^0)} (no birth-cohort effect), or
#' * informative truncation: \eqn{A^0 + \delta} given \eqn{Z^0 = z} is Weibull
#'   with shape \eqn{a_1'} and scale \eqn{a_2' z^{-1/a_1'}} (birth-cohort
#'   effect).
#'
#' Individuals are enrolled only if alive at entry (\eqn{A^0 \le T^0}); the
#' first `n` enrolled individuals form the sample. Each subject has a baseline
#' visit at \eqn{A} and up to \eqn{k} follow-up visits at
#' \eqn{A + v_1 + \dots + v_l}; follow-up stops at death or after the
#' \eqn{C}-th follow-up visit, where \eqn{C = \min\{k, \lfloor W \rfloor\}}{C = min(k, floor(W))}
#' and \eqn{W} is Weibull, so that \eqn{P(C = l) = F_W(l + 1) - F_W(l)} for
#' \eqn{l < k}.
#'
#' @param n Number of enrolled subjects (the sample size).
#' @param shape Shape function: `"monotone"` ([shape_monotone()]), `"ushape"`
#'   ([shape_ushape()]), or a vectorized function that is positive and
#'   integrates to one on \eqn{[0, 1]}.
#' @param visit_gaps `"equal"` for the equally spaced design of the paper
#'   (\eqn{k = 5} follow-ups with gap \eqn{v = 1/\lfloor 6 n^{0.16}\rfloor}),
#'   `"unequal"` for the unequally spaced design
#'   \eqn{(v_1, \dots, v_5) = (5, 5, 6, 5, 6) \times 0.02}, or a numeric vector
#'   of gaps (its length is the number of follow-up visits \eqn{k}).
#' @param truncation `"random"` or `"informative"`.
#' @param trunc_par Parameters of the entry-age distribution: for random
#'   truncation a named vector `c(min = , max = )` (default `c(-0.3, 1.9)`);
#'   for informative truncation `c(shape = , scale = , shift = )` (default
#'   `c(3.2, 1.2, 0.3)`).
#' @param death_par Weibull parameters `c(shape = , scale = )` of the age at
#'   death (default `c(4.7, 0.85)`).
#' @param dropout_par Weibull parameters `c(shape = , scale = )` of \eqn{W}
#'   (default `c(3, 3)`; the paper uses scale 5 in its small-sample design).
#' @param frailty_par Mean and variance `c(mean = , var = )` of the log-normal
#'   latent variable \eqn{Z^0} (default `c(1, 0.25)`).
#' @param noise_var Variance of the log-normal noise \eqn{e_t}, whose mean is
#'   fixed at one so that \eqn{E\{Y^0(t) \mid Z^0\} = Z^0 f(t)} (default 0.05;
#'   use 0 for noise-free data).
#' @param pop_size Number of individuals drawn from the underlying population
#'   before truncation. The default `max(10000, 4 * n)` gives the population
#'   size used in the paper for \eqn{n \le 2500}.
#' @param seed Seed, or `NULL` to use the current random number stream. The
#'   caller's random number stream is left unchanged when a seed is supplied.
#'
#' @return An object of class `"trajchain_sim"`, a list with components
#'   \describe{
#'     \item{`Y`}{\eqn{n \times (k+1)} matrix of observed marker values, `NA`
#'       after death or dropout;}
#'     \item{`entry_age`, `visit_gaps`, `tau0`, `tau1`}{entry ages, visit gaps
#'       and the age interval \eqn{[0, 1]}, ready to be passed to
#'       [trajchain()];}
#'     \item{`latent`}{data frame with the latent variable `Z`, the age at
#'       death `death_age` and the dropout index `dropout` (\eqn{C});}
#'     \item{`Y_complete`}{the marker at every scheduled visit, observed or
#'       not;}
#'     \item{`shape_fun`}{the true shape function;}
#'     \item{`prop_enrolled`, `mean_visits`, `settings`}{the proportion of the
#'       population alive at entry, the mean number of observed measurements per
#'       subject, and all simulation parameters.}
#'   }
#'
#' @template ref
#' @seealso [sim_truth()] for the true curves of these designs.
#'
#' @examples
#' sim <- simulate_cohort(n = 1000, seed = 1)
#' sim
#' head(sim$Y)
#'
#' # unequally spaced visits with informative truncation
#' sim2 <- simulate_cohort(n = 1000, shape = "ushape", visit_gaps = "unequal",
#'                         truncation = "informative", seed = 2)
#' sim2$visit_gaps
#'
#' @export
simulate_cohort <- function(n = 1000, shape = "monotone", visit_gaps = "equal",
                            truncation = c("random", "informative"),
                            trunc_par = NULL,
                            death_par = c(shape = 4.7, scale = 0.85),
                            dropout_par = c(shape = 3, scale = 3),
                            frailty_par = c(mean = 1, var = 0.25),
                            noise_var = 0.05,
                            pop_size = max(10000, 4 * n), seed = NULL) {
  n <- .check_count(n, "n", min = 2L)
  pop_size <- .check_count(pop_size, "pop_size", min = 2L)
  f <- .resolve_shape(shape)
  truncation <- match.arg(truncation)
  if (is.null(trunc_par)) {
    trunc_par <- if (truncation == "random") c(min = -0.3, max = 1.9) else
      c(shape = 3.2, scale = 1.2, shift = 0.3)
  }
  if (truncation == "random") {
    trunc_par <- .check_par(trunc_par, "trunc_par", c("min", "max"))
    if (trunc_par[["min"]] >= trunc_par[["max"]]) .stopf("`trunc_par`: need min < max.")
  } else {
    trunc_par <- .check_par(trunc_par, "trunc_par", c("shape", "scale", "shift"))
  }
  death_par <- .check_par(death_par, "death_par", c("shape", "scale"))
  dropout_par <- .check_par(dropout_par, "dropout_par", c("shape", "scale"))
  frailty_par <- .check_par(frailty_par, "frailty_par", c("mean", "var"))
  if (frailty_par[["mean"]] <= 0 || frailty_par[["var"]] < 0) {
    .stopf("`frailty_par`: the mean must be positive and the variance non-negative.")
  }
  noise_var <- .check_number(noise_var, "noise_var", nonnegative = TRUE)

  if (is.character(visit_gaps)) {
    visit_gaps <- match.arg(visit_gaps, c("equal", "unequal"))
    design <- visit_gaps
    gaps <- if (design == "equal") rep(1 / floor(6 * (n^0.16)), 5L) else c(5, 5, 6, 5, 6) * 0.02
  } else {
    if (!is.numeric(visit_gaps) || !length(visit_gaps) || any(!is.finite(visit_gaps)) ||
        any(visit_gaps <= 0)) {
      .stopf("`visit_gaps` must be \"equal\", \"unequal\" or a vector of positive numbers.")
    }
    design <- "custom"
    gaps <- as.numeric(visit_gaps)
  }
  k <- length(gaps)
  equal <- all(gaps == gaps[1L])
  v <- gaps[1L]
  time_visit <- if (equal) v * (0:k) else cumsum(c(0, gaps))

  seed <- .check_seed(seed)
  if (!is.null(seed)) {
    old <- .get_rng_state()
    on.exit(.restore_rng_state(old), add = TRUE)
    set.seed(seed)
  }

  # 1) population: latent variable, entry age and age at death
  pz <- .lognormal_par(frailty_par[["mean"]], frailty_par[["var"]])
  Z0 <- stats::rlnorm(pop_size, meanlog = pz$meanlog, sdlog = pz$sdlog)
  if (truncation == "random") {
    A0 <- stats::runif(pop_size, min = trunc_par[["min"]], max = trunc_par[["max"]])
  } else {
    a1 <- trunc_par[["shape"]]
    A0 <- stats::rweibull(pop_size, shape = a1,
                          scale = trunc_par[["scale"]] * Z0^(-1 / a1)) - trunc_par[["shift"]]
  }
  t1 <- death_par[["shape"]]
  T0 <- stats::rweibull(pop_size, shape = t1, scale = death_par[["scale"]] * Z0^(-1 / t1))

  # 2) enrollment: alive at entry
  idx_all <- which(A0 <= T0)
  if (length(idx_all) < n) {
    .stopf(paste0(
      "Only %d of the %d simulated individuals are alive at entry, fewer than n = %d; ",
      "increase `pop_size`."), length(idx_all), pop_size, n)
  }
  sel <- idx_all[seq_len(n)]
  Z <- Z0[sel]
  A <- A0[sel]
  death_age <- T0[sel]

  # 3) dropout index C in {0, ..., k}
  C <- pmin(k, pmax(0L, as.integer(floor(
    stats::rweibull(n, shape = dropout_par[["shape"]], scale = dropout_par[["scale"]])))))

  # 4) marker at all scheduled visits
  age_mat <- matrix(rep(A, times = k + 1L), nrow = n, ncol = k + 1L) +
    matrix(rep(time_visit, each = n), nrow = n, ncol = k + 1L)
  F_mat <- f(age_mat)
  if (length(F_mat) != length(age_mat) || any(!is.finite(F_mat))) {
    .stopf("`shape` must be a vectorized function returning finite values.")
  }
  F_mat <- matrix(F_mat, nrow = n, ncol = k + 1L)
  Z_mat <- matrix(Z, nrow = n, ncol = k + 1L)
  pe <- .lognormal_par(1, noise_var)
  E_mat <- matrix(stats::rlnorm(n * (k + 1L), meanlog = pe$meanlog, sdlog = pe$sdlog),
                  nrow = n, ncol = k + 1L)
  Y_complete <- F_mat * Z_mat * E_mat

  # 5) observation: visit l is observed iff l <= C and it occurs before death
  T_tilde <- pmin(A + time_visit[C + 1L], death_age)
  observed <- age_mat <= matrix(T_tilde, nrow = n, ncol = k + 1L)
  observed[, 1L] <- TRUE
  Y <- Y_complete
  Y[!observed] <- NA_real_
  colnames(Y) <- colnames(Y_complete) <- paste0("visit_", 0:k)

  structure(list(
    Y = Y, entry_age = A, visit_gaps = gaps, tau0 = 0, tau1 = 1,
    latent = data.frame(Z = Z, death_age = death_age, dropout = C),
    Y_complete = Y_complete,
    shape_fun = f,
    prop_enrolled = length(idx_all) / pop_size,
    mean_visits = mean(rowSums(observed)),
    settings = list(n = n, design = design, truncation = truncation,
                    trunc_par = trunc_par, death_par = death_par,
                    dropout_par = dropout_par, frailty_par = frailty_par,
                    noise_var = noise_var, pop_size = pop_size, seed = seed)
  ), class = "trajchain_sim")
}

#' True curves of the simulation designs
#'
#' Computes, by numerical integration over the log-normal latent variable, the
#' true shape function, \eqn{\mu_{Z^0}(t) = E(Z^0 \mid A^0 = t, T^0 \ge t)}{mu(t) = E(Z0 | A0 = t, T0 >= t)},
#' the birth-cohort-survivor bias \eqn{\mu_{Z^0}(t)/\mu_{Z^0}(0)} and the
#' unconditional mean \eqn{\mu_{Z^0}(0) f(t)} for the designs of
#' [simulate_cohort()]. Under random truncation \eqn{\mu_{Z^0}(t) = E(Z^0 \mid T^0 \ge t)},
#' the pure survivor bias, and \eqn{\mu_{Z^0}(0) = E(Z^0)}.
#'
#' @inheritParams simulate_cohort
#' @param t Ages on the unit scale at which to evaluate the curves.
#'
#' @return A data frame with columns `t`, `shape`, `mu`, `bias` and `mean`.
#' @template ref
#' @examples
#' truth <- sim_truth(seq(0, 1, by = 0.1))
#' truth
#'
#' # informative truncation: birth-cohort and survivor effects combined
#' sim_truth(c(0, 0.5, 1), truncation = "informative")
#' @export
sim_truth <- function(t, shape = "monotone", truncation = c("random", "informative"),
                      trunc_par = NULL,
                      death_par = c(shape = 4.7, scale = 0.85),
                      frailty_par = c(mean = 1, var = 0.25)) {
  if (!is.numeric(t) || any(!is.finite(t))) .stopf("`t` must be a vector of finite numbers.")
  f <- .resolve_shape(shape)
  truncation <- match.arg(truncation)
  if (is.null(trunc_par) && truncation == "informative") {
    trunc_par <- c(shape = 3.2, scale = 1.2, shift = 0.3)
  }
  if (truncation == "informative") {
    trunc_par <- .check_par(trunc_par, "trunc_par", c("shape", "scale", "shift"))
  }
  death_par <- .check_par(death_par, "death_par", c("shape", "scale"))
  frailty_par <- .check_par(frailty_par, "frailty_par", c("mean", "var"))
  pz <- .lognormal_par(frailty_par[["mean"]], frailty_par[["var"]])

  # E[g(Z)] for log-normal Z
  expect <- function(g) {
    if (pz$sdlog == 0) return(g(exp(pz$meanlog)))
    lo <- pz$meanlog - 12 * pz$sdlog
    hi <- pz$meanlog + 12 * pz$sdlog
    stats::integrate(function(x) g(exp(x)) * stats::dnorm(x, pz$meanlog, pz$sdlog),
                     lo, hi, rel.tol = 1e-10, subdivisions = 500L)$value
  }
  # Given Z = z: P(T0 >= t) = exp(-z (t / scale)^shape); under informative
  # truncation the entry-age density at t is proportional to
  # z exp(-z ((t + shift) / scale_A)^shape_A).
  mu_at <- function(tt) {
    sT <- if (tt > 0) (tt / death_par[["scale"]])^death_par[["shape"]] else 0
    if (truncation == "random") {
      num <- expect(function(z) z * exp(-z * sT))
      den <- expect(function(z) exp(-z * sT))
    } else {
      x <- tt + trunc_par[["shift"]]
      if (x <= 0) return(NA_real_)
      sA <- (x / trunc_par[["scale"]])^trunc_par[["shape"]]
      num <- expect(function(z) z^2 * exp(-z * (sT + sA)))
      den <- expect(function(z) z * exp(-z * (sT + sA)))
    }
    num / den
  }
  mu <- vapply(t, mu_at, numeric(1))
  mu0 <- mu_at(0)
  data.frame(t = t, shape = f(t), mu = mu, bias = mu / mu0, mean = mu0 * f(t))
}
