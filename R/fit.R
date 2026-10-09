#' Estimate marker trajectories in the presence of death and delayed entry
#'
#' Fits the nonparametric estimators of Hu, Wu, Zhao and Sun for longitudinal
#' cohort studies in which individuals enter at different ages (delayed entry,
#' i.e. left truncation) and are followed until death, dropout or the end of
#' the visit schedule. The function estimates
#' * the **shape function** \eqn{f(t)}, the common age-related pattern of the
#'   trajectory;
#' * \eqn{\mu_{Z^0}(t) = E(Z^0 \mid A^0 = t, T^0 \ge t)}{mu(t) = E(Z0 | A0 = t, T0 >= t)},
#'   the mean of the latent variable among survivors entering at age \eqn{t};
#' * the **birth-cohort-survivor bias**
#'   \eqn{\mu_{Z^0}(t)/\mu_{Z^0}(\tau_0)}{mu(t) / mu(tau0)};
#' * the **unconditional mean** trajectory
#'   \eqn{\mu_{Z^0}(\tau_0) f(t)}{mu(tau0) f(t)}, the mean marker level in the
#'   absence of death for the reference subgroup entering at the starting age;
#' * naive Nadaraya--Watson curves, reported for comparison: `nw_baseline`
#'   smooths the baseline measurements against the entry age and estimates
#'   \eqn{E\{Y(A) \mid A = t\} = \mu_{Z^0}(t) f(t)}, and `nw_pooled` smooths all
#'   observed measurements against age. Both mix the shape with survivor and
#'   birth-cohort selection.
#'
#' @section Model and assumptions:
#' Let \eqn{\tau_0} be the starting age of the target population (\eqn{t = 0}
#' in the paper; age 30 in the ABC-DS application). For an individual alive at
#' \eqn{\tau_0}, let \eqn{Y^0(t)} be the marker at age \eqn{t}, \eqn{T^0} the age
#' at death, \eqn{A^0} the age at potential study entry and \eqn{Z^0 > 0} a
#' latent variable. The model is
#' \deqn{E\{Y^0(t) \mid Z^0, T^0 = s\} = Z^0 f(t), \qquad \tau_0 \le t \le s,}{
#'   E{Y0(t) | Z0, T0 = s} = Z0 f(t),  tau0 <= t <= s,}
#' with \eqn{\int_{\tau_0}^{\tau_1} f(t)\,dt = 1}. Individuals are sampled only
#' if alive at entry (\eqn{A^0 \le T^0}), and \eqn{C} denotes the last follow-up
#' visit before loss to follow-up for reasons other than death. The estimators
#' rely on
#' * (A1) \eqn{A^0} is independent of \eqn{\{Y^0(t), t \le s\}} given
#'   \eqn{Z^0} and \eqn{T^0 = s}: the entry age may depend on the latent
#'   variable (birth-cohort effects), but not otherwise on the marker process;
#' * (A2) \eqn{C} is independent of \eqn{(T, Y(\cdot))} given \eqn{(A, Z)}:
#'   dropout may depend on the entry age and the latent variable, but not on the
#'   marker values themselves or on the time of death beyond what \eqn{Z}
#'   explains.
#'
#' The distribution of \eqn{Z^0} and its association with death, entry age and
#' dropout are left unspecified. Under (A1)--(A2) these nuisance factors cancel
#' in the ratios \eqn{f(t)/f(t - v)}, which are therefore identified from pairs
#' of consecutive measurements. Dropout that depends on the marker values (for
#' example, sicker participants leaving the study) violates (A2); the paper's
#' sensitivity analysis found the estimators stable under a mild violation and
#' degraded under a severe one.
#'
#' @section Data:
#' `Y` holds one row per subject and one column per scheduled visit, the first
#' column being the baseline measurement at entry age \eqn{A_i}; visit \eqn{l}
#' is scheduled at age \eqn{A_i + v_1 + \dots + v_l}. Unobserved visits (after
#' death or dropout) are `NA`. Survival times and the reason for dropout are
#' not needed. Subjects may enter before \eqn{\tau_0}; all observed visits are
#' used. Use [visit_matrix()] to build `Y` from long-format data.
#'
#' @section Estimation of the shape function:
#' The age interval \eqn{[\tau_0, \tau_1]} is divided into \eqn{m} subintervals
#' of length \eqn{\delta = (\tau_1 - \tau_0)/m}, whose midpoints
#' \eqn{u_j = \tau_0 + (j - 0.5)\delta} form the grid.
#' * `method = "chain"` (Section 3 of the paper) requires equally spaced visits
#'   with \eqn{\delta} equal to the visit gap. The ratios
#'   \eqn{\hat r(u_j, u_{j-1})} of [ratio_estimate()] are chained,
#'   \deqn{\hat f(u_j) = \frac{\prod_{k=2}^{j} \hat r(u_k, u_{k-1})}
#'     {\delta \{\sum_{q=2}^{m} \prod_{k=2}^{q} \hat r(u_k, u_{k-1}) + 1\}},}{
#'     f(u_j) = prod_{k=2..j} r(u_k, u_{k-1}) / [delta {sum_{q=2..m} prod_{k=2..q} r(u_k, u_{k-1}) + 1}],}
#'   so that the midpoint Riemann sum of \eqn{\hat f} equals one.
#' * `method = "penalized"` (Section 5) handles widely spaced visits (a grid
#'   finer than the visit gap) and unequally spaced visits (grid step equal to
#'   the greatest common divisor of the gaps). For every gap \eqn{v_l} and every
#'   pair of grid points with \eqn{u_j - u_{j'} = v_l}, the ratio
#'   \eqn{\hat r(u_j, u_{j'})} is estimated; the grid values \eqn{g(u_j)}
#'   minimize
#'   \deqn{\sum_{(j, j')} \hat w_{j j'} \{g(u_j) - \hat r(u_j, u_{j'}) g(u_{j'})\}^2
#'     + \lambda m^3 \sum_{j=2}^{m-1} \{g(u_{j-1}) - 2 g(u_j) + g(u_{j+1})\}^2,}{
#'     sum w_jj' {g(u_j) - r(u_j, u_j') g(u_j')}^2 + lambda m^3 sum {g(u_{j-1}) - 2 g(u_j) + g(u_{j+1})}^2,}
#'   where \eqn{\hat w_{jj'} = n_{jj'} I(n_{jj'} \ge c^*) / \sum n I(n \ge c^*)},
#'   \eqn{n_{jj'}} is the number of observations contributing to the ratio and
#'   \eqn{c^* = } `threshold` \eqn{\times\, n h}. Both \eqn{h} in \eqn{c^*} and
#'   the penalty are computed on the unit time scale
#'   \eqn{(t - \tau_0)/(\tau_1 - \tau_0)}, so that `lambda` is dimensionless and
#'   comparable across studies. Ratios whose kernel window contains no observed
#'   pair are undefined and are dropped. With `constraint = "anchor"` (the
#'   authors' implementation) the quadratic form is minimized with \eqn{g(u_1)}
#'   fixed and the solution is rescaled to integrate to one; with
#'   `constraint = "integral"` it is solved as a quadratic program under the
#'   constraints \eqn{\delta \sum_j g(u_j) = 1} and \eqn{g(u_j) \ge 0}, as written
#'   in Section 5 of the paper. With equally spaced visits, a grid step equal
#'   to the visit gap, `lambda = 0` and `threshold = 0`, the anchored penalized
#'   estimator reduces to ratio chaining.
#'
#' In both cases \eqn{\hat f} is extended to \eqn{[\tau_0, \tau_1]} by linear
#' interpolation between grid points and linear extrapolation of the first and
#' last segments.
#'
#' @section Estimation of the mean trajectories:
#' With \eqn{\hat f} available, \eqn{\mu_{Z^0}(t)} is estimated by
#' * `mu_method = "baseline"`: equation (5),
#'   \eqn{\hat\mu_{Z^0}(t) = \sum_i K_{h'}(t - A_i) Y_i(A_i) / \{\hat f(t) \sum_i K_{h'}(t - A_i)\}},
#'   which uses the baseline measurements only and is valid under (A1)--(A2),
#'   i.e. in the presence of birth-cohort effects (used for the ABC-DS analysis);
#' * `mu_method = "pooled"`: equation (6), the same estimator pooling all
#'   observed visits. It requires the stronger assumption (A2') that dropout is
#'   independent of \eqn{(T, Y(\cdot), Z)} given \eqn{A}, and it approximates
#'   \eqn{\mu_{Z^0}(t)} when the visit gaps are small; in the presence of
#'   birth-cohort effects and wider gaps it mixes subjects who entered at
#'   different ages. It was used in the simulations of the paper, which have no
#'   birth-cohort effect.
#'
#' Both versions are always computed and stored; `mu_method` only selects the
#' one shown by default in `print()`, `summary()` and `plot()`. The unconditional
#' mean and the bias curve are anchored at \eqn{\tau_0} through
#' \eqn{\hat\mu_{Z^0}(\tau_0)}, so they inherit the uncertainty of the
#' estimates at the starting age; `trajchain()` warns when \eqn{\hat f} is not
#' positive somewhere in \eqn{[\tau_0, \tau_1]}, when
#' \eqn{\hat\mu_{Z^0}(\tau_0)} cannot be computed, or when fewer than 10
#' measurements lie within one bandwidth \eqn{h'} of \eqn{\tau_0}.
#'
#' @section Tuning parameters:
#' By default the ratio bandwidth is \eqn{h = \hat\sigma_A n^{-1/6}} (the
#' "small bandwidth" of the paper) and the mean bandwidth is
#' \eqn{h' = 2.34\,\hat\sigma_A n^{-1/5}}, where \eqn{\hat\sigma_A} is the sample
#' standard deviation of the entry ages; these are the bandwidths of the
#' authors' simulation and data-analysis code. Use `bandwidth = "cv"` (or the
#' output of [cv_bandwidth()]) to choose the constant \eqn{\kappa} in
#' \eqn{h = \kappa \hat\sigma_A n^{-1/6}} by five-fold cross-validation, and
#' `lambda = "cv"` (or the output of [cv_lambda()]) to choose the penalty of
#' the penalized estimator.
#'
#' @inheritParams ratio_estimate
#' @param visit_gaps Scheduled time between consecutive visits, on the same
#'   scale as `entry_age`: a single number for equally spaced visits, or a
#'   vector \eqn{(v_1, \dots, v_k)} of length `ncol(Y) - 1`. Use the scheduled
#'   (rounded) gaps, not observed average gaps.
#' @param tau0 Starting age \eqn{\tau_0} of the trajectory (the target
#'   population consists of individuals alive at this age).
#' @param tau1,m End of the age interval \eqn{\tau_1}, or the number of grid
#'   points \eqn{m}; supply one of them. \eqn{\tau_1 - \tau_0} must be an
#'   integer multiple of the grid step.
#' @param grid_step Spacing \eqn{\delta} of the grid. The default is the
#'   greatest common divisor of `visit_gaps` (the visit gap itself when visits
#'   are equally spaced). Use a value \eqn{v/q} to refine the grid for widely
#'   spaced visits (Section 5 of the paper).
#' @param method `"chain"`, `"penalized"`, or `"auto"` (the default), which uses
#'   `"chain"` whenever every visit gap equals the grid step and `"penalized"`
#'   otherwise.
#' @param mu_method Estimator of \eqn{\mu_{Z^0}(t)} reported by default:
#'   `"baseline"` (equation 5) or `"pooled"` (equation 6). See the section
#'   *Estimation of the mean trajectories*.
#' @param bandwidth Bandwidth \eqn{h} of the ratio estimator on the age scale:
#'   `"rot"` (default, \eqn{\hat\sigma_A n^{-1/6}}), `"cv"` (constant chosen
#'   by [cv_bandwidth()] with its default settings), a `trajchain_cv` object
#'   returned by [cv_bandwidth()] (to control the cross-validation), or a
#'   positive number.
#' @param bandwidth_mu Bandwidth \eqn{h'} of the mean estimators: `"rot"`
#'   (default, \eqn{2.34\,\hat\sigma_A n^{-1/5}}) or a positive number.
#' @param lambda Penalty \eqn{\lambda} of the penalized estimator: `"cv"`
#'   (default; five-fold cross-validation over the grid used in the paper, see
#'   [cv_lambda()]), a `trajchain_cv` object returned by [cv_lambda()], or a
#'   non-negative number. Ignored by the chain estimator.
#' @param threshold The constant in \eqn{c^* = \mathrm{threshold} \times n h},
#'   where \eqn{h} is the ratio bandwidth on the unit time scale,
#'   \eqn{h/(\tau_1 - \tau_0)}: ratio estimates based on fewer than \eqn{c^*}
#'   observations receive zero weight in the penalized estimator. The paper
#'   uses 0.6.
#' @param constraint How the penalized problem is normalized: `"anchor"`
#'   (default; the solution used for the results in the paper) or `"integral"`
#'   (quadratic program with the integral and non-negativity constraints).
#' @param n_grid Number of equally spaced ages in \eqn{[\tau_0, \tau_1]} at which
#'   the curves are evaluated (default 801).
#' @param keep_data Keep a copy of the data in the fitted object. Required by
#'   [trajchain_boot()]; [predict.trajchain()] uses it to evaluate the mean
#'   curves exactly (otherwise it interpolates `curves`).
#'
#' @return An object of class `"trajchain"`, a list with components
#'   \describe{
#'     \item{`curves`}{data frame with the estimates on `n_grid` equally spaced
#'       ages: `t`, `shape` (\eqn{\hat f}), and for each estimator of
#'       \eqn{\mu_{Z^0}} (suffix `_baseline` or `_pooled`) the unconditional mean
#'       `mean_*`, the birth-cohort-survivor bias `bias_*`, \eqn{\hat\mu_{Z^0}(t)}
#'       `mu_*`, and the Nadaraya--Watson curve `nw_*`;}
#'     \item{`grid`}{data frame with the grid points `t`, the estimated shape
#'       function `shape` and, for the chain estimator, the chained ratios
#'       `ratio` \eqn{= \hat r(u_j, u_{j-1})};}
#'     \item{`pairs`}{for the penalized estimator, a data frame describing every
#'       ratio constraint (ages, gap, estimated ratio, number of contributing
#'       observations `n`, whether it was kept, and its weight);}
#'     \item{`mu0`}{\eqn{\hat\mu_{Z^0}(\tau_0)} for both estimators;}
#'     \item{`shape_fun`}{a function returning \eqn{\hat f(t)} for
#'       \eqn{t \in [\tau_0, \tau_1]} (`NA` outside);}
#'     \item{`bandwidth`, `lambda`, `threshold`, `method`, `mu_method`,
#'       `constraint`, `kernel`}{the settings used;}
#'     \item{`design`}{a summary of the design (interval, grid, visit gaps,
#'       numbers of subjects and measurements);}
#'     \item{`cv`}{the cross-validation objects, when `"cv"` was requested;}
#'     \item{`data`, `settings`, `call`}{the data (if `keep_data = TRUE`), the
#'       internal settings used for refitting, and the call.}
#'   }
#'
#' @template ref
#' @seealso [trajchain_boot()] for bootstrap standard errors and confidence
#'   bands, [cv_bandwidth()] and [cv_lambda()] for tuning, [plot.trajchain()]
#'   and [plot_trajectories()] for graphics, [simulate_cohort()] for data from
#'   the simulation designs of the paper.
#'
#' @examples
#' # Equally spaced visits (Section 3): simulated data from the paper's design
#' sim <- simulate_cohort(n = 1000, shape = "monotone", seed = 2026)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  mu_method = "pooled")
#' fit
#' plot(fit)
#'
#' # compare with the true shape function
#' curve(shape_monotone(x), 0, 1, lty = 2, ylab = "f(t)", xlab = "t")
#' lines(fit$curves$t, fit$curves$shape, col = "#2a78d6", lwd = 2)
#'
#' \donttest{
#' # Unequally spaced visits (Section 5), penalty chosen by cross-validation
#' sim3 <- simulate_cohort(n = 1000, visit_gaps = "unequal", seed = 7)
#' fit3 <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, tau0 = 0, tau1 = 1,
#'                   mu_method = "pooled")
#' fit3$lambda
#' plot(fit3, which = c("shape", "mean"))
#' }
#'
#' @export
trajchain <- function(Y, entry_age, visit_gaps, tau0, tau1 = NULL, m = NULL,
                      grid_step = NULL,
                      method = c("auto", "chain", "penalized"),
                      mu_method = c("baseline", "pooled"),
                      bandwidth = "rot", bandwidth_mu = "rot",
                      lambda = "cv", threshold = 0.6,
                      constraint = c("anchor", "integral"),
                      kernel = "epanechnikov", n_grid = 801, keep_data = TRUE) {
  cl <- match.call()
  method <- match.arg(method)
  mu_method <- match.arg(mu_method)
  constraint <- match.arg(constraint)
  .check_flag(keep_data, "keep_data")
  n_grid <- .check_count(n_grid, "n_grid", min = 2L)
  kern <- kernel_function(kernel)
  design <- .prepare_design(Y, entry_age, visit_gaps, tau0, tau1 = tau1, m = m,
                            grid_step = grid_step)

  chain_ok <- all(design$gap_units == 1L)
  if (method == "auto") method <- if (chain_ok) "chain" else "penalized"
  if (method == "chain" && !chain_ok) {
    .stopf(paste0(
      "method = \"chain\" requires every visit gap to equal the grid step (%s); ",
      "the gaps are %s grid steps. Use method = \"penalized\"."),
      .format_num(design$grid_step), paste(design$gap_units, collapse = ", "))
  }

  bw <- .resolve_ratio_bandwidth(bandwidth, design, kern)
  h_mu <- .resolve_mean_bandwidth(bandwidth_mu, design)

  lam <- NA_real_
  cv_lam <- NULL
  if (method == "penalized") {
    if (design$m < 3L) .stopf("The penalized estimator needs at least three grid points (m >= 3).")
    if (design$m > 2000L) {
      .stopf("The penalized estimator with m = %d grid points is too large; use a coarser `grid_step`.",
             design$m)
    }
    threshold <- .check_number(threshold, "threshold", nonnegative = TRUE)
    if (inherits(lambda, "trajchain_cv")) {
      if (!identical(lambda$type, "lambda")) {
        .stopf("`lambda` must be a cross-validation object returned by cv_lambda().")
      }
      lam <- lambda$lambda
      cv_lam <- lambda
    } else if (is.character(lambda)) {
      match.arg(lambda, "cv")
      cv_lam <- .cv_lambda_core(design, h = bw$h, threshold = threshold,
                                constraint = constraint, kern = kern)
      lam <- cv_lam$lambda
    } else {
      lam <- .check_number(lambda, "lambda", nonnegative = TRUE)
    }
  } else if (!(is.character(lambda) && identical(lambda, "cv"))) {
    .warnf("`lambda` is ignored by the chain estimator.")
  }

  settings <- list(visit_gaps = design$visit_gaps, tau0 = design$tau0,
                   tau1 = design$tau1, grid_step = design$grid_step, m = design$m,
                   method = method, h = bw$h, h_mu = h_mu, lambda = lam,
                   threshold = if (method == "penalized") threshold else NA_real_,
                   constraint = constraint, kern = kern, n_grid = n_grid)
  core <- .fit_core(design, settings)
  .warn_anchor(core, design, mu_method, h_mu)

  grid <- data.frame(t = design$t_grid, shape = core$phi / design$tau)
  if (method == "chain") grid$ratio <- core$shape$ratio
  pairs <- if (method == "penalized") {
    tab <- core$shape$table
    data.frame(age_from = design$t_grid[tab$i1], age_to = design$t_grid[tab$i2],
               gap = tab$gap * design$grid_step, ratio = tab$num / tab$den,
               n = tab$n, kept = tab$kept, weight = tab$weight)
  } else {
    NULL
  }

  structure(list(
    curves = core$curves,
    grid = grid,
    pairs = pairs,
    mu0 = core$mu0,
    shape_fun = .make_shape_fun(core$phi_fun, design$tau0, design$tau1),
    method = method,
    mu_method = mu_method,
    bandwidth = c(ratio = bw$h, mean = h_mu),
    bandwidth_constant = bw$constant,
    lambda = lam,
    threshold = settings$threshold,
    constraint = if (method == "penalized") constraint else NA_character_,
    kernel = .kernel_label(kernel),
    design = list(tau0 = design$tau0, tau1 = design$tau1, m = design$m,
                  grid_step = design$grid_step, visit_gaps = design$visit_gaps,
                  n = design$n, n_obs = design$n_obs,
                  n_visits_max = design$p, n_intermittent = design$n_intermittent),
    cv = list(bandwidth = bw$cv, lambda = cv_lam),
    data = if (keep_data) list(Y = design$Y, entry_age = design$entry_age) else NULL,
    settings = settings,
    call = cl
  ), class = "trajchain")
}

# ---------------------------------------------------------------------------
# Core fitting routine shared by trajchain(), the bootstrap and the CV code.
# All quantities are numeric here (bandwidths on the age scale).
# ---------------------------------------------------------------------------
.fit_core <- function(design, settings) {
  h_u <- settings$h / design$tau
  if (settings$method == "chain") {
    sh <- .shape_chain(design, h_u, settings$kern)
  } else {
    sys <- .penalized_system(design, h_u, settings$threshold, settings$kern)
    sh <- .penalized_solve(sys, design$m, settings$lambda, settings$constraint)
    sh$table <- sys$table
  }
  mc <- .mean_curves(design, sh$phi, settings$h_mu / design$tau, settings$kern,
                     settings$n_grid)
  list(phi = sh$phi, shape = sh, curves = mc$curves, mu0 = mc$mu0, phi_fun = mc$phi_fun)
}

# Ratio chaining (Section 3). Returns the shape on the unit scale (phi) and the
# chained ratios.
.shape_chain <- function(design, h_u, kern) {
  pg <- design$pairs[["1"]]
  if (is.null(pg) || !length(pg$u)) {
    .stopf("No pair of consecutive observed visits is available to estimate the ratios.")
  }
  m <- design$m
  s_eval <- design$s_grid[-1L]
  S <- .ksum(s_eval, pg$u, cbind(pg$ycur, pg$ylag), h_u, kern)
  num <- S[, 1L]
  den <- S[, 2L]
  ratio <- num / den
  ratio[den == 0] <- NA_real_
  bad <- which(!is.finite(ratio))
  if (length(bad)) {
    ages <- design$t_grid[-1L][bad]
    .stopf(paste0(
      "The ratio estimate is undefined at grid age(s) %s: no pair of consecutive ",
      "observed visits lies within one bandwidth (h = %s) of these ages. Increase ",
      "`bandwidth` or narrow the age interval [tau0, tau1]."),
      paste(.format_num(utils::head(ages, 5)), collapse = ", "),
      .format_num(h_u * design$tau))
  }
  cp <- cumprod(ratio)
  D <- (1 / m) * (sum(cp) + 1)
  phi <- c(1, cp) / D
  list(phi = phi, ratio = c(NA_real_, ratio))
}

# Penalized estimator (Section 5): the weighted least-squares system.
.penalized_system <- function(design, h_u, threshold, kern) {
  m <- design$m
  s_grid <- design$s_grid
  tabs <- list()
  for (key in names(design$pairs)) {
    pg <- design$pairs[[key]]
    g <- pg$gap
    if (g >= m) next
    i1 <- seq_len(m - g)
    i2 <- i1 + g
    s_r <- s_grid[i2]
    S <- .ksum(s_r, pg$u, cbind(pg$ycur, pg$ylag), h_u, kern)
    tabs[[key]] <- data.frame(i1 = i1, i2 = i2, gap = g, num = S[, 1L], den = S[, 2L],
                              n = .count_window(s_r, pg$u, h_u))
  }
  if (!length(tabs)) {
    .stopf("No ratio between grid points can be estimated: every visit gap exceeds the age interval.")
  }
  tab <- do.call(rbind, tabs)
  rownames(tab) <- NULL
  cstar <- threshold * design$n * h_u
  keep <- tab$n >= cstar
  W_sum <- sum(tab$n[keep])
  if (!is.finite(W_sum) || W_sum <= 0) {
    .stopf(paste0(
      "No ratio estimate is based on at least c* = threshold * n * h = %s observations; ",
      "decrease `threshold` or increase `bandwidth`."), .format_num(cstar))
  }
  w <- ifelse(keep, sqrt(tab$n / W_sum), 0)
  Amat <- matrix(0, nrow = nrow(tab), ncol = m)
  rr <- seq_len(nrow(tab))
  Amat[cbind(rr, tab$i1)] <- -(w * tab$num) / tab$den
  Amat[cbind(rr, tab$i2)] <- w
  ok <- rowSums(!is.finite(Amat)) == 0L
  tab$kept <- keep & ok
  tab$weight <- ifelse(tab$kept, w^2, 0)
  Amat <- Amat[ok, , drop = FALSE]
  if (!nrow(Amat)) .stopf("All ratio estimates are undefined; increase `bandwidth`.")
  list(Amat = Amat, table = tab, cstar = cstar, W_sum = W_sum)
}

.second_difference <- function(m) {
  D <- matrix(0, nrow = m - 2L, ncol = m)
  idx <- seq_len(m - 2L)
  D[cbind(idx, idx)] <- 1
  D[cbind(idx, idx + 1L)] <- -2
  D[cbind(idx, idx + 2L)] <- 1
  D
}

# Solve the penalized problem for one lambda. Returns phi on the unit scale
# (mean one over the grid). With constraint = "anchor", g' H g is minimized
# subject to g_1 = 1 (closed form g = H^{-1} e_1 / (H^{-1})_{11}, computed from
# the reduced system so that lambda = 0 is allowed whenever the ratios link all
# grid points) and g is then rescaled to integrate to one.
.penalized_solve <- function(sys, m, lambda, constraint = "anchor") {
  D <- .second_difference(m)
  H <- crossprod(sys$Amat) + (lambda * (m^3)) * crossprod(D)
  if (constraint == "anchor") {
    rest <- tryCatch(solve(H[-1L, -1L, drop = FALSE], -H[-1L, 1L]),
                     error = function(e) NULL)
    x <- c(1, rest)
    if (is.null(rest) || any(!is.finite(x)) || sum(x) == 0) {
      .stopf(paste0(
        "The penalized problem is singular: the estimable ratios do not link all ",
        "grid points. Use a positive `lambda` or a larger `bandwidth`."))
    }
    phi <- as.numeric(m * x / sum(x))
  } else {
    scale <- mean(diag(H))
    sol <- tryCatch(
      quadprog::solve.QP(Dmat = 2 * H / scale, dvec = rep(0, m),
                         Amat = cbind(rep(1, m), diag(m)), bvec = c(m, rep(0, m)),
                         meq = 1L),
      error = function(e) NULL)
    if (is.null(sol)) {
      .stopf("The quadratic program could not be solved; use a positive `lambda` or a larger `bandwidth`.")
    }
    phi <- pmax(sol$solution, 0)
    phi <- m * phi / sum(phi)
  }
  list(phi = phi)
}

# Linear interpolation on [0, 1] with linear extrapolation of the end segments.
.interp_unit <- function(x, y) {
  L <- length(x)
  sl <- (y[2L] - y[1L]) / (x[2L] - x[1L])
  y0 <- y[1L] + sl * (0 - x[1L])
  sr <- (y[L] - y[L - 1L]) / (x[L] - x[L - 1L])
  y1 <- y[L] + sr * (1 - x[L])
  stats::approxfun(x = c(0, x, 1), y = c(y0, y, y1), method = "linear", rule = 2,
                   ties = "ordered")
}

# Warnings about the quantities anchored at tau0 (unconditional mean and bias).
.warn_anchor <- function(core, design, mu_method, h_mu) {
  shape <- core$curves$shape
  if (any(!is.finite(shape) | shape <= 0)) {
    .warnf(paste0(
      "The estimated shape function is not positive at some ages in [tau0, tau1] ",
      "(possibly through linear extrapolation beyond the first or last grid point); ",
      "the mean curves are NA there."))
  }
  if (!is.finite(core$mu0[[mu_method]])) {
    .warnf(paste0(
      "mu(tau0) could not be estimated with mu_method = \"%s\", so the unconditional mean ",
      "and the birth-cohort-survivor bias are NA. Check the shape estimate near tau0 and ",
      "the data available near tau0."), mu_method)
  } else {
    src <- if (mu_method == "baseline") design$baseline$u else design$pooled$u
    near <- sum(abs(src) * design$tau <= h_mu)
    if (near < 10) {
      .warnf(paste0(
        "Only %d %s within one bandwidth (h' = %s) of tau0 = %s; the unconditional mean ",
        "and the birth-cohort-survivor bias, which are anchored at tau0, may be unreliable. ",
        "Consider a larger tau0."), near,
        if (mu_method == "baseline") "subjects entered the study" else "measurements lie",
        .format_num(h_mu), .format_num(design$tau0))
    }
  }
  invisible(NULL)
}

.make_shape_fun <- function(phi_fun, tau0, tau1) {
  force(phi_fun)
  tau <- tau1 - tau0
  function(t) {
    out <- phi_fun((t - tau0) / tau) / tau
    out[t < tau0 | t > tau1] <- NA_real_
    out
  }
}

# Nadaraya-Watson numerator/denominator -> NW curve, mu, bias and mean.
.mu_from_sums <- function(S, f) {
  eps <- 1e-12
  num <- S[, 1L]
  den <- S[, 2L]
  den_ok <- is.finite(den) & den > eps
  f_ok <- is.finite(f) & f > eps
  nw <- rep(NA_real_, length(num))
  nw[den_ok] <- num[den_ok] / den[den_ok]
  mu <- rep(NA_real_, length(num))
  ok <- den_ok & f_ok
  mu[ok] <- num[ok] / (f[ok] * den[ok])
  list(nw = nw, mu = mu, bias = mu / mu[1L], mean = mu[1L] * f)
}

# Mean trajectories (Section 4) on a dense grid of n_grid ages.
.mean_curves <- function(design, phi, h_mu_u, kern, n_grid) {
  tau0 <- design$tau0
  tau <- design$tau
  t_dense <- seq(tau0, design$tau1, length.out = n_grid)
  s_dense <- (t_dense - tau0) / tau
  phi_fun <- .interp_unit(design$s_grid, phi)
  f_dense <- phi_fun(s_dense) / tau
  b <- .mu_from_sums(.ksum(s_dense, design$baseline$u, cbind(design$baseline$y, 1),
                           h_mu_u, kern), f_dense)
  p <- .mu_from_sums(.ksum(s_dense, design$pooled$u, cbind(design$pooled$y, 1),
                           h_mu_u, kern), f_dense)
  curves <- data.frame(t = t_dense, shape = f_dense,
                       mean_baseline = b$mean, bias_baseline = b$bias,
                       mu_baseline = b$mu, nw_baseline = b$nw,
                       mean_pooled = p$mean, bias_pooled = p$bias,
                       mu_pooled = p$mu, nw_pooled = p$nw)
  list(curves = curves, mu0 = c(baseline = b$mu[1L], pooled = p$mu[1L]),
       phi_fun = phi_fun)
}
