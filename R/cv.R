#' Cross-validated choice of the bandwidth constant
#'
#' Chooses the constant \eqn{\kappa} in the ratio bandwidth
#' \eqn{h = \kappa\,\hat\sigma_A n^{-1/6}}{h = kappa sd(A) n^(-1/6)} by
#' \eqn{F}-fold cross-validation over subjects, using the criterion proposed in
#' the supplementary material of Hu et al. For a visit gap \eqn{v} and each
#' fixed age \eqn{u}, the population minimizer over \eqn{r} of
#' \eqn{E[Y(u - v) r^2 - 2 Y(u) r \mid dN(u) = 1]} is the ratio
#' \eqn{r(u, u - v) = f(u)/f(u - v)}. The criterion is the out-of-sample version
#' of this loss,
#' \deqn{CV(\kappa) = \frac{1}{F}\sum_{j=1}^{F} \frac{1}{n_j} \sum_{i \in I_j}
#'   \sum_{v} \int_{\tau_0 + v}^{\tau_1} \left[ Y_i(u - v)\{\hat r^{(-j)}(u, u - v)\}^2
#'   - 2 Y_i(u)\, \hat r^{(-j)}(u, u - v) \right] dN_i^{(v)}(u),}{
#'   CV(kappa) = (1/F) sum_j (1/n_j) sum_{i in I_j} sum_v int [ Y_i(u - v) r_(-j)(u, u - v)^2
#'   - 2 Y_i(u) r_(-j)(u, u - v) ] dN_i^(v)(u),}
#' where \eqn{\hat r^{(-j)}} is the ratio estimator computed without fold
#' \eqn{j}, with bandwidth \eqn{\kappa\,\hat\sigma_{A} n_{-j}^{-1/6}}, the sum
#' over \eqn{v} runs over the distinct visit gaps (a single term for equally
#' spaced visits), and the integral is a sum over the observed visits of the
#' test subjects whose preceding visit lies \eqn{v} earlier. The criterion
#' involves only the ratio estimators and is therefore free of the penalty
#' \eqn{\lambda} of the penalized estimator. Test points at which the training
#' ratio is undefined (possible for very small constants) are dropped; their
#' number is reported in the `n_dropped` column of `table`. A flat criterion
#' means that the data hardly discriminate between bandwidths; it is then
#' advisable to check the sensitivity of the conclusions to \eqn{\kappa}.
#'
#' @inheritParams trajchain
#' @param candidates Candidate values of the constant \eqn{\kappa}. The default
#'   is the set used in the paper.
#' @param n_folds Number of folds \eqn{F}.
#' @param n_repeats Number of independent random partitions into folds; the
#'   criterion is averaged over partitions (seeds `seed`, `seed + 1`, ...). Use
#'   more than one partition for small samples, where a single partition adds
#'   noticeable noise.
#' @param seed Seed for the random partition, or `NULL` to use the current
#'   random number stream. The caller's random number stream is left unchanged
#'   when a seed is supplied.
#' @param sd_scope Scale \eqn{\hat\sigma_A} used in the training-fold
#'   bandwidth: the standard deviation of the training subjects' entry ages
#'   (`"fold"`, default) or of all subjects (`"full"`).
#'
#' @return An object of class `"trajchain_cv"` (with `type = "bandwidth"`), a
#'   list with components
#'   \describe{
#'     \item{`constant`}{the selected constant \eqn{\hat\kappa};}
#'     \item{`bandwidth`}{the corresponding full-data bandwidth
#'       \eqn{\hat\kappa\,\hat\sigma_A n^{-1/6}}, ready to be passed to
#'       [trajchain()];}
#'     \item{`table`}{data frame with each candidate, its criterion, its
#'       full-data bandwidth and the average number of test points dropped
#'       because the training ratio was undefined;}
#'     \item{`cv_by_repeat`, `cv_folds`}{criterion for every partition, and for
#'       every fold of the first partition;}
#'     \item{`n_test`, `settings`}{average number of test points per partition,
#'       and the settings used.}
#'   }
#'
#' @template ref
#' @seealso [trajchain()], which accepts the result as its `bandwidth` argument;
#'   [plot.trajchain_cv()] to inspect the criterion.
#'
#' @examples
#' sim <- simulate_cohort(n = 600, seed = 3)
#' cvb <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
#' cvb
#' plot(cvb)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  bandwidth = cvb)
#' fit$bandwidth
#'
#' @export
cv_bandwidth <- function(Y, entry_age, visit_gaps, tau0, tau1 = NULL, m = NULL,
                         grid_step = NULL,
                         candidates = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4),
                         n_folds = 5, n_repeats = 1, seed = 1,
                         sd_scope = c("fold", "full"), kernel = "epanechnikov") {
  sd_scope <- match.arg(sd_scope)
  kern <- kernel_function(kernel)
  design <- .prepare_design(Y, entry_age, visit_gaps, tau0, tau1 = tau1, m = m,
                            grid_step = grid_step)
  out <- .cv_bandwidth_core(design, kern = kern, candidates = candidates,
                            n_folds = n_folds, n_repeats = n_repeats, seed = seed,
                            sd_scope = sd_scope)
  out$call <- match.call()
  out
}

.cv_bandwidth_core <- function(design, kern,
                               candidates = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4),
                               n_folds = 5L, n_repeats = 1L, seed = 1L,
                               sd_scope = "fold") {
  if (!is.numeric(candidates) || !length(candidates) || any(!is.finite(candidates)) ||
      any(candidates <= 0)) {
    .stopf("`candidates` must be positive numbers.")
  }
  n_folds <- .check_count(n_folds, "n_folds", min = 2L)
  n_repeats <- .check_count(n_repeats, "n_repeats", min = 1L)
  seed <- .check_seed(seed)
  n <- design$n
  if (n < n_folds) .stopf("Cross-validation needs at least `n_folds` = %d subjects.", n_folds)
  A <- design$entry_age
  sd_full <- stats::sd(A)
  tau <- design$tau
  m <- design$m
  tol <- 1e-8
  nc <- length(candidates)

  if (!is.null(seed)) {
    old <- .get_rng_state()
    on.exit(.restore_rng_state(old), add = TRUE)
  }
  cv_rep <- matrix(NA_real_, nrow = n_repeats, ncol = nc)
  drop_rep <- matrix(0, nrow = n_repeats, ncol = nc)
  n_test <- numeric(n_repeats)
  folds_first <- NULL
  for (r in seq_len(n_repeats)) {
    if (!is.null(seed)) set.seed(seed + r - 1L)
    fold_id <- sample(rep_len(seq_len(n_folds), n))
    cv_mat <- matrix(NA_real_, nrow = n_folds, ncol = nc)
    for (j in seq_len(n_folds)) {
      test <- fold_id == j
      n_j <- sum(test)
      n_mj <- n - n_j
      scale_j <- if (sd_scope == "fold") stats::sd(A[!test]) else sd_full
      parts <- list()
      for (key in names(design$pairs)) {
        pg <- design$pairs[[key]]
        tr <- !test[pg$id]
        te <- test[pg$id] & pg$u >= pg$gap / m - tol & pg$u <= 1 + tol
        if (!any(tr) || !any(te)) next
        parts[[key]] <- list(D = outer(pg$u[te], pg$u[tr], "-"),
                             ycur_tr = pg$ycur[tr], ylag_tr = pg$ylag[tr],
                             ycur_te = pg$ycur[te], ylag_te = pg$ylag[te])
        n_test[r] <- n_test[r] + sum(te)
      }
      if (!length(parts)) next
      for (ci in seq_len(nc)) {
        h_u <- candidates[ci] * scale_j * n_mj^(-1 / 6) / tau
        tot <- 0
        bad <- 0
        any_ok <- FALSE
        for (pt in parts) {
          K <- kern(pt$D / h_u)
          num <- as.vector(K %*% pt$ycur_tr)
          den <- as.vector(K %*% pt$ylag_tr)
          r_te <- ifelse(den > 0, num / den, NA_real_)
          contrib <- pt$ylag_te * r_te^2 - 2 * pt$ycur_te * r_te
          ok <- is.finite(contrib)
          bad <- bad + sum(!ok)
          if (any(ok)) {
            tot <- tot + sum(contrib[ok])
            any_ok <- TRUE
          }
        }
        drop_rep[r, ci] <- drop_rep[r, ci] + bad
        cv_mat[j, ci] <- if (any_ok) tot / n_j else NA_real_
      }
    }
    cv_rep[r, ] <- colMeans(cv_mat, na.rm = TRUE)
    if (r == 1L) folds_first <- cv_mat
  }
  cv <- colMeans(cv_rep, na.rm = TRUE)
  if (all(!is.finite(cv))) {
    .stopf("The cross-validation criterion could not be computed for any candidate constant.")
  }
  best <- which.min(cv)
  h_unit_const <- bw_rot(A, -1 / 6)
  structure(list(
    type = "bandwidth",
    constant = candidates[best],
    bandwidth = candidates[best] * h_unit_const,
    candidates = candidates,
    cv = cv,
    table = data.frame(constant = candidates, cv = cv,
                       bandwidth = candidates * h_unit_const,
                       n_dropped = colMeans(drop_rep)),
    cv_by_repeat = cv_rep,
    cv_folds = folds_first,
    n_test = mean(n_test),
    settings = list(n_folds = n_folds, n_repeats = n_repeats, seed = seed,
                    sd_scope = sd_scope)
  ), class = "trajchain_cv")
}

#' Cross-validated choice of the penalty of the penalized estimator
#'
#' Chooses the penalty \eqn{\lambda} of the penalized shape estimator
#' (`method = "penalized"` in [trajchain()]) by \eqn{F}-fold cross-validation
#' over subjects. For each fold, the shape function \eqn{\hat f^{(-j)}} is
#' fitted on the training subjects and the ratio \eqn{\hat r_{(j)}} is
#' estimated on the test subjects; the criterion
#' \deqn{\sum_{v} \int_{\tau_0 + v}^{\tau_1} w_v(t) \left\{\hat f^{(-j)}(t) -
#'   \hat r_{(j)}(t, t - v)\, \hat f^{(-j)}(t - v)\right\}^2 dt}{
#'   sum_v int w_v(t) {f_(-j)(t) - r_(j)(t, t - v) f_(-j)(t - v)}^2 dt}
#' measures how well the training fit reproduces the ratios seen in held-out
#' subjects. Here \eqn{w_v(t) = n_v(t) I\{n_v(t) \ge c^*\}}, \eqn{n_v(t)} is the
#' number of training observations with preceding gap \eqn{v} within one
#' bandwidth of \eqn{t}, and the integral is approximated by a Riemann sum with
#' step `dt` on the unit time scale. The criterion is averaged over folds and
#' minimized over `lambda_grid`.
#'
#' @inheritParams trajchain
#' @inheritParams cv_bandwidth
#' @param bandwidth Ratio bandwidth \eqn{h} on the age scale: `"rot"`, `"cv"`,
#'   the result of [cv_bandwidth()], or a positive number. Within each training
#'   fold it is rescaled by the ratio of the training-fold to the full-sample
#'   standard deviation of the entry ages.
#' @param lambda_grid Candidate values of \eqn{\lambda}. A value whose fit fails
#'   in any training fold gets an infinite criterion. The default is the
#'   logarithmic grid of 25 values from \eqn{10^{-6}} to \eqn{10^{-2}} used in
#'   the paper.
#' @param dt Step of the Riemann sum on the unit time scale
#'   \eqn{(t - \tau_0)/(\tau_1 - \tau_0)}; the default is \eqn{1/(4m)}.
#' @param seed Seed for the random partition into folds, or `NULL` to use the
#'   current random number stream.
#'
#' @return An object of class `"trajchain_cv"` (with `type = "lambda"`), a list
#'   with components `lambda` (the selected value), `lambda_grid`, `cv` (the
#'   criterion averaged over folds; `Inf` where the fit failed), `scores` (the
#'   criterion for every value and fold), `bandwidth` and `settings`.
#'
#' @template ref
#' @seealso [trajchain()], which accepts the result as its `lambda` argument.
#'
#' @examples
#' \donttest{
#' sim <- simulate_cohort(n = 1000, visit_gaps = "unequal", seed = 11)
#' cvl <- cv_lambda(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
#' cvl
#' plot(cvl)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  lambda = cvl, mu_method = "pooled")
#' }
#'
#' @export
cv_lambda <- function(Y, entry_age, visit_gaps, tau0, tau1 = NULL, m = NULL,
                      grid_step = NULL, bandwidth = "rot",
                      lambda_grid = exp(seq(log(1e-6), log(1e-2), length.out = 25)),
                      threshold = 0.6, n_folds = 5, seed = 1, dt = NULL,
                      constraint = c("anchor", "integral"),
                      kernel = "epanechnikov") {
  constraint <- match.arg(constraint)
  kern <- kernel_function(kernel)
  design <- .prepare_design(Y, entry_age, visit_gaps, tau0, tau1 = tau1, m = m,
                            grid_step = grid_step)
  if (design$m < 3L) .stopf("The penalized estimator needs at least three grid points (m >= 3).")
  bw <- .resolve_ratio_bandwidth(bandwidth, design, kern)
  out <- .cv_lambda_core(design, h = bw$h, lambda_grid = lambda_grid,
                         threshold = threshold, n_folds = n_folds, seed = seed,
                         dt = dt, constraint = constraint, kern = kern)
  out$call <- match.call()
  out
}

.cv_lambda_core <- function(design, h, kern,
                            lambda_grid = exp(seq(log(1e-6), log(1e-2), length.out = 25)),
                            threshold = 0.6, n_folds = 5L, seed = 1L, dt = NULL,
                            constraint = "anchor") {
  if (!is.numeric(lambda_grid) || !length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0)) {
    .stopf("`lambda_grid` must contain non-negative numbers.")
  }
  threshold <- .check_number(threshold, "threshold", nonnegative = TRUE)
  n_folds <- .check_count(n_folds, "n_folds", min = 2L)
  seed <- .check_seed(seed)
  m <- design$m
  if (is.null(dt)) dt <- 1 / (4 * m)
  dt <- .check_number(dt, "dt", positive = TRUE)
  n <- design$n
  if (n < n_folds) .stopf("Cross-validation needs at least `n_folds` = %d subjects.", n_folds)
  A <- design$entry_age
  sd_full <- stats::sd(A)
  tau <- design$tau

  if (!is.null(seed)) {
    old <- .get_rng_state()
    on.exit(.restore_rng_state(old), add = TRUE)
    set.seed(seed)
  }
  fold_id <- sample(rep(seq_len(n_folds), length.out = n))
  scores <- matrix(NA_real_, nrow = length(lambda_grid), ncol = n_folds)
  for (fk in seq_len(n_folds)) {
    test_rows <- which(fold_id == fk)
    train_rows <- which(fold_id != fk)
    train <- .subset_design(design, train_rows)
    test <- .subset_design(design, test_rows)
    h_u <- (h * stats::sd(A[train_rows]) / sd_full) / tau
    sys <- tryCatch(.penalized_system(train, h_u, threshold, kern),
                    error = function(e) NULL)
    if (is.null(sys)) {
      scores[, fk] <- Inf
      next
    }
    cstar <- threshold * train$n * h_u
    terms <- list()
    for (key in names(train$pairs)) {
      v_ch <- train$pairs[[key]]$gap * (1 / m)
      if (v_ch >= 1) next   # gap longer than the age interval: no ratio to compare
      tg <- seq(v_ch, 1, by = dt)
      if (length(tg) < 2L) next
      pte <- test$pairs[[key]]
      S <- .ksum(tg, pte$u, cbind(pte$ycur, pte$ylag), h_u, kern)
      r_te <- ifelse(is.finite(S[, 2L]) & abs(S[, 2L]) > 0, S[, 1L] / S[, 2L], NA_real_)
      n_t <- .count_window(tg, train$pairs[[key]]$u, h_u)
      terms[[key]] <- list(tg = tg, lag = tg - v_ch, r = r_te, w = n_t * (n_t >= cstar))
    }
    for (li in seq_along(lambda_grid)) {
      phi <- tryCatch(.penalized_solve(sys, m, lambda_grid[li], constraint)$phi,
                      error = function(e) NULL)
      if (is.null(phi)) {
        scores[li, fk] <- Inf
        next
      }
      fun <- .interp_unit(train$s_grid, phi)
      total <- 0
      for (tm in terms) {
        resid <- fun(tm$tg) - tm$r * fun(tm$lag)
        val <- tm$w * resid^2
        ok <- is.finite(val)
        total <- total + sum(val[ok]) * dt
      }
      scores[li, fk] <- total
    }
  }
  cv <- rowMeans(scores, na.rm = TRUE)
  cv[!is.finite(cv)] <- Inf
  if (all(!is.finite(cv))) {
    .stopf(paste0("For every candidate lambda, the penalized estimator could not be fitted ",
                  "in at least one training fold; increase `bandwidth` or decrease `threshold`."))
  }
  best <- which.min(cv)
  structure(list(
    type = "lambda",
    lambda = lambda_grid[best],
    lambda_grid = lambda_grid,
    cv = cv,
    scores = scores,
    bandwidth = h,
    settings = list(n_folds = n_folds, seed = seed, dt = dt, threshold = threshold,
                    constraint = constraint)
  ), class = "trajchain_cv")
}
