#' Subject-level bootstrap standard errors and confidence bands
#'
#' Resamples subjects with replacement and refits the estimators of a
#' [trajchain()] fit on each bootstrap sample, holding the bandwidths, the
#' penalty and all other settings fixed at their full-sample values. Pointwise
#' bootstrap standard errors and percentile confidence bands are computed for
#' every curve in `fit$curves` (shape function, unconditional mean,
#' birth-cohort-survivor bias, \eqn{\hat\mu_{Z^0}} and the Nadaraya--Watson
#' curve, for both mean estimators).
#'
#' @param fit A `"trajchain"` object fitted with `keep_data = TRUE`.
#' @param B Number of bootstrap replicates.
#' @param seed Seed for the resampling, or `NULL` to use the current random
#'   number stream. The bootstrap samples depend only on `seed` (not on
#'   `cores`), and the caller's random number stream is left unchanged when a
#'   seed is supplied.
#' @param conf Confidence level of the pointwise bands.
#' @param cores Number of cores. Values larger than one use
#'   [parallel::mclapply()] and are ignored on Windows.
#' @param keep_curves Keep the \eqn{B \times} `n_grid` matrices of bootstrap
#'   curves (needed for custom summaries).
#' @param progress Print a progress message every 100 replicates (only used
#'   with `cores = 1`).
#'
#' @details A replicate can fail, for instance when a resampled data set leaves
#'   a grid age without observed visit pairs within one bandwidth, which makes
#'   the ratio estimator undefined. Failed replicates are dropped; their number
#'   is reported in `n_failed`, and `n_ok` gives, for every curve and age, the
#'   number of replicates with a finite value. By construction the
#'   birth-cohort-survivor bias equals one at \eqn{\tau_0} in every replicate,
#'   so its band has zero width there; this reflects the normalization, not
#'   precise estimation of \eqn{\mu_{Z^0}(\tau_0)}, whose uncertainty is
#'   carried by the bias at other ages and by the unconditional mean.
#'
#' @return An object of class `"trajchain_boot"`, a list with components
#'   \describe{
#'     \item{`est`, `se`, `lower`, `upper`, `boot_mean`, `n_ok`}{data frames
#'       with column `t` and one column per curve of `fit$curves`: the
#'       full-sample estimates, the bootstrap standard errors, the pointwise
#'       percentile confidence limits, the bootstrap means and the number of
#'       finite bootstrap values;}
#'     \item{`n_incomplete`}{for every curve, the number of replicates with at
#'       least one non-finite value;}
#'     \item{`n_failed`}{number of replicates in which the fit failed;}
#'     \item{`curves`}{if `keep_curves = TRUE`, a list of \eqn{B \times}
#'       `n_grid` matrices of bootstrap curves;}
#'     \item{`B`, `conf`, `seed`, `n`, `bandwidth`, `lambda`, `mu_method`,
#'       `design`}{settings and a summary of the original fit.}
#'   }
#'
#' @template ref
#' @seealso [summary.trajchain_boot()] for standard errors at selected ages,
#'   [plot.trajchain_boot()] and [plot_trajectories()] for confidence bands.
#'
#' @examples
#' sim <- simulate_cohort(n = 500, seed = 4)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  n_grid = 201)
#' bt <- trajchain_boot(fit, B = 50, seed = 1)
#' bt
#' summary(bt, ages = c(0.2, 0.5, 0.8))
#' plot(bt, which = c("shape", "bias"))
#'
#' @export
trajchain_boot <- function(fit, B = 1000, seed = NULL, conf = 0.95, cores = 1,
                           keep_curves = TRUE, progress = FALSE) {
  if (!inherits(fit, "trajchain")) .stopf("`fit` must be an object returned by trajchain().")
  if (is.null(fit$data)) {
    .stopf("`fit` does not contain the data; refit with trajchain(..., keep_data = TRUE).")
  }
  B <- .check_count(B, "B", min = 2L)
  conf <- .check_number(conf, "conf")
  if (conf <= 0 || conf >= 1) .stopf("`conf` must lie strictly between 0 and 1.")
  cores <- .check_count(cores, "cores", min = 1L)
  .check_flag(keep_curves, "keep_curves")
  .check_flag(progress, "progress")
  seed <- .check_seed(seed)

  st <- fit$settings
  Y <- fit$data$Y
  A <- fit$data$entry_age
  n <- nrow(Y)

  if (!is.null(seed)) {
    old <- .get_rng_state()
    on.exit(.restore_rng_state(old), add = TRUE)
    set.seed(seed)
  }
  idx_list <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))

  one_fit <- function(idx) {
    tryCatch({
      d <- .prepare_design(Y[idx, , drop = FALSE], A[idx], st$visit_gaps, st$tau0,
                           tau1 = st$tau1, grid_step = st$grid_step, warn = FALSE)
      .fit_core(d, st)$curves
    }, error = function(e) NULL)
  }
  use_parallel <- cores > 1L && .Platform$OS.type != "windows"
  if (use_parallel) {
    fits <- parallel::mclapply(idx_list, one_fit, mc.cores = cores)
  } else {
    fits <- vector("list", B)
    for (b in seq_len(B)) {
      fits[[b]] <- one_fit(idx_list[[b]])
      if (progress && b %% 100L == 0L) message(sprintf("bootstrap replicate %d / %d", b, B))
    }
  }

  est <- fit$curves
  vars <- setdiff(names(est), "t")
  nd <- nrow(est)
  curves <- lapply(stats::setNames(vars, vars), function(v) matrix(NA_real_, B, nd))
  failed <- 0L
  for (b in seq_len(B)) {
    fb <- fits[[b]]
    if (is.null(fb) || !is.data.frame(fb)) {
      failed <- failed + 1L
      next
    }
    for (v in vars) curves[[v]][b, ] <- fb[[v]]
  }

  al <- (1 - conf) / 2
  as_df <- function(fun) {
    out <- data.frame(t = est$t)
    for (v in vars) out[[v]] <- fun(curves[[v]])
    out
  }
  structure(list(
    t = est$t,
    est = est,
    se = as_df(function(M) apply(M, 2L, .safe_sd)),
    lower = as_df(function(M) apply(M, 2L, .safe_quantile, p = al)),
    upper = as_df(function(M) apply(M, 2L, .safe_quantile, p = 1 - al)),
    boot_mean = as_df(function(M) colMeans(M, na.rm = TRUE)),
    n_ok = as_df(function(M) colSums(is.finite(M))),
    n_incomplete = vapply(curves, function(M) sum(!apply(is.finite(M), 1L, all)), integer(1)),
    n_failed = failed,
    curves = if (keep_curves) curves else NULL,
    B = B, conf = conf, seed = seed, n = n,
    bandwidth = fit$bandwidth, lambda = fit$lambda, method = fit$method,
    mu_method = fit$mu_method, design = fit$design
  ), class = "trajchain_boot")
}
