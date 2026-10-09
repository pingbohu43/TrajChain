# S3 methods: print, summary, predict, as.data.frame.

.curve_column <- function(what, mu_method) {
  what <- match.arg(what, c("shape", "mean", "bias", "mu", "nw"))
  if (what == "shape") "shape" else paste0(what, "_", mu_method)
}

.curve_labels <- c(shape = "Shape function", mean = "Unconditional mean",
                   bias = "Birth-cohort-survivor bias", mu = "Mean of latent variable",
                   nw = "Simple kernel smoothing")

.default_ages <- function(tau0, tau1, n = 6L) seq(tau0, tau1, length.out = n)

.method_label <- function(x) {
  if (identical(x$method, "chain")) {
    "ratio chaining (Section 3)"
  } else {
    sprintf("penalized ratio matching (Section 5), lambda = %s, constraint = \"%s\"",
            .format_num(x$lambda, 3), x$constraint)
  }
}

.mu_label <- function(mu_method) {
  if (identical(mu_method, "baseline")) {
    "baseline measurements (equation 5)"
  } else {
    "all observed visits (equation 6)"
  }
}

.gap_label <- function(gaps) {
  k <- length(gaps)
  if (all(abs(gaps - gaps[1L]) <= 1e-12 * gaps[1L])) {
    sprintf("%s, equally spaced (%d follow-up visit%s)", .format_num(gaps[1L]), k,
            if (k == 1L) "" else "s")
  } else {
    sprintf("%s (%d follow-up visits)", paste(.format_num(gaps), collapse = ", "), k)
  }
}

#' Print, summarize and extract a trajchain fit
#'
#' Methods for objects returned by [trajchain()]. `print()` and `summary()`
#' describe the design and report the estimates at a few ages; `predict()`
#' evaluates the estimated curves at arbitrary ages; `as.data.frame()` returns
#' the curves on the dense age grid.
#'
#' @param x,object A `"trajchain"` object.
#' @param ages Ages at which to report the estimates (default: six equally
#'   spaced ages spanning \eqn{[\tau_0, \tau_1]}).
#' @param mu_method `"baseline"` or `"pooled"`: which estimator of
#'   \eqn{\mu_{Z^0}} to report. Defaults to the one chosen in [trajchain()].
#' @param newdata Numeric vector of ages at which to evaluate the curve.
#' @param type The curve to evaluate: `"shape"` (\eqn{\hat f}), `"mean"`
#'   (unconditional mean), `"bias"` (birth-cohort-survivor bias), `"mu"`
#'   (\eqn{\hat\mu_{Z^0}}) or `"nw"` (Nadaraya--Watson curve).
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `"summary.trajchain"` whose `table` component is a data frame of
#'   estimates at `ages`. `predict()` returns a numeric vector (`NA` outside
#'   \eqn{[\tau_0, \tau_1]}). `as.data.frame()` returns `x$curves`.
#'
#' @details `predict()` recomputes the kernel sums at the requested ages when
#'   the fit contains the data (`keep_data = TRUE`), so its values agree exactly
#'   with `x$curves` at the grid ages; otherwise it interpolates `x$curves`
#'   linearly.
#'
#' @examples
#' sim <- simulate_cohort(n = 500, seed = 5)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
#' summary(fit, ages = c(0.1, 0.5, 0.9))
#' predict(fit, c(0.25, 0.75), type = "shape")
#' predict(fit, c(0.25, 0.75), type = "bias", mu_method = "pooled")
#' head(as.data.frame(fit))
#'
#' @name trajchain-methods
NULL

#' @rdname trajchain-methods
#' @export
print.trajchain <- function(x, digits = 4, ...) {
  d <- x$design
  cat("<TrajChain fit>\n")
  cat(sprintf("  Shape estimator : %s\n", .method_label(x)))
  cat(sprintf("  Mean estimator  : %s\n", .mu_label(x$mu_method)))
  cat(sprintf("  Subjects        : %d, with %d measurements (%.2f per subject)\n",
              d$n, d$n_obs, d$n_obs / d$n))
  cat(sprintf("  Visit gaps      : %s\n", .gap_label(d$visit_gaps)))
  cat(sprintf("  Age interval    : [%s, %s], m = %d grid points (step %s)\n",
              .format_num(d$tau0), .format_num(d$tau1), d$m, .format_num(d$grid_step)))
  cat(sprintf("  Bandwidths      : h = %s (ratio), h' = %s (mean)\n",
              .format_num(x$bandwidth[["ratio"]]), .format_num(x$bandwidth[["mean"]])))
  cat("\n")
  tab <- summary(x)$table
  print(tab, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @rdname trajchain-methods
#' @export
summary.trajchain <- function(object, ages = NULL, mu_method = NULL, ...) {
  mu_method <- if (is.null(mu_method)) object$mu_method else
    match.arg(mu_method, c("baseline", "pooled"))
  d <- object$design
  if (is.null(ages)) ages <- .default_ages(d$tau0, d$tau1)
  if (!is.numeric(ages) || any(!is.finite(ages))) .stopf("`ages` must be finite numbers.")
  cur <- object$curves
  get <- function(col) stats::approx(cur$t, cur[[col]], xout = ages, rule = 1)$y
  tab <- data.frame(age = ages,
                    shape = get("shape"),
                    mean = get(paste0("mean_", mu_method)),
                    bias = get(paste0("bias_", mu_method)),
                    nw = get(paste0("nw_", mu_method)))
  structure(list(table = tab, mu_method = mu_method, fit = object),
            class = "summary.trajchain")
}

#' @rdname trajchain-methods
#' @export
print.summary.trajchain <- function(x, digits = 4, ...) {
  f <- x$fit
  d <- f$design
  cat("<Summary of TrajChain fit>\n")
  cat(sprintf("  Shape estimator : %s\n", .method_label(f)))
  cat(sprintf("  Mean estimator  : %s\n", .mu_label(x$mu_method)))
  cat(sprintf("  Subjects        : %d, with %d measurements (%.2f per subject)\n",
              d$n, d$n_obs, d$n_obs / d$n))
  if (d$n_intermittent > 0) {
    cat(sprintf("                    %d subject(s) with intermittent missing visits\n",
                d$n_intermittent))
  }
  cat(sprintf("  Visit gaps      : %s\n", .gap_label(d$visit_gaps)))
  cat(sprintf("  Age interval    : [%s, %s], m = %d grid points (step %s)\n",
              .format_num(d$tau0), .format_num(d$tau1), d$m, .format_num(d$grid_step)))
  cat(sprintf("  Bandwidths      : h = %s (ratio), h' = %s (mean), kernel: %s\n",
              .format_num(f$bandwidth[["ratio"]]), .format_num(f$bandwidth[["mean"]]),
              f$kernel))
  if (!is.null(f$cv$bandwidth)) {
    cat(sprintf("                    ratio bandwidth constant %s chosen by %d-fold CV\n",
                .format_num(f$cv$bandwidth$constant), f$cv$bandwidth$settings$n_folds))
  }
  if (identical(f$method, "penalized")) {
    kept <- sum(f$pairs$kept)
    cat(sprintf("  Ratio pairs     : %d of %d kept (c* = threshold * n * h, threshold = %s)\n",
                kept, nrow(f$pairs), .format_num(f$threshold)))
    if (!is.null(f$cv$lambda)) {
      cat(sprintf("                    lambda chosen by %d-fold CV over %d values\n",
                  f$cv$lambda$settings$n_folds, length(f$cv$lambda$lambda_grid)))
    }
  }
  cat(sprintf("  mu(tau0)        : %s (baseline), %s (pooled)\n",
              .format_num(f$mu0[["baseline"]]), .format_num(f$mu0[["pooled"]])))
  cat("\nEstimates (shape = f(t); mean = unconditional mean; bias = birth-cohort-survivor bias;\n")
  cat("nw = Nadaraya-Watson smoothing of the observed markers):\n")
  print(x$table, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @rdname trajchain-methods
#' @export
predict.trajchain <- function(object, newdata, type = c("shape", "mean", "bias", "mu", "nw"),
                              mu_method = NULL, ...) {
  type <- match.arg(type)
  mu_method <- if (is.null(mu_method)) object$mu_method else
    match.arg(mu_method, c("baseline", "pooled"))
  if (missing(newdata)) return(object$curves[[.curve_column(type, mu_method)]])
  if (!is.numeric(newdata) || any(!is.finite(newdata))) {
    .stopf("`newdata` must be a numeric vector of finite ages.")
  }
  t <- as.numeric(newdata)
  f <- object$shape_fun(t)
  if (type == "shape") return(f)
  d <- object$design
  inside <- t >= d$tau0 & t <= d$tau1
  if (is.null(object$data)) {
    cur <- object$curves
    return(stats::approx(cur$t, cur[[.curve_column(type, mu_method)]], xout = t, rule = 1)$y)
  }
  st <- object$settings
  des <- .prepare_design(object$data$Y, object$data$entry_age, st$visit_gaps, st$tau0,
                         tau1 = st$tau1, grid_step = st$grid_step, warn = FALSE)
  src <- if (mu_method == "baseline") des$baseline else des$pooled
  s <- (t - d$tau0) / (d$tau1 - d$tau0)
  S <- .ksum(s, src$u, cbind(src$y, 1), st$h_mu / des$tau, st$kern)
  res <- .mu_from_sums(S, f)
  mu0 <- object$mu0[[mu_method]]
  out <- switch(type,
                nw = res$nw,
                mu = res$mu,
                bias = res$mu / mu0,
                mean = mu0 * f)
  out[!inside] <- NA_real_
  out
}

#' @rdname trajchain-methods
#' @export
as.data.frame.trajchain <- function(x, ...) x$curves

#' Summaries of bootstrap results
#'
#' Methods for objects returned by [trajchain_boot()]. `summary()` reports the
#' estimate, bootstrap standard error and pointwise confidence limits of
#' selected curves at selected ages (values between grid ages are obtained by
#' linear interpolation of the pointwise quantities; ages outside
#' \eqn{[\tau_0, \tau_1]} give `NA`). `as.data.frame()` returns
#' all curves in long format, convenient for plotting with other graphics
#' packages.
#'
#' @param x,object A `"trajchain_boot"` object.
#' @param ages Ages at which to report the results (default: six equally spaced
#'   ages spanning \eqn{[\tau_0, \tau_1]}).
#' @param which Curves to report: any of `"shape"`, `"mean"`, `"bias"`, `"mu"`
#'   and `"nw"`.
#' @param mu_method `"baseline"` or `"pooled"`; defaults to the estimator
#'   chosen when fitting.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `summary()` returns a data frame with columns `quantity`, `age`,
#'   `estimate`, `se`, `lower` and `upper`. `as.data.frame()` returns a data
#'   frame with columns `curve`, `t`, `estimate`, `se`, `lower`, `upper` and
#'   `n_ok`. `print()` returns `x` invisibly.
#'
#' @examples
#' sim <- simulate_cohort(n = 400, seed = 6)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  n_grid = 101)
#' bt <- trajchain_boot(fit, B = 30, seed = 2)
#' summary(bt, ages = c(0.25, 0.5, 0.75), which = "bias")
#' head(as.data.frame(bt))
#'
#' @name trajchain_boot-methods
NULL

#' @rdname trajchain_boot-methods
#' @export
print.trajchain_boot <- function(x, digits = 4, ...) {
  cat("<TrajChain bootstrap>\n")
  cat(sprintf("  %d subject-level bootstrap replicates (%d failed), n = %d subjects\n",
              x$B, x$n_failed, x$n))
  cat(sprintf("  %g%% pointwise percentile bands; bandwidths held fixed at h = %s, h' = %s\n",
              100 * x$conf, .format_num(x$bandwidth[["ratio"]]),
              .format_num(x$bandwidth[["mean"]])))
  if (identical(x$method, "penalized")) {
    cat(sprintf("  penalty held fixed at lambda = %s\n", .format_num(x$lambda)))
  }
  cat("\nStandard errors at selected ages:\n")
  print(summary(x), digits = digits, row.names = FALSE)
  invisible(x)
}

#' @rdname trajchain_boot-methods
#' @export
summary.trajchain_boot <- function(object, ages = NULL, which = c("shape", "mean", "bias"),
                                   mu_method = NULL, ...) {
  mu_method <- if (is.null(mu_method)) object$mu_method else
    match.arg(mu_method, c("baseline", "pooled"))
  which <- match.arg(which, c("shape", "mean", "bias", "mu", "nw"), several.ok = TRUE)
  if (is.null(ages)) ages <- .default_ages(object$design$tau0, object$design$tau1)
  if (!is.numeric(ages) || any(!is.finite(ages))) .stopf("`ages` must be finite numbers.")
  gv <- function(df, col) stats::approx(df$t, df[[col]], xout = ages, rule = 1)$y
  out <- do.call(rbind, lapply(which, function(w) {
    col <- .curve_column(w, mu_method)
    data.frame(quantity = w, age = ages,
               estimate = gv(object$est, col), se = gv(object$se, col),
               lower = gv(object$lower, col), upper = gv(object$upper, col))
  }))
  rownames(out) <- NULL
  out
}

#' @rdname trajchain_boot-methods
#' @export
as.data.frame.trajchain_boot <- function(x, ...) {
  vars <- setdiff(names(x$est), "t")
  out <- do.call(rbind, lapply(vars, function(v) {
    data.frame(curve = v, t = x$t, estimate = x$est[[v]], se = x$se[[v]],
               lower = x$lower[[v]], upper = x$upper[[v]], n_ok = x$n_ok[[v]])
  }))
  rownames(out) <- NULL
  out
}

#' @export
print.trajchain_cv <- function(x, digits = 4, ...) {
  if (identical(x$type, "bandwidth")) {
    cat("<TrajChain cross-validation: bandwidth constant>\n")
    cat(sprintf("  %d-fold CV over subjects, %d partition(s), %s training-fold scale\n",
                x$settings$n_folds, x$settings$n_repeats,
                if (x$settings$sd_scope == "fold") "fold-specific" else "full-sample"))
    cat(sprintf("  selected constant %s, i.e. h = %s * sd(A) * n^(-1/6) = %s\n\n",
                .format_num(x$constant), .format_num(x$constant), .format_num(x$bandwidth)))
    print(x$table, digits = digits, row.names = FALSE)
  } else {
    cat("<TrajChain cross-validation: penalty lambda>\n")
    cat(sprintf("  %d-fold CV over subjects, %d candidate values, ratio bandwidth h = %s\n",
                x$settings$n_folds, length(x$lambda_grid), .format_num(x$bandwidth)))
    cat(sprintf("  selected lambda = %s\n\n", .format_num(x$lambda)))
    print(data.frame(lambda = x$lambda_grid, cv = x$cv), digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.trajchain_sim <- function(x, ...) {
  s <- x$settings
  cat("<TrajChain simulated cohort>\n")
  cat(sprintf("  n = %d subjects enrolled (%.1f%% of the %d simulated individuals were alive at entry)\n",
              s$n, 100 * x$prop_enrolled, s$pop_size))
  cat(sprintf("  truncation: %s; visits: %s\n", s$truncation, .gap_label(x$visit_gaps)))
  cat(sprintf("  mean number of observed measurements per subject: %.2f\n", x$mean_visits))
  cat("  pass Y, entry_age and visit_gaps to trajchain() with tau0 = 0 and tau1 = 1\n")
  invisible(x)
}
