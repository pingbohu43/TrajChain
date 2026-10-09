# Graphics (base R).

# Categorical colors in a fixed order (color-vision-deficiency checked).
.trajchain_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                        "#008300", "#4a3aa7", "#e34948")

#' Plot estimated trajectories for one or several groups
#'
#' Draws side-by-side panels of the estimated shape function, unconditional
#' mean, Nadaraya--Watson curve and birth-cohort-survivor bias (the layout of
#' Fig. 5 of Hu et al.), for one fit or for several groups (for example, women
#' and men analyzed separately). Bootstrap results add pointwise confidence
#' bands.
#'
#' @param x A `"trajchain"` or `"trajchain_boot"` object, or a named list of
#'   such objects (one per group; the names are used in the legend).
#' @param which Panels to draw, in order: any of `"shape"` (shape function),
#'   `"mean"` (unconditional mean), `"nw"` (simple kernel smoothing of the
#'   observed markers), `"bias"` (birth-cohort-survivor bias) and `"mu"`
#'   (\eqn{\hat\mu_{Z^0}}).
#' @param mu_method `"baseline"` or `"pooled"`; defaults to the estimator
#'   chosen when fitting the first object.
#' @param band Draw confidence bands for bootstrap objects.
#' @param ci Type of band: pointwise `"percentile"` intervals (default) or
#'   `"normal"` intervals, estimate \eqn{\pm z_{1-\alpha/2}} times the bootstrap
#'   standard error.
#' @param truth Optional data frame of true curves, as returned by
#'   [sim_truth()], drawn as a dashed black line (useful for simulations; its
#'   `t` must be on the age scale of the fits).
#' @param col,lty,lwd Colors, line types and line width, one per group.
#' @param xlab Label of the horizontal axis.
#' @param ylab Labels of the vertical axes, one per panel (defaults to the
#'   names used in the paper).
#' @param legend Draw a legend above the panels. Defaults to `TRUE` when there
#'   is more than one group or a `truth` curve.
#' @param include_zero Include zero in the vertical range of every panel.
#' @param mfrow Panel layout `c(rows, columns)`; defaults to a single row.
#' @param ... Further graphical arguments passed to [graphics::plot()].
#'
#' @return `NULL`, invisibly. Called for its side effect.
#'
#' @seealso [plot.trajchain()], [plot.trajchain_boot()].
#'
#' @examples
#' simA <- simulate_cohort(n = 500, shape = "monotone", seed = 1)
#' simB <- simulate_cohort(n = 500, shape = "ushape", seed = 2)
#' fitA <- trajchain(simA$Y, simA$entry_age, simA$visit_gaps, tau0 = 0, tau1 = 1)
#' fitB <- trajchain(simB$Y, simB$entry_age, simB$visit_gaps, tau0 = 0, tau1 = 1)
#' plot_trajectories(list(Monotone = fitA, "U-shaped" = fitB),
#'                   which = c("shape", "mean"))
#'
#' @export
plot_trajectories <- function(x, which = c("shape", "mean", "nw", "bias"),
                              mu_method = NULL, band = TRUE,
                              ci = c("percentile", "normal"), truth = NULL,
                              col = NULL, lty = NULL, lwd = 2, xlab = "Age",
                              ylab = NULL, legend = NULL, include_zero = TRUE,
                              mfrow = NULL, ...) {
  objs <- if (inherits(x, c("trajchain", "trajchain_boot"))) list(x) else x
  if (!is.list(objs) || !length(objs) ||
      !all(vapply(objs, inherits, logical(1), what = c("trajchain", "trajchain_boot")))) {
    .stopf("`x` must be a trajchain or trajchain_boot object, or a list of them.")
  }
  G <- length(objs)
  nms <- names(objs)
  if (is.null(nms)) nms <- if (G == 1L) "Estimate" else paste("Group", seq_len(G))
  nms[!nzchar(nms)] <- paste("Group", which(!nzchar(nms)))
  which <- match.arg(which, c("shape", "mean", "nw", "bias", "mu"), several.ok = TRUE)
  ci <- match.arg(ci)
  mu_method <- if (is.null(mu_method)) objs[[1L]]$mu_method else
    match.arg(mu_method, c("baseline", "pooled"))
  if (is.null(col)) col <- .trajchain_palette[(seq_len(G) - 1L) %% length(.trajchain_palette) + 1L]
  if (is.null(lty)) lty <- seq_len(G)
  col <- rep_len(col, G)
  lty <- rep_len(lty, G)
  lwd <- rep_len(lwd, G)
  if (is.null(ylab)) ylab <- unname(.curve_labels[which])
  ylab <- rep_len(ylab, length(which))
  if (!is.null(truth) && !(is.data.frame(truth) && "t" %in% names(truth))) {
    .stopf("`truth` must be a data frame with a column `t`, as returned by sim_truth().")
  }
  if (is.null(legend)) legend <- G > 1L || !is.null(truth)
  if (is.null(mfrow)) mfrow <- c(1L, length(which))

  op <- graphics::par(mfrow = mfrow, mar = c(4, 4.5, 1, 1),
                      oma = c(0, 0, if (legend) 2.4 else 0, 0))
  on.exit(graphics::par(op), add = TRUE)

  for (k in seq_along(which)) {
    w <- which[k]
    colname <- .curve_column(w, mu_method)
    pieces <- lapply(objs, function(o) {
      est <- if (inherits(o, "trajchain")) o$curves else o$est
      out <- list(t = est$t, y = est[[colname]], lo = NULL, hi = NULL)
      if (band && inherits(o, "trajchain_boot")) {
        if (ci == "percentile") {
          out$lo <- o$lower[[colname]]
          out$hi <- o$upper[[colname]]
        } else {
          z <- stats::qnorm(1 - (1 - o$conf) / 2)
          out$lo <- est[[colname]] - z * o$se[[colname]]
          out$hi <- est[[colname]] + z * o$se[[colname]]
        }
      }
      out
    })
    tr <- if (!is.null(truth) && w %in% names(truth)) truth[[w]] else NULL
    yr <- c(if (include_zero) 0, unlist(lapply(pieces, function(p) c(p$y, p$lo, p$hi))), tr)
    if (!any(is.finite(yr))) {
      graphics::plot.new()
      graphics::title(ylab = ylab[k], xlab = xlab)
      next
    }
    xr <- range(unlist(lapply(pieces, `[[`, "t")), finite = TRUE)
    graphics::plot(NA, xlim = xr, ylim = range(yr, finite = TRUE), xlab = xlab,
                   ylab = ylab[k], ...)
    graphics::grid(col = "grey90", lty = 1)
    for (g in seq_len(G)) {
      p <- pieces[[g]]
      if (!is.null(p$lo)) {
        ok <- is.finite(p$lo) & is.finite(p$hi)
        if (any(ok)) {
          graphics::polygon(c(p$t[ok], rev(p$t[ok])), c(p$lo[ok], rev(p$hi[ok])),
                            col = grDevices::adjustcolor(col[g], alpha.f = 0.18),
                            border = NA)
        }
      }
    }
    if (w == "bias") graphics::abline(h = 1, col = "grey40", lty = 3)
    if (!is.null(tr)) graphics::lines(truth$t, tr, col = "black", lty = 2, lwd = 1.5)
    for (g in seq_len(G)) {
      graphics::lines(pieces[[g]]$t, pieces[[g]]$y, col = col[g], lty = lty[g], lwd = lwd[g])
    }
  }

  if (legend) {
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    lab <- nms
    lcol <- col
    llty <- lty
    llwd <- lwd
    if (!is.null(truth)) {
      lab <- c(lab, "True")
      lcol <- c(lcol, "black")
      llty <- c(llty, 2)
      llwd <- c(llwd, 1.5)
    }
    graphics::legend("top", legend = lab, col = lcol, lty = llty, lwd = llwd, bty = "n",
                     horiz = TRUE)
  }
  invisible(NULL)
}

#' Plot a trajchain fit or its bootstrap confidence bands
#'
#' `plot()` methods for [trajchain()] and [trajchain_boot()] results. They call
#' [plot_trajectories()] for a single group; bootstrap results add pointwise
#' confidence bands.
#'
#' @param x A `"trajchain"` or `"trajchain_boot"` object.
#' @inheritParams plot_trajectories
#' @param ... Further arguments passed to [plot_trajectories()].
#'
#' @return `NULL`, invisibly.
#' @examples
#' sim <- simulate_cohort(n = 500, seed = 8)
#' fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
#'                  mu_method = "pooled")
#' plot(fit)
#' plot(fit, which = c("shape", "bias"), truth = sim_truth(fit$curves$t))
#' @name plot-trajchain
NULL

#' @rdname plot-trajchain
#' @export
plot.trajchain <- function(x, which = c("shape", "mean", "nw", "bias"), mu_method = NULL,
                           ...) {
  plot_trajectories(x, which = which, mu_method = mu_method, ...)
}

#' @rdname plot-trajchain
#' @export
plot.trajchain_boot <- function(x, which = c("shape", "mean", "nw", "bias"), mu_method = NULL,
                                ci = c("percentile", "normal"), ...) {
  plot_trajectories(x, which = which, mu_method = mu_method, band = TRUE,
                    ci = match.arg(ci), ...)
}

#' Plot a cross-validation criterion
#'
#' Plots the cross-validation criterion returned by [cv_bandwidth()] or
#' [cv_lambda()] against the candidate values (logarithmic axis), together with
#' the criterion for individual partitions or folds (gray) and the selected
#' value (dashed line). A U-shaped curve with an interior minimum indicates
#' that the candidate set covers the optimum.
#'
#' @param x A `"trajchain_cv"` object.
#' @param n_curves Maximum number of individual partitions or folds to draw.
#' @param ... Further graphical arguments passed to [graphics::plot()].
#' @return `NULL`, invisibly.
#' @examples
#' sim <- simulate_cohort(n = 500, seed = 9)
#' cvb <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
#' plot(cvb)
#' @export
plot.trajchain_cv <- function(x, n_curves = 30, ...) {
  main_col <- .trajchain_palette[1L]
  if (identical(x$type, "bandwidth")) {
    xx <- x$candidates
    yy <- x$cv
    M <- if (x$settings$n_repeats > 1L) x$cv_by_repeat else x$cv_folds
    lab_sp <- if (x$settings$n_repeats > 1L) "individual partitions" else "individual folds"
    xlab <- expression("constant " * kappa ~ "in" ~ h == kappa %.% hat(sigma)[A] %.% n^{-1/6})
    sel <- x$constant
  } else {
    xx <- x$lambda_grid
    yy <- x$cv
    M <- t(x$scores)
    lab_sp <- "individual folds"
    xlab <- expression("penalty " * lambda)
    sel <- x$lambda
  }
  pos <- xx > 0
  if (!all(pos)) .warnf("Non-positive candidate values are not shown on the logarithmic axis.")
  M <- M[seq_len(min(nrow(M), n_curves)), , drop = FALSE]
  yr <- range(c(yy[pos], M[, pos]), finite = TRUE)
  few <- identical(x$type, "bandwidth") && sum(pos) <= 12L
  graphics::plot(xx[pos], yy[pos], type = "n", log = "x", ylim = yr, xlab = xlab,
                 ylab = "CV criterion", xaxt = if (few) "n" else "s", ...)
  if (few) graphics::axis(1, at = xx[pos], labels = .format_num(xx[pos], 3))
  graphics::grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  for (s in seq_len(nrow(M))) graphics::lines(xx[pos], M[s, pos], col = "grey80")
  graphics::lines(xx[pos], yy[pos], type = "b", pch = 19, col = main_col, lwd = 2)
  graphics::abline(v = sel, lty = 2, col = main_col)
  graphics::legend("topright", cex = 0.85, bg = "white", box.col = NA, inset = 0.01,
                   legend = c(lab_sp, "CV criterion", sprintf("selected: %s", .format_num(sel))),
                   col = c("grey70", main_col, main_col), lty = c(1, 1, 2),
                   pch = c(NA, 19, NA), lwd = c(1, 2, 1))
  invisible(NULL)
}
