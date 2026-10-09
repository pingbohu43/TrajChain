sim3 <- simulate_cohort(n = 400, visit_gaps = "unequal", seed = 501)

test_that("ill-posed penalized problems give informative errors", {
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 0),
               "singular")
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 0,
                         constraint = "integral"), "could not be solved")
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 1e-4,
                         threshold = 1e6), "decrease `threshold`")
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = -1),
               "non-negative")
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, m = 2, lambda = 1),
               "at least three grid points")
})

test_that("summaries describe penalized fits and cross-validation", {
  fit <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, bandwidth = "cv",
                   lambda = "cv", n_grid = 101)
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("chosen by 5-fold CV", out)))
  expect_true(any(grepl("lambda chosen by 5-fold CV over 25 values", out)))
  expect_true(any(grepl("Ratio pairs", out)))
  expect_output(print(fit), "penalized ratio matching")
  bt <- trajchain_boot(fit, B = 5, seed = 1)
  expect_output(print(bt), "penalty held fixed")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fit$cv$lambda))
  expect_error(summary(fit, ages = NA), "finite")
  expect_error(summary(bt, ages = "a"), "finite")
})

test_that("summary reports intermittent missingness", {
  sim <- simulate_cohort(n = 300, seed = 502)
  Y <- sim$Y
  rows <- which(rowSums(!is.na(Y)) >= 4)[1:3]
  Y[rows, 2] <- NA
  fit <- suppressWarnings(trajchain(Y, sim$entry_age, sim$visit_gaps, 0, 1))
  expect_output(print(summary(fit)), "intermittent missing visits")
})

test_that("argument checks catch invalid input", {
  sim <- simulate_cohort(n = 200, seed = 503)
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, keep_data = "yes"),
               "TRUE or FALSE")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_grid = 1.5),
               "whole number")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = -1),
               "positive")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth_mu = NA),
               "single finite number")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, kernel = 1),
               "kernel name or a function")
  expect_error(bw_rot("a"), "numeric")
  expect_error(bw_rot(1:10, constant = 0), "positive")
  expect_error(ratio_estimate(sim$Y, sim$entry_age, sim$visit_gaps, t = NA), "finite ages")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1,
                         bandwidth = cv_lambda(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1,
                                               lambda_grid = c(1e-5, 1e-4))),
               "cv_bandwidth")
})

test_that("simulation helpers validate their parameters", {
  expect_error(simulate_cohort(n = 50, frailty_par = c(mean = -1, var = 0.2)), "positive")
  expect_error(sim_truth(NA), "finite numbers")
  tr <- sim_truth(c(-0.5, 0.5), truncation = "informative")
  expect_true(is.na(tr$mu[1]))
  expect_equal(sim_truth(c(0, 0.5), frailty_par = c(mean = 1, var = 0))$mu, c(1, 1))
  s <- simulate_cohort(n = 50, frailty_par = c(mean = 2, var = 0), noise_var = 0,
                       seed = 1)
  expect_equal(s$latent$Z, rep(2, 50))
})

test_that("observed (unrounded) visit gaps are caught before building a huge grid", {
  expect_error(trajchain(sim3$Y, sim3$entry_age, c(0.1, 0.1, 0.121, 0.1, 0.12), 0, 1),
               "much smaller than the gaps")
  fit <- trajchain(sim3$Y, sim3$entry_age, c(0.1, 0.1, 0.12, 0.1, 0.12), 0, 1, lambda = 1e-4)
  expect_equal(fit$design$m, 50L)
})

test_that("quantities anchored at tau0 come with warnings when unreliable", {
  d <- subset(simcohort, group == "A")
  Y <- visit_matrix(d, "id", "visit", "marker")
  A <- d$entry_age[match(rownames(Y), d$id)]
  expect_warning(trajchain(Y, A, 16 / 12, tau0 = 22, m = 22), "may be unreliable")
  expect_no_warning(trajchain(Y, A, 16 / 12, tau0 = 30, m = 22))

  beta <- 25
  steep <- function(t) beta * exp(beta * t) / (exp(beta) - 1)
  s25 <- simulate_cohort(n = 800, shape = steep, noise_var = 0, seed = 1)
  w <- character(0)
  fit <- withCallingHandlers(
    trajchain(s25$Y, s25$entry_age, s25$visit_gaps, 0, 1),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not positive", w)))
  expect_true(any(grepl("could not be estimated", w)))
  expect_true(all(is.na(fit$curves$bias_baseline)))
})
