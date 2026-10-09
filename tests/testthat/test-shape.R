beta <- 1.2
f_exp <- function(t) beta * exp(beta * t) / (exp(beta) - 1)

noise_free <- function(visit_gaps, seed) {
  simulate_cohort(n = 800, shape = f_exp, visit_gaps = visit_gaps,
                  noise_var = 0, seed = seed)
}

test_that("ratio chaining recovers an exponential shape exactly without noise", {
  # f(t) / f(t - v) = exp(beta v) for every t, so every kernel ratio is exact.
  sim <- noise_free("equal", seed = 1)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
  u <- fit$grid$t
  target <- exp(beta * u) / mean(exp(beta * u))
  expect_equal(fit$grid$shape, target, tolerance = 1e-12)
  expect_equal(fit$grid$ratio[-1], rep(exp(beta * sim$visit_gaps[1]), length(u) - 1),
               tolerance = 1e-12)
})

test_that("the penalized estimator recovers an exponential shape without noise", {
  # every ratio constraint is satisfied exactly by exp(beta t); with a negligible
  # penalty the solution differs from it only by the tiny smoothing bias
  sim <- noise_free("unequal", seed = 2)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
                   lambda = 1e-10)
  u <- fit$grid$t
  target <- exp(beta * u) / mean(exp(beta * u))
  expect_equal(fit$grid$shape, target, tolerance = 1e-4)
  fit_qp <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
                      lambda = 1e-10, constraint = "integral")
  expect_equal(fit_qp$grid$shape, target, tolerance = 1e-4)
})

test_that("the shape estimate integrates to one (midpoint rule)", {
  sim <- simulate_cohort(n = 500, seed = 3)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  expect_equal(sum(fit$grid$shape) / fit$design$m, 1)
  sim3 <- simulate_cohort(n = 500, visit_gaps = "unequal", seed = 3)
  for (cons in c("anchor", "integral")) {
    fit3 <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 1e-4,
                      constraint = cons)
    expect_equal(sum(fit3$grid$shape) / fit3$design$m, 1, label = cons)
  }
  fit_qp <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 1e-4,
                      constraint = "integral")
  expect_true(all(fit_qp$grid$shape >= 0))
})

test_that("results are invariant to the age scale", {
  sim <- simulate_cohort(n = 500, seed = 4)
  fit01 <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  a <- 40
  b <- 30
  fit_age <- trajchain(sim$Y, b + a * sim$entry_age, a * sim$visit_gaps, tau0 = b,
                       tau1 = b + a)
  expect_equal(fit_age$curves$t, b + a * fit01$curves$t)
  expect_equal(fit_age$curves$shape, fit01$curves$shape / a, tolerance = 1e-10)
  for (v in c("mean_baseline", "bias_baseline", "nw_baseline", "mean_pooled",
              "bias_pooled", "nw_pooled")) {
    expect_equal(fit_age$curves[[v]], fit01$curves[[v]], tolerance = 1e-10, label = v)
  }
  expect_equal(fit_age$bandwidth, a * fit01$bandwidth)
})

test_that("estimates are close to the truth in large samples", {
  sim <- simulate_cohort(n = 2000, seed = 5)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, mu_method = "pooled")
  truth <- sim_truth(fit$curves$t)
  inner <- fit$curves$t >= 0.1 & fit$curves$t <= 0.9
  expect_lt(max(abs(fit$curves$shape - truth$shape)[inner]), 0.15)
  expect_lt(mean(abs(fit$curves$shape - truth$shape)), 0.05)
  expect_lt(mean(abs(fit$curves$bias_pooled - truth$bias)[inner]), 0.06)
  expect_lt(mean(abs(fit$curves$mean_pooled - truth$mean)), 0.06)

  sim3 <- simulate_cohort(n = 2000, shape = "ushape", visit_gaps = "unequal", seed = 6)
  fit3 <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 3e-5)
  truth3 <- sim_truth(fit3$curves$t, shape = "ushape")
  expect_lt(mean(abs(fit3$curves$shape - truth3$shape)), 0.05)
})

test_that("an undefined ratio gives an informative error", {
  sim <- simulate_cohort(n = 300, seed = 7)
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = 1e-4),
               "ratio estimate is undefined")
})

test_that("ratio_estimate() agrees with the ratios chained by trajchain()", {
  sim <- simulate_cohort(n = 400, seed = 8)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  r <- ratio_estimate(sim$Y, sim$entry_age, sim$visit_gaps, t = fit$grid$t[-1],
                      bandwidth = fit$bandwidth[["ratio"]])
  expect_equal(r$ratio, fit$grid$ratio[-1])
  expect_true(all(r$n_pairs > 0))
  sim3 <- simulate_cohort(n = 400, visit_gaps = "unequal", seed = 8)
  expect_error(ratio_estimate(sim3$Y, sim3$entry_age, sim3$visit_gaps, t = 0.5), "unequally")
  r6 <- ratio_estimate(sim3$Y, sim3$entry_age, sim3$visit_gaps, t = 0.5, gap = 0.12)
  expect_true(is.finite(r6$ratio))
  expect_error(ratio_estimate(sim3$Y, sim3$entry_age, sim3$visit_gaps, t = 0.5, gap = 0.2),
               "does not match")
})

test_that("bandwidth rules follow sd(A) * n^alpha", {
  set.seed(1)
  A <- runif(150)
  expect_equal(bw_rot(A), sd(A) * 150^(-1 / 6))
  expect_equal(bw_rot(A, -1 / 5, 2.34), 2.34 * sd(A) * 150^(-1 / 5))
  expect_error(bw_rot(1), "two finite")
  expect_error(bw_rot(rep(2, 5)), "zero spread")
  sim <- simulate_cohort(n = 300, seed = 9)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = 0.2,
                   bandwidth_mu = 0.3)
  expect_equal(unname(fit$bandwidth), c(0.2, 0.3))
})

test_that("kernels other than Epanechnikov can be used", {
  sim <- simulate_cohort(n = 400, seed = 10)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, kernel = "biweight")
  expect_identical(fit$kernel, "biweight")
  expect_true(all(is.finite(fit$grid$shape)))
})

test_that("the anchored penalized estimator with lambda = 0 reduces to ratio chaining", {
  # Section 5 of the paper: ratio chaining is the special case q = 1,
  # lambda = 0 and c* = 0 of the penalized estimator
  sim <- simulate_cohort(n = 600, seed = 11)
  f_chain <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  f_pen <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, method = "penalized",
                     lambda = 0, threshold = 0)
  expect_equal(f_pen$grid$shape, f_chain$grid$shape, tolerance = 1e-10)
  expect_equal(f_pen$curves, f_chain$curves, tolerance = 1e-10)
})

test_that("visit gaps longer than the age interval are skipped", {
  sim <- simulate_cohort(n = 800, visit_gaps = c(0.1, 0.1, 0.5), seed = 3)
  fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 0.4)
  expect_identical(fit$method, "penalized")
  expect_equal(unique(fit$pairs$gap), 0.1)
  expect_true(fit$lambda %in% fit$cv$lambda$lambda_grid)
})
