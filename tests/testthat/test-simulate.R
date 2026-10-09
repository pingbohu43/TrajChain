test_that("the shape functions are the densities of the paper", {
  expect_equal(integrate(shape_monotone, 0, 1)$value, 1, tolerance = 1e-10)
  expect_equal(integrate(shape_ushape, 0, 1)$value, 1, tolerance = 1e-10)
  expect_equal(shape_ushape(0.3), 0.8)
  expect_equal(shape_monotone(0.5), 8 / 22 + 7 / 11)
})

test_that("simulated cohorts have the expected structure", {
  sim <- simulate_cohort(n = 500, seed = 1)
  expect_s3_class(sim, "trajchain_sim")
  expect_equal(dim(sim$Y), c(500, 6))
  expect_equal(sim$visit_gaps, rep(1 / floor(6 * 500^0.16), 5))
  expect_false(anyNA(sim$Y[, 1]))
  obs <- !is.na(sim$Y)
  expect_true(all(apply(obs, 1, function(o) all(diff(o) <= 0))))  # monotone dropout
  expect_true(all(sim$entry_age <= sim$latent$death_age))         # alive at entry
  expect_true(all(rowSums(obs) - 1 <= sim$latent$dropout))
  expect_equal(sim$mean_visits, mean(rowSums(obs)))
  expect_equal(sim$Y[obs], sim$Y_complete[obs])
  expect_true(all(sim$Y > 0, na.rm = TRUE))
  sim3 <- simulate_cohort(n = 200, visit_gaps = "unequal", truncation = "informative", seed = 2)
  expect_equal(sim3$visit_gaps, c(0.1, 0.1, 0.12, 0.1, 0.12))
  expect_identical(sim3$settings$truncation, "informative")
  simc <- simulate_cohort(n = 200, visit_gaps = c(0.1, 0.2), seed = 3)
  expect_equal(ncol(simc$Y), 3)
})

test_that("a seed makes the simulation reproducible without touching the global RNG", {
  set.seed(10)
  before <- .Random.seed
  a <- simulate_cohort(n = 100, seed = 5)
  expect_identical(.Random.seed, before)
  b <- simulate_cohort(n = 100, seed = 5)
  expect_identical(a$Y, b$Y)
})

test_that("invalid simulation settings are rejected", {
  expect_error(simulate_cohort(n = 100, trunc_par = c(min = 2, max = 1)), "min < max")
  expect_error(simulate_cohort(n = 100, death_par = c(4.7, 0.85)), "named numeric")
  expect_error(simulate_cohort(n = 5000, pop_size = 2000), "increase `pop_size`")
  expect_error(simulate_cohort(n = 100, visit_gaps = c(0.1, -1)), "positive numbers")
  expect_error(simulate_cohort(n = 100, shape = function(t) 1), "vectorized")
})

test_that("sim_truth() agrees with Monte Carlo", {
  set.seed(1)
  N <- 4e5
  Z <- rlnorm(N, meanlog = log(1) - 0.5 * log1p(0.25), sdlog = sqrt(log1p(0.25)))
  Tdeath <- rweibull(N, shape = 4.7, scale = 0.85 * Z^(-1 / 4.7))
  tt <- c(0, 0.5, 0.9)
  mc <- sapply(tt, function(t) mean(Z[Tdeath >= t]))
  tr <- sim_truth(tt)
  expect_equal(tr$mu, mc, tolerance = 0.01)
  expect_equal(tr$mu[1], 1, tolerance = 1e-8)
  expect_equal(tr$bias, tr$mu / tr$mu[1])
  expect_equal(tr$mean, tr$mu[1] * shape_monotone(tt))

  # informative truncation: E(Z | A = t, T >= t), estimated with a narrow window
  A <- rweibull(N, shape = 3.2, scale = 1.2 * Z^(-1 / 3.2)) - 0.3
  mc2 <- sapply(c(0.3, 0.6), function(t) {
    keep <- abs(A - t) < 0.01 & Tdeath >= t
    mean(Z[keep])
  })
  tr2 <- sim_truth(c(0.3, 0.6), truncation = "informative")
  expect_equal(tr2$mu, mc2, tolerance = 0.03)
})
