sim <- simulate_cohort(n = 400, seed = 201)

test_that("cv_bandwidth() returns a valid choice and is reproducible", {
  cvb <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 3)
  expect_s3_class(cvb, "trajchain_cv")
  expect_identical(cvb$type, "bandwidth")
  expect_true(cvb$constant %in% cvb$candidates)
  expect_equal(cvb$bandwidth, cvb$constant * bw_rot(sim$entry_age))
  expect_equal(cvb$cv[which(cvb$candidates == cvb$constant)], min(cvb$cv))
  expect_equal(cvb, cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 3),
               ignore_attr = TRUE, ignore_function_env = TRUE)
  expect_output(print(cvb), "selected constant")
})

test_that("CV leaves the global random number stream untouched when a seed is given", {
  set.seed(42)
  before <- .Random.seed
  cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 1, candidates = c(0.5, 1, 2))
  expect_identical(.Random.seed, before)
})

test_that("repeated partitions are averaged", {
  cv3 <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 3, n_repeats = 3,
                      candidates = c(0.5, 1, 2))
  expect_equal(dim(cv3$cv_by_repeat), c(3, 3))
  expect_equal(cv3$cv, colMeans(cv3$cv_by_repeat))
  cv1 <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 3,
                      candidates = c(0.5, 1, 2))
  expect_equal(cv3$cv_by_repeat[1, ], cv1$cv)
})

test_that("trajchain(bandwidth = 'cv') uses the cross-validated constant", {
  cvb <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  fit1 <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = "cv")
  fit2 <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = cvb)
  expect_equal(fit1$bandwidth[["ratio"]], cvb$bandwidth)
  expect_equal(fit1$curves, fit2$curves)
  expect_equal(fit1$bandwidth_constant, cvb$constant)
})

test_that("invalid CV settings are rejected", {
  expect_error(cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, candidates = c(-1, 1)),
               "positive")
  expect_error(cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_folds = 1),
               "whole number")
  expect_error(cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, seed = 1.5), "seed")
})

test_that("cv_lambda() selects a penalty for unequally spaced visits", {
  sim3 <- simulate_cohort(n = 600, visit_gaps = "unequal", seed = 202)
  grid <- c(1e-6, 1e-5, 1e-4, 1e-3)
  cvl <- cv_lambda(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda_grid = grid, seed = 2)
  expect_s3_class(cvl, "trajchain_cv")
  expect_identical(cvl$type, "lambda")
  expect_true(cvl$lambda %in% grid)
  expect_equal(dim(cvl$scores), c(4, 5))
  expect_output(print(cvl), "selected lambda")
  fit <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = cvl)
  expect_equal(fit$lambda, cvl$lambda)
  expect_identical(fit$cv$lambda, cvl)
  expect_error(trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1,
                         lambda = cv_bandwidth(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1)),
               "cv_lambda")
  expect_error(cv_lambda(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda_grid = -1),
               "non-negative")
})
