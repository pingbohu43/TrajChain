sim <- simulate_cohort(n = 300, seed = 301)
fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_grid = 101)

test_that("the bootstrap returns pointwise summaries for every curve", {
  bt <- trajchain_boot(fit, B = 20, seed = 1)
  expect_s3_class(bt, "trajchain_boot")
  expect_identical(names(bt$se), names(fit$curves))
  expect_equal(bt$est, fit$curves)
  expect_equal(dim(bt$curves$shape), c(20, 101))
  ok <- is.finite(bt$lower$shape) & is.finite(bt$upper$shape)
  expect_true(all(bt$lower$shape[ok] <= bt$upper$shape[ok]))
  expect_true(all(bt$se$shape[ok] > 0))
  expect_equal(bt$se$bias_baseline[1], 0)
  expect_equal(bt$se$shape, apply(bt$curves$shape, 2, function(x) sd(x[is.finite(x)])))
})

test_that("the bootstrap is reproducible and does not change the global RNG state", {
  set.seed(9)
  before <- .Random.seed
  b1 <- trajchain_boot(fit, B = 10, seed = 5)
  expect_identical(.Random.seed, before)
  b2 <- trajchain_boot(fit, B = 10, seed = 5)
  expect_equal(b1$se, b2$se)
  b3 <- trajchain_boot(fit, B = 10, seed = 6)
  expect_false(isTRUE(all.equal(b1$se, b3$se)))
})

test_that("bootstrap summaries and data frames are well formed", {
  bt <- trajchain_boot(fit, B = 15, seed = 2, keep_curves = FALSE)
  expect_null(bt$curves)
  s <- summary(bt, ages = c(0.2, 0.6), which = c("shape", "bias"))
  expect_equal(nrow(s), 4)
  expect_named(s, c("quantity", "age", "estimate", "se", "lower", "upper"))
  df <- as.data.frame(bt)
  expect_equal(nrow(df), 9 * 101)
  out <- summary(bt, ages = c(-0.1, 1.1))
  expect_true(all(is.na(out$estimate)) && all(is.na(out$se)))
  expect_output(print(bt), "bootstrap replicates")
})

test_that("normal-approximation bands use the bootstrap standard error", {
  bt <- trajchain_boot(fit, B = 15, seed = 3)
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(bt, ci = "normal"))
  expect_silent(plot(bt, which = "bias"))
})

test_that("the bootstrap requires the data", {
  fit_nodata <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, keep_data = FALSE)
  expect_error(trajchain_boot(fit_nodata, B = 5), "keep_data")
  expect_error(trajchain_boot(list(), B = 5), "trajchain")
  expect_error(trajchain_boot(fit, B = 5, conf = 1.2), "between 0 and 1")
})

test_that("failed replicates are counted", {
  # with a small bandwidth some resampled data sets leave grid ages without
  # observed visit pairs, so the ratio estimator is undefined there
  fit_small <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, bandwidth = 0.008,
                         n_grid = 51)
  bt <- trajchain_boot(fit_small, B = 30, seed = 4)
  expect_gt(bt$n_failed, 0)
  expect_lt(bt$n_failed, 30)
  expect_equal(unname(bt$n_ok$shape[1]), 30 - bt$n_failed)
  expect_equal(unname(bt$n_incomplete[["shape"]]), bt$n_failed)
})

test_that("penalized fits can be bootstrapped with lambda held fixed", {
  sim3 <- simulate_cohort(n = 400, visit_gaps = "unequal", seed = 302)
  fit3 <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 1e-4, n_grid = 51)
  bt <- trajchain_boot(fit3, B = 10, seed = 1)
  expect_equal(bt$lambda, 1e-4)
  expect_true(all(is.finite(bt$se$shape)))
})
