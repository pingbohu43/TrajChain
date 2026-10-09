sim <- simulate_cohort(n = 400, seed = 401)
fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_grid = 201)

test_that("print and summary describe the fit", {
  expect_output(print(fit), "TrajChain fit")
  s <- summary(fit, ages = c(0.1, 0.5))
  expect_s3_class(s, "summary.trajchain")
  expect_equal(s$table$age, c(0.1, 0.5))
  expect_output(print(s), "Summary of TrajChain fit")
  s2 <- summary(fit, mu_method = "pooled")
  expect_identical(s2$mu_method, "pooled")
  expect_equal(nrow(s2$table), 6)
})

test_that("predict() reproduces the curves and interpolates the shape", {
  ages <- fit$curves$t[c(1, 37, 101, 201)]
  for (type in c("shape", "mean", "bias", "mu", "nw")) {
    for (mm in c("baseline", "pooled")) {
      col <- if (type == "shape") "shape" else paste0(type, "_", mm)
      expect_equal(predict(fit, ages, type = type, mu_method = mm),
                   fit$curves[[col]][c(1, 37, 101, 201)], tolerance = 1e-12,
                   label = paste(type, mm))
    }
  }
  expect_equal(predict(fit, type = "bias"), fit$curves$bias_baseline)
  expect_true(all(is.na(predict(fit, c(-0.1, 1.1), type = "mean"))))
  expect_true(all(is.na(predict(fit, c(-0.1, 1.1)))))
  fit_nd <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_grid = 201,
                      keep_data = FALSE)
  expect_equal(predict(fit_nd, ages, type = "nw"), fit$curves$nw_baseline[c(1, 37, 101, 201)])
  expect_error(predict(fit, "a"), "finite ages")
})

test_that("shape_fun() interpolates the grid values", {
  expect_equal(fit$shape_fun(fit$grid$t), fit$grid$shape)
  expect_true(is.na(fit$shape_fun(1.2)))
})

test_that("as.data.frame() returns the curves", {
  expect_identical(as.data.frame(fit), fit$curves)
})

test_that("plots are produced without error", {
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(fit))
  expect_silent(plot(fit, which = c("shape", "bias"), mu_method = "pooled",
                     truth = sim_truth(seq(0, 1, by = 0.05))))
  fit2 <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, n_grid = 101,
                    bandwidth = 0.15)
  expect_silent(plot_trajectories(list(A = fit, B = fit2), which = "shape",
                                  col = c("red", "blue")))
  bt <- trajchain_boot(fit, B = 10, seed = 1)
  expect_silent(plot_trajectories(list(Boot = bt, Fit = fit2), which = c("mean", "nw")))
  expect_silent(plot(cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1,
                                  candidates = c(0.5, 1, 2))))
  expect_error(plot_trajectories(list(1, 2)), "trajchain")
  expect_error(plot_trajectories(fit, truth = 1:3), "truth")
})

test_that("simulated-data objects print", {
  expect_output(print(sim), "simulated cohort")
})
