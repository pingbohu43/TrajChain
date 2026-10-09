# Regression tests against reference values computed with the authors' original
# research scripts (see data-raw/regression-fixtures.R). They guarantee that
# TrajChain reproduces the estimators used for the results in the paper.
fx <- readRDS(test_path("fixtures", "original-code.rds"))

test_that("Case 1 (equally spaced visits) reproduces compute_fhat()", {
  d <- fx$case1
  fit <- trajchain(d$Y, d$A, d$v, tau0 = 0, tau1 = 1)
  expect_equal(unname(fit$bandwidth), c(d$h1, d$h2), tolerance = 1e-14)
  expect_equal(fit$grid$shape, d$fhat_grid, tolerance = 1e-12)
  expect_equal(fit$grid$ratio[-1], d$rhat, tolerance = 1e-12)
  expect_equal(fit$curves$shape[d$idx], d$fhat, tolerance = 1e-12)
  expect_equal(fit$curves$mu_pooled[d$idx], d$g, tolerance = 1e-12)
  expect_equal(fit$curves$mean_pooled[d$idx], d$EY, tolerance = 1e-12)
  expect_equal(fit$mu0[["pooled"]], d$EZ0, tolerance = 1e-12)
})

test_that("bandwidth CV reproduces cv_select_c() of the simulation code", {
  d <- fx$case1
  cvb <- cv_bandwidth(d$Y, d$A, d$v, tau0 = 0, tau1 = 1, seed = 99)
  expect_equal(cvb$cv, d$cv, tolerance = 1e-12)
  expect_equal(cvb$cv_folds, d$cv_folds, tolerance = 1e-12)
  expect_equal(cvb$constant, d$c_opt)
})

test_that("Case 3 (unequally spaced visits) reproduces optimize_fhat_with_pairs_cvxr()", {
  d <- fx$case3
  fit <- trajchain(d$Y, d$A, d$gaps, tau0 = 0, tau1 = 1, lambda = d$lambda)
  expect_equal(unname(fit$bandwidth), c(d$h1, d$h2), tolerance = 1e-14)
  expect_equal(fit$pairs$n, d$counts)
  expect_equal(fit$grid$shape, d$fhat_grid, tolerance = 1e-8)
  expect_equal(fit$curves$shape[d$idx], d$fhat, tolerance = 1e-8)
  expect_equal(fit$curves$mu_pooled[d$idx], d$g, tolerance = 1e-8)
  expect_equal(fit$curves$mean_pooled[d$idx], d$EY, tolerance = 1e-8)
})

test_that("penalty CV reproduces cv_select_lambda()", {
  d <- fx$case3
  cvl <- cv_lambda(d$Y, d$A, d$gaps, tau0 = 0, tau1 = 1, lambda_grid = d$lambda_grid, seed = 5)
  expect_equal(cvl$scores, d$cv_scores, tolerance = 1e-7)
  expect_equal(cvl$lambda, d$best_lambda)
})

test_that("bandwidth CV with unequal gaps reproduces the Case 3 cv_select_c()", {
  d <- fx$case3
  cvb <- cv_bandwidth(d$Y, d$A, d$gaps, tau0 = 0, tau1 = 1, seed = 77)
  expect_equal(cvb$cv, d$cv_c, tolerance = 1e-12)
  expect_equal(cvb$constant, d$c_opt)
})

test_that("the data-analysis pipeline is reproduced on the age scale", {
  d <- fx$data_analysis
  fit <- trajchain(d$Y, d$A, visit_gaps = 16 / 12, tau0 = d$tau0, m = d$m)
  expect_equal(unname(fit$bandwidth), c(d$h1, d$h2), tolerance = 1e-14)
  cur <- fit$curves[d$idx, ]
  expect_equal(cur$shape, d$fhat, tolerance = 1e-12)
  for (v in c("bias_baseline", "mean_baseline", "nw_baseline", "bias_pooled",
              "mean_pooled", "nw_pooled")) {
    expect_equal(cur[[v]], d[[v]], tolerance = 1e-12, label = v)
  }
  expect_equal(fit$mu0, d$EZ0, tolerance = 1e-12)
})

test_that("data-analysis CV with full-sample scale reproduces cv_select_c_data_analysis()", {
  d <- fx$data_analysis
  cvb <- cv_bandwidth(d$Y, d$A, 16 / 12, tau0 = d$tau0, m = d$m, seed = 2026,
                      n_repeats = 2, sd_scope = "full")
  expect_equal(cvb$cv, d$cv, tolerance = 1e-12)
  expect_equal(cvb$cv_by_repeat, d$cv_by_repeat, tolerance = 1e-12)
  expect_equal(cvb$constant, d$c_opt)
})

test_that("the bootstrap reproduces bootstrap_fhat_data_analysis()", {
  d <- fx$data_analysis
  fit <- trajchain(d$Y, d$A, visit_gaps = 16 / 12, tau0 = d$tau0, m = d$m)
  bt <- trajchain_boot(fit, B = d$boot$B, seed = d$boot$seed)
  map <- c(fhat = "shape", EY = "mean_baseline", bias = "bias_baseline",
           NW = "nw_baseline")
  for (k in names(map)) {
    expect_equal(bt$se[[map[[k]]]][d$idx], d$boot$se[[k]], tolerance = 1e-10, label = k)
    expect_equal(bt$lower[[map[[k]]]][d$idx], d$boot$lower[[k]], tolerance = 1e-10, label = k)
    expect_equal(bt$upper[[map[[k]]]][d$idx], d$boot$upper[[k]], tolerance = 1e-10, label = k)
    expect_equal(unname(bt$n_incomplete[[map[[k]]]]), unname(d$boot$n_incomplete[[k]]))
  }
})

test_that("simulate_cohort() reproduces the original data generators", {
  for (case in c("case1", "case3")) {
    d <- fx[[case]]
    sim <- simulate_cohort(n = d$n, shape = if (case == "case1") "monotone" else "ushape",
                           visit_gaps = if (case == "case1") "equal" else "unequal",
                           seed = d$seed)
    expect_equal(sim$entry_age, d$A, tolerance = 0, label = case)
    expect_equal(unname(sim$Y), unname(d$Y), tolerance = 1e-12, label = case)
  }
  g <- fx$gen_informative_unequal
  sim <- simulate_cohort(n = 300, shape = "ushape", visit_gaps = "unequal",
                         truncation = "informative", seed = g$seed)
  expect_equal(sim$entry_age, g$A, tolerance = 0)
  expect_equal(sim$latent$dropout, g$C)
  expect_equal(unname(sim$Y), unname(g$Y), tolerance = 1e-12)
})
