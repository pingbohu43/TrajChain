sim <- simulate_cohort(n = 300, seed = 101)

test_that("invalid data are rejected with informative errors", {
  Y <- sim$Y
  A <- sim$entry_age
  v <- sim$visit_gaps
  expect_error(trajchain(Y[, 1, drop = FALSE], A, v, 0, 1), "at least two columns")
  Yneg <- Y
  Yneg[1, 1] <- -1
  expect_error(trajchain(Yneg, A, v, 0, 1), "non-negative")
  Yna <- Y
  Yna[2, 1] <- NA
  expect_error(trajchain(Yna, A, v, 0, 1), "baseline")
  Ychr <- Y
  storage.mode(Ychr) <- "character"
  expect_error(trajchain(Ychr, A, v, 0, 1), "numeric matrix")
  Yinf <- Y
  Yinf[3, 2] <- Inf
  expect_error(trajchain(Yinf, A, v, 0, 1), "infinite")
  expect_error(trajchain(Y, A[-1], v, 0, 1), "one value per row")
  expect_error(trajchain(Y, replace(A, 1, NA), v, 0, 1), "missing or infinite")
  expect_error(trajchain(Y, A, c(v[1], v[1]), 0, 1), "length 1")
  expect_error(trajchain(Y, A, -1, 0, 1), "positive")
  expect_error(trajchain(Y, A, v, 0), "tau1")
  expect_error(trajchain(Y, A, v, 0, tau1 = 0.97), "integer multiple")
  expect_error(trajchain(Y, A, v, 0, tau1 = 1, m = 17), "inconsistent")
  expect_error(trajchain(Y, A, v, 1, tau1 = 0.5), "larger than")
  expect_error(trajchain(Y, A, v[1] * 1.5, 0, 1, grid_step = v[1]), "integer multiple")
})

test_that("the grid can be given through tau1 or m", {
  v <- sim$visit_gaps[1]
  f1 <- trajchain(sim$Y, sim$entry_age, v, tau0 = 0, tau1 = 1)
  f2 <- trajchain(sim$Y, sim$entry_age, v, tau0 = 0, m = round(1 / v))
  expect_equal(f1$curves, f2$curves)
  expect_equal(f1$design$m, round(1 / v))
  expect_equal(f1$grid$t, (seq_len(f1$design$m) - 0.5) / f1$design$m)
})

test_that("a data frame is accepted for Y", {
  f1 <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  f2 <- trajchain(as.data.frame(sim$Y), sim$entry_age, sim$visit_gaps, 0, 1)
  expect_equal(f1$curves, f2$curves)
})

test_that("the method is chosen from the visit design", {
  f_eq <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1)
  expect_identical(f_eq$method, "chain")
  expect_error(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, method = "chain",
                         grid_step = sim$visit_gaps[1] / 2), "requires every visit gap")
  # widely spaced visits: refining the grid switches to the penalized estimator
  f_wide <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1,
                      grid_step = sim$visit_gaps[1] / 2, lambda = 1e-4)
  expect_identical(f_wide$method, "penalized")
  expect_equal(f_wide$design$m, 2 * f_eq$design$m)
  expect_warning(trajchain(sim$Y, sim$entry_age, sim$visit_gaps, 0, 1, lambda = 1),
                 "ignored by the chain estimator")
})

test_that("intermittent missing visits trigger a warning and are handled pairwise", {
  Y <- sim$Y
  rows <- which(rowSums(!is.na(Y)) >= 4)[1:5]
  Y[rows, 2] <- NA
  expect_warning(fit <- trajchain(Y, sim$entry_age, sim$visit_gaps, 0, 1),
                 "intermittent")
  expect_equal(fit$design$n_intermittent, 5)
  expect_true(all(is.finite(fit$grid$shape)))
})

test_that("the greatest common divisor of visit gaps defines the grid", {
  expect_equal(TrajChain:::.gcd_numeric(c(0.1, 0.1, 0.12)), 0.02, tolerance = 1e-12)
  expect_equal(TrajChain:::.gcd_numeric(c(16, 16) / 12), 16 / 12)
  expect_equal(TrajChain:::.gcd_numeric(c(6, 9, 15)), 3)
  sim3 <- simulate_cohort(n = 300, visit_gaps = "unequal", seed = 1)
  fit <- trajchain(sim3$Y, sim3$entry_age, sim3$visit_gaps, 0, 1, lambda = 1e-4)
  expect_equal(fit$design$m, 50L)
  expect_equal(fit$design$grid_step, 0.02)
})

test_that("visit_matrix() reshapes long data", {
  long <- data.frame(id = c("b", "b", "a", "a", "a", "c"),
                     visit = c(1, 2, 1, 2, 3, 1),
                     y = c(1, 2, 3, 4, 5, 6))
  M <- visit_matrix(long, "id", "visit", "y")
  expect_equal(dim(M), c(3, 3))
  expect_equal(rownames(M), c("b", "a", "c"))
  expect_equal(colnames(M), paste0("visit_", 1:3))
  expect_equal(unname(M["a", ]), c(3, 4, 5))
  expect_equal(unname(M["c", ]), c(6, NA, NA))
  expect_error(visit_matrix(long, "id", "visit", "nope"), "not found")
  expect_error(visit_matrix(rbind(long, long[1, ]), "id", "visit", "y"), "more than one row")
  long$visit[1] <- 1.5
  expect_error(visit_matrix(long, "id", "visit", "y"), "whole visit numbers")
})

test_that("visit_matrix() and trajchain() agree with the example data", {
  d <- subset(simcohort, group == "A")
  Y <- visit_matrix(d, "id", "visit", "marker")
  expect_equal(nrow(Y), 300)
  expect_equal(sum(!is.na(Y)), nrow(d))
  A <- d$entry_age[match(rownames(Y), d$id)]
  fit <- trajchain(Y, A, visit_gaps = 16 / 12, tau0 = 30, m = 22)
  expect_equal(fit$design$tau1, 30 + 22 * 16 / 12)
  expect_equal(sum(fit$grid$shape) * fit$design$grid_step, 1)
})
