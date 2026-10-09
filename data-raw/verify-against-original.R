# Numerical comparison of TrajChain with the original research code.
#
# Usage (from the package root, with TrajChain installed):
#   TRAJCHAIN_CODE_DIR=/path/to/biometrika-revision/code Rscript data-raw/verify-against-original.R
# where the code directory contains independent-truncation-code/,
# dependent-truncation-code/ and data-analysis/. Every line of the output should
# end in "OK". Needs the Matrix package (used by the original Case 3 code).
CODE <- Sys.getenv("TRAJCHAIN_CODE_DIR", "..")
suppressMessages(library(Matrix))
suppressMessages(library(TrajChain))

# Evaluate the function definitions of an original script (library(), source()
# and top-level analysis code are skipped), optionally the whole script.
load_defs <- function(file, env, all = FALSE) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[!grepl("^\\s*(library\\(|source\\(|FUN_PATH|CV_PATH)", lines)]
  for (e in parse(text = lines)) {
    is_fun <- is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
    is_shape <- is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) == "mon_shape"
    if (all || is_fun || is_shape) eval(e, env)
  }
  invisible(env)
}
IND <- file.path(CODE, "independent-truncation-code")
DEP <- file.path(CODE, "dependent-truncation-code")
DA <- file.path(CODE, "data-analysis")
env_sim <- new.env()
load_defs(file.path(IND, "assumption-B-functions.R"), env_sim, all = TRUE)
load_defs(file.path(IND, "assumption-B-cv-lambda.R"), env_sim, all = TRUE)
env_dep <- new.env()
load_defs(file.path(DEP, "functions.R"), env_dep, all = TRUE)
env_da <- new.env()
load_defs(file.path(DA, "data-analysis-Case-1-including-birth-cohort-CV-bandwidth.R"), env_da)
load_defs(file.path(DA, "data-analysis-Case-1-including-birth-cohort-bootstrap.R"), env_da)
env_da$n_dense <- 801
env_c1 <- new.env(parent = env_sim)
load_defs(file.path(IND, "assumption-B-Case-1-Mon-CV-small-bandwidth.R"), env_c1)
env_c3 <- new.env(parent = env_sim)
load_defs(file.path(IND, "assumption-B-Case-3-Mon-CV-both-bandwidth-lambda.R"), env_c3)

maxrel <- function(a, b) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  stopifnot(length(a) == length(b))
  if (!identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(a)
  if (!any(ok)) return(0)
  max(abs(a[ok] - b[ok]) / pmax(1e-300, abs(b[ok]), 1e-8))
}
report <- function(label, val, tol = 1e-10) {
  cat(sprintf("%-70s max rel diff = %.3e  %s\n", label, val, if (val <= tol) "OK" else "**** FAIL ****"))
  invisible(val <= tol)
}

P <- list(n = 10000, k = 5, z1 = 1, z2 = 0.25, e1 = 1, e2 = 0.05, t1 = 4.7, t2 = 0.85,
          c1 = 3, c2 = 3)
mon <- env_sim$mon_shape
ush <- env_sim$U_quad_pdf

cat("\n=== A. Data generators ===\n")
for (SEED in c(1, 17)) for (rs in c(1000, 2000)) {
  m <- floor(6 * rs^0.16); v <- 1 / m
  g0 <- env_sim$generate_data(SEED, P$n, rs, mon, P$k, v, 0, P$z1, P$z2, -0.3, 1.9,
                              P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  s0 <- simulate_cohort(n = rs, shape = mon, visit_gaps = "equal", seed = SEED)
  ind <- outer(g0$A, (0:P$k) * v, "+") <= matrix(g0$T_tilde, rs, P$k + 1)
  Yo <- g0$Y; Yo[!ind] <- NA
  report(sprintf("random/equal  seed %d n %d: A, Z, T, C", SEED, rs),
         max(maxrel(s0$entry_age, g0$A), maxrel(s0$latent$Z, g0$Z),
             maxrel(s0$latent$death_age, g0$T.to.death), maxrel(s0$latent$dropout, g0$C)), 0)
  report(sprintf("random/equal  seed %d n %d: observed Y (with NA pattern)", SEED, rs),
         maxrel(unname(s0$Y), Yo), 0)
  report(sprintf("random/equal  seed %d n %d: mean visits", SEED, rs), maxrel(s0$mean_visits, g0$ave_visit_num), 0)
}
for (SEED in c(3, 4)) {
  g0 <- env_sim$generate_data_different_visit_interval(SEED, P$n, 1000, ush, P$k, 0.02, 0, 5, 6,
         P$z1, P$z2, -0.3, 1.9, P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  s0 <- simulate_cohort(n = 1000, shape = ush, visit_gaps = "unequal", seed = SEED)
  ind <- g0$visit_age <= matrix(g0$T_tilde, 1000, P$k + 1)
  Yo <- g0$Y; Yo[!ind] <- NA
  report(sprintf("random/unequal seed %d: A, C, observed Y", SEED),
         max(maxrel(s0$entry_age, g0$A), maxrel(s0$latent$dropout, g0$C), maxrel(unname(s0$Y), Yo)), 0)
}
for (SEED in c(5, 6)) {
  m <- floor(6 * 1000^0.16); v <- 1 / m
  g0 <- env_dep$generate_data(SEED, P$n, 1000, mon, P$k, v, 0.3, P$z1, P$z2, 3.2, 1.2,
                              P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  s0 <- simulate_cohort(n = 1000, shape = mon, truncation = "informative", seed = SEED)
  ind <- outer(g0$A, (0:P$k) * v, "+") <= matrix(g0$T_tilde, 1000, P$k + 1)
  Yo <- g0$Y; Yo[!ind] <- NA
  report(sprintf("informative/equal seed %d: A, T, observed Y", SEED),
         max(maxrel(s0$entry_age, g0$A), maxrel(s0$latent$death_age, g0$T.to.death),
             maxrel(unname(s0$Y), Yo)), 0)
  g1 <- env_dep$generate_data_different_visit_interval(SEED, P$n, 1000, mon, P$k, 0.02, 0.3, 5, 6,
         P$z1, P$z2, 3.2, 1.2, P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  s1 <- simulate_cohort(n = 1000, shape = mon, visit_gaps = "unequal", truncation = "informative", seed = SEED)
  ind <- g1$visit_age <= matrix(g1$T_tilde, 1000, P$k + 1)
  Yo <- g1$Y; Yo[!ind] <- NA
  report(sprintf("informative/unequal seed %d: observed Y", SEED), maxrel(unname(s1$Y), Yo), 0)
}

cat("\n=== B. Case 1 (equal spacing), simulation pipeline compute_fhat() ===\n")
for (SEED in c(1, 2)) for (rs in c(1000, 2000)) {
  m <- floor(6 * rs^0.16)
  gen <- env_sim$generate_data(SEED, P$n, rs, mon, P$k, 1 / m, 0, P$z1, P$z2, -0.3, 1.9,
                               P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  sdA <- sd(gen$A)
  h1 <- env_sim$rot_bw_from_A(sdA, real_size = rs, -1/6, kernel = "epanechnikov")
  h2 <- env_sim$rot_bw_from_A(sdA, real_size = rs, -1/5, kernel = "epanechnikov")
  out <- env_sim$compute_fhat(SEED = SEED, m = m, h1 = h1, h2 = 2.34 * h2, n = P$n, real_size = rs,
                              f = mon, k = P$k, delta = 0, z1 = P$z1, z2 = P$z2, a1 = -0.3, a2 = 1.9,
                              t1 = P$t1, t2 = P$t2, e1 = P$e1, e2 = P$e2, c1 = P$c1, c2 = P$c2,
                              Ker = env_sim$K_epan, return_gen = FALSE, n_dense = 801)
  s0 <- simulate_cohort(n = rs, shape = mon, seed = SEED)
  fit <- trajchain(s0$Y, s0$entry_age, s0$visit_gaps, tau0 = 0, tau1 = 1)   # default bandwidths
  report(sprintf("seed %d n %d: default bandwidths equal h1, 2.34*h2", SEED, rs),
         maxrel(fit$bandwidth, c(h1, 2.34 * h2)), 1e-14)
  report(sprintf("seed %d n %d: grid fhat", SEED, rs), maxrel(fit$grid$shape, out$results$fhat))
  report(sprintf("seed %d n %d: grid ratios", SEED, rs), maxrel(fit$grid$ratio[-1], out$results$rhat[-1]))
  report(sprintf("seed %d n %d: dense fhat", SEED, rs), maxrel(fit$curves$shape, out$dense$fhat))
  report(sprintf("seed %d n %d: dense g = mu_pooled", SEED, rs), maxrel(fit$curves$mu_pooled, out$dense$g))
  report(sprintf("seed %d n %d: dense EY = mean_pooled", SEED, rs), maxrel(fit$curves$mean_pooled, out$dense$EY))
  report(sprintf("seed %d n %d: EZ0", SEED, rs), maxrel(fit$mu0[["pooled"]], out$EZ0))
}

cat("\n=== C. Case 1 data-analysis pipeline compute_fhat_data_analysis() on the age scale ===\n")
for (SEED in c(11, 12)) {
  m <- 22; tau0 <- 30; v <- 16 / 12; tau1 <- tau0 + m * v
  s0 <- simulate_cohort(n = 300, shape = mon, visit_gaps = rep(1 / 22, 2), truncation = "informative",
                        seed = SEED)
  A <- tau0 + (tau1 - tau0) * s0$entry_age
  Y <- unname(s0$Y)
  h1 <- env_da$rot_bw_from_A(A, -1/6); h2 <- 2.34 * env_da$rot_bw_from_A(A, -1/5)
  res <- env_da$compute_fhat_data_analysis(Y, h1, h2, Ker = env_da$K_epan, m = m, A = A,
                                           tau0 = tau0, tau1 = tau1)
  fit <- trajchain(Y, A, visit_gaps = v, tau0 = tau0, m = m)
  sb <- res$survival_birth_cohort_result; nb <- res$result_no_birth_cohort
  report(sprintf("seed %d: bandwidths", SEED), maxrel(fit$bandwidth, c(h1, h2)), 1e-14)
  report(sprintf("seed %d: t grid", SEED), maxrel(fit$curves$t, sb$t), 1e-14)
  report(sprintf("seed %d: fhat", SEED), maxrel(fit$curves$shape, sb$fhat))
  report(sprintf("seed %d: survival_birth_cohort_bias = bias_baseline", SEED), maxrel(fit$curves$bias_baseline, sb$survival_birth_cohort_bias))
  report(sprintf("seed %d: EY_survival_birth_cohort = mean_baseline", SEED), maxrel(fit$curves$mean_baseline, sb$EY_survival_birth_cohort))
  report(sprintf("seed %d: NW_survival_birth_cohort = nw_baseline", SEED), maxrel(fit$curves$nw_baseline, sb$NW_survival_birth_cohort))
  report(sprintf("seed %d: surv_bias = bias_pooled", SEED), maxrel(fit$curves$bias_pooled, nb$surv_bias))
  report(sprintf("seed %d: EY = mean_pooled", SEED), maxrel(fit$curves$mean_pooled, nb$EY))
  report(sprintf("seed %d: NW_no_birth_cohort = nw_pooled", SEED), maxrel(fit$curves$nw_pooled, nb$NW_no_birth_cohort))
  report(sprintf("seed %d: EZ0 (both)", SEED), maxrel(fit$mu0, c(res$EZ0_survival_birth_cohort, res$EZ0)))
  report(sprintf("seed %d: g_fun_baseline at 3 ages = mu_baseline via predict", SEED),
         maxrel(predict(fit, c(33, 45, 57), type = "mu"), sapply(c(33, 45, 57), res$g_fun_baseline)))
}

cat("\n=== D. Bandwidth cross-validation ===\n")
for (SEED in c(1, 2)) {
  rs <- 1000; m <- floor(6 * rs^0.16); cand <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4)
  gen <- env_sim$generate_data(SEED, P$n, rs, mon, P$k, 1 / m, 0, P$z1, P$z2, -0.3, 1.9,
                               P$t1, P$t2, P$e1, P$e2, P$c1, P$c2)
  cv0 <- env_c1$cv_select_c(gen, m = m, cand_c = cand, n_folds = 5L, Ker = env_sim$K_epan,
                            fold_seed = 1000000L + SEED)
  s0 <- simulate_cohort(n = rs, shape = mon, seed = SEED)
  cv1 <- cv_bandwidth(s0$Y, s0$entry_age, s0$visit_gaps, tau0 = 0, tau1 = 1, candidates = cand,
                      seed = 1000000L + SEED, sd_scope = "fold")
  report(sprintf("Case 1 sim cv_select_c seed %d: CV(c)", SEED), maxrel(cv1$cv, cv0$cv))
  report(sprintf("Case 1 sim cv_select_c seed %d: per-fold CV", SEED), maxrel(cv1$cv_folds, cv0$cv_folds))
  report(sprintf("Case 1 sim cv_select_c seed %d: c_opt, n_dropped", SEED),
         max(maxrel(cv1$constant, cv0$c_opt), maxrel(cv1$table$n_dropped, cv0$n_dropped)), 0)
}
for (SEED in c(11, 12)) for (nr in c(1, 3)) {
  m <- 22; tau0 <- 30; v <- 16 / 12; tau1 <- tau0 + m * v
  s0 <- simulate_cohort(n = 300, shape = mon, visit_gaps = rep(1 / 22, 2), truncation = "informative", seed = SEED)
  A <- tau0 + (tau1 - tau0) * s0$entry_age; Y <- unname(s0$Y)
  cand <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4)
  cv0 <- env_da$cv_select_c_data_analysis(Y, A, m, tau0, tau1, cand, n_folds = 5, Ker = env_da$K_epan,
                                          fold_seed = 2026, n_repeats = nr)
  cv1 <- cv_bandwidth(Y, A, v, tau0 = tau0, m = m, candidates = cand, seed = 2026, n_repeats = nr,
                      sd_scope = "full")
  report(sprintf("data-analysis CV seed %d repeats %d: CV(c)", SEED, nr), maxrel(cv1$cv, cv0$cv))
  report(sprintf("data-analysis CV seed %d repeats %d: per-repeat CV", SEED, nr), maxrel(cv1$cv_by_repeat, cv0$cv_by_repeat))
  report(sprintf("data-analysis CV seed %d repeats %d: c_opt, n_dropped, n_test", SEED, nr),
         max(maxrel(cv1$constant, cv0$c_opt), maxrel(cv1$table$n_dropped, cv0$cv_table$n_dropped),
             maxrel(cv1$n_test, cv0$n_test_points)), 0)
}
for (SEED in c(3, 4)) {
  gen <- env_sim$generate_data_different_visit_interval(SEED = SEED, n = P$n, real_size = 1000, f = mon,
           k = P$k, v = 1 / 50, delta = 0, ch1 = 5, ch2 = 6, z1 = P$z1, z2 = P$z2, a1 = -0.3, a2 = 1.9,
           t1 = P$t1, t2 = P$t2, e1 = P$e1, e2 = P$e2, c1 = P$c1, c2 = P$c2)
  cand <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4)
  cv0 <- env_c3$cv_select_c(gen, m = 50, cand_c = cand, ch1 = 5, ch2 = 6, n_folds = 5L,
                            Ker = env_sim$K_epan, fold_seed = 1000000L + SEED, tau = 1)
  s0 <- simulate_cohort(n = 1000, shape = mon, visit_gaps = "unequal", seed = SEED)
  cv1 <- cv_bandwidth(s0$Y, s0$entry_age, s0$visit_gaps, tau0 = 0, tau1 = 1, candidates = cand,
                      seed = 1000000L + SEED)
  report(sprintf("Case 3 cv_select_c seed %d: CV(c)", SEED), maxrel(cv1$cv, cv0$cv))
  report(sprintf("Case 3 cv_select_c seed %d: per-fold CV", SEED), maxrel(cv1$cv_folds, cv0$cv_folds))
  report(sprintf("Case 3 cv_select_c seed %d: c_opt, n_dropped", SEED),
         max(maxrel(cv1$constant, cv0$c_opt), maxrel(cv1$table$n_dropped, cv0$n_dropped)), 0)
}

cat("\n=== E. Case 3 (unequal spacing) optimize_fhat_with_pairs_cvxr() ===\n")
for (SEED in c(3, 4)) for (lam in c(1e-6, 3.16e-5, 1e-3)) for (f_case in c("mon", "ush")) {
  ff <- if (f_case == "mon") mon else ush
  gen <- env_sim$generate_data_different_visit_interval(SEED = SEED, n = P$n, real_size = 1000, f = ff,
           k = P$k, v = 1 / 50, delta = 0, ch1 = 5, ch2 = 6, z1 = P$z1, z2 = P$z2, a1 = -0.3, a2 = 1.9,
           t1 = P$t1, t2 = P$t2, e1 = P$e1, e2 = P$e2, c1 = P$c1, c2 = P$c2)
  sdA <- sd(gen$A)
  h1 <- env_sim$rot_bw_from_A(sdA, real_size = 1000, -1/6, kernel = "epanechnikov")
  h2 <- env_sim$rot_bw_from_A(sdA, real_size = 1000, -1/5, kernel = "epanechnikov")
  out <- env_sim$optimize_fhat_with_pairs_cvxr(m = 50, w1 = 0, threshold = 0.6, SEED = SEED, h1 = h1,
           h2 = 2.34 * h2, Ker = env_sim$K_epan, n = P$n, real_size = 1000, f = ff, k = P$k, delta = 0,
           z1 = P$z1, z2 = P$z2, a1 = -0.3, a2 = 1.9, t1 = P$t1, t2 = P$t2, e1 = P$e1, e2 = P$e2,
           c1 = P$c1, c2 = P$c2, n_dense = 801, ridge = lam, ch1 = 5, ch2 = 6)
  s0 <- simulate_cohort(n = 1000, shape = ff, visit_gaps = "unequal", seed = SEED)
  fit <- trajchain(s0$Y, s0$entry_age, s0$visit_gaps, tau0 = 0, tau1 = 1, lambda = lam)
  lab <- sprintf("seed %d lambda %.2g %s", SEED, lam, f_case)
  report(paste(lab, ": counts n per pair"), maxrel(fit$pairs$n, c(out$number_5_vec, out$number_6_vec)), 0)
  report(paste(lab, ": grid fhat"), maxrel(fit$grid$shape, out$results$fhat), 1e-8)
  report(paste(lab, ": dense fhat"), maxrel(fit$curves$shape, out$dense$fhat), 1e-8)
  report(paste(lab, ": dense g = mu_pooled"), maxrel(fit$curves$mu_pooled, out$dense$g), 1e-8)
  report(paste(lab, ": dense EY = mean_pooled"), maxrel(fit$curves$mean_pooled, out$dense$EY), 1e-8)
}

cat("\n=== F. Penalty cross-validation cv_select_lambda() ===\n")
LG <- exp(seq(log(1e-6), log(1e-2), length.out = 25))
for (SEED in c(3, 8)) {
  gen <- env_sim$generate_data_different_visit_interval(SEED = SEED, n = P$n, real_size = 1000, f = mon,
           k = P$k, v = 1 / 50, delta = 0, ch1 = 5, ch2 = 6, z1 = P$z1, z2 = P$z2, a1 = -0.3, a2 = 1.9,
           t1 = P$t1, t2 = P$t2, e1 = P$e1, e2 = P$e2, c1 = P$c1, c2 = P$c2)
  h1 <- env_sim$rot_bw_from_A(sd(gen$A), real_size = 1000, -1/6, kernel = "epanechnikov")
  cv0 <- env_sim$cv_select_lambda(gen = gen, lambda_grid = LG, m = 50, ch1 = 5, ch2 = 6, h = h1,
           threshold = 0.6, Ker = env_sim$K_epan, n_folds = 5, seed_fold = SEED, cstar_n = "train",
           dt = 1 / (4 * 50), tau = 1)
  s0 <- simulate_cohort(n = 1000, shape = mon, visit_gaps = "unequal", seed = SEED)
  cv1 <- cv_lambda(s0$Y, s0$entry_age, s0$visit_gaps, tau0 = 0, tau1 = 1, bandwidth = h1, seed = SEED)
  report(sprintf("seed %d: fold scores", SEED), maxrel(cv1$scores, cv0$fold_scores), 1e-7)
  report(sprintf("seed %d: mean CV", SEED), maxrel(cv1$cv, cv0$cv_table$cv_score), 1e-7)
  report(sprintf("seed %d: selected lambda", SEED), maxrel(cv1$lambda, cv0$best_lambda), 0)
}

cat("\n=== G. Bootstrap bootstrap_fhat_data_analysis() ===\n")
for (SEED in c(11)) {
  m <- 22; tau0 <- 30; v <- 16 / 12; tau1 <- tau0 + m * v
  s0 <- simulate_cohort(n = 300, shape = mon, visit_gaps = rep(1 / 22, 2), truncation = "informative", seed = SEED)
  A <- tau0 + (tau1 - tau0) * s0$entry_age; Y <- unname(s0$Y)
  h1 <- env_da$rot_bw_from_A(A, -1/6); h2 <- 2.34 * env_da$rot_bw_from_A(A, -1/5)
  b0 <- env_da$bootstrap_fhat_data_analysis(Y, A, m = m, tau0 = tau0, tau1 = tau1, h1 = h1, h2 = h2,
                                            B = 60L, boot_seed = 2026, verbose = FALSE)
  fit <- trajchain(Y, A, visit_gaps = v, tau0 = tau0, m = m)
  b1 <- trajchain_boot(fit, B = 60, seed = 2026)
  map <- c(fhat = "shape", EY = "mean_baseline", bias = "bias_baseline", NW = "nw_baseline")
  for (k in names(map)) {
    report(sprintf("seed %d %s: se", SEED, k), maxrel(b1$se[[map[[k]]]], b0$se[[k]]), 1e-8)
    report(sprintf("seed %d %s: lower/upper", SEED, k),
           max(maxrel(b1$lower[[map[[k]]]], b0$ci_lower[[k]]), maxrel(b1$upper[[map[[k]]]], b0$ci_upper[[k]])), 1e-8)
    report(sprintf("seed %d %s: boot mean, n_ok", SEED, k),
           max(maxrel(b1$boot_mean[[map[[k]]]], b0$boot_mean[[k]]), maxrel(b1$n_ok[[map[[k]]]], b0$n_ok[[k]])), 1e-8)
  }
  report(sprintf("seed %d: n_incomplete", SEED),
         maxrel(b1$n_incomplete[unname(map)], b0$n_incomplete[names(map)]), 0)
  tab0 <- env_da$boot_se_table(b0, ages = seq(30, 55, by = 5), which = c("fhat", "EY", "bias"))
  tab1 <- summary(b1, ages = seq(30, 55, by = 5), which = c("shape", "mean", "bias"))
  report(sprintf("seed %d: boot_se_table vs summary()", SEED), maxrel(tab1[, -1], tab0[, -1]), 1e-8)
  b2 <- trajchain_boot(fit, B = 60, seed = 2026, cores = 2)
  report(sprintf("seed %d: cores = 2 identical to cores = 1", SEED), maxrel(b2$se[, -1], b1$se[, -1]), 0)
}
cat("\nDone.\n")
