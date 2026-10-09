# Creates tests/testthat/fixtures/original-code.rds: reference values computed
# with the authors' original research scripts (not distributed with the
# package), against which the regression tests in
# tests/testthat/test-regression.R compare TrajChain.
#
# Set ORIG to a folder containing copies of
#   assumption-B-functions.R          (independent truncation)  -> indep_functions.R
#   functions.R                        (dependent truncation)    -> dep_functions.R
#   assumption-B-cv-lambda.R                                     -> cv_lambda.R
#   data-analysis-Case-1-including-birth-cohort-CV-bandwidth.R   -> da_cv.R
#   data-analysis-Case-1-including-birth-cohort-bootstrap.R      -> da_boot.R
#   assumption-B-Case-1-Mon-CV-small-bandwidth.R                 -> sim_case1_cv.R
#   assumption-B-Case-3-Mon-CV-both-bandwidth-lambda.R           -> sim_case3_cv.R
# with the library(), source() and FUN_PATH lines removed.
ORIG <- Sys.getenv("TRAJCHAIN_ORIG")
stopifnot(nzchar(ORIG))
suppressMessages(library(Matrix))

load_defs <- function(file, env) {
  for (e in parse(file)) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e, env)
  }
  invisible(env)
}
env_sim <- new.env()
sys.source(file.path(ORIG, "indep_functions.R"), envir = env_sim)
sys.source(file.path(ORIG, "cv_lambda.R"), envir = env_sim)
env_dep <- new.env()
sys.source(file.path(ORIG, "dep_functions.R"), envir = env_dep)
env_da <- new.env()
load_defs(file.path(ORIG, "da_cv.R"), env_da)
load_defs(file.path(ORIG, "da_boot.R"), env_da)
env_da$n_dense <- 801
env_c1 <- new.env(parent = env_sim)
load_defs(file.path(ORIG, "sim_case1_cv.R"), env_c1)
env_c3 <- new.env(parent = env_sim)
load_defs(file.path(ORIG, "sim_case3_cv.R"), env_c3)

mon <- env_sim$mon_shape
ush <- env_sim$U_quad_pdf
K <- env_sim$K_epan
idx <- seq(1, 801, by = 40)
cand <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4)
mask <- function(Y, ind) {
  Y[!ind] <- NA
  Y
}
fx <- list()

## 1) Case 1, simulation pipeline (random truncation, equal spacing), n = 300
rs <- 300
m <- floor(6 * rs^0.16)
v <- 1 / m
SEED <- 21
gen <- env_sim$generate_data(SEED, 10000, rs, mon, 5, v, 0, 1, 0.25, -0.3, 1.9, 4.7, 0.85,
                             1, 0.05, 3, 3)
ind <- outer(gen$A, (0:5) * v, "+") <= matrix(gen$T_tilde, rs, 6)
h1 <- env_sim$rot_bw_from_A(sd(gen$A), rs, -1 / 6)
h2 <- env_sim$rot_bw_from_A(sd(gen$A), rs, -1 / 5)
out <- env_sim$compute_fhat(SEED = SEED, m = m, h1 = h1, h2 = 2.34 * h2, n = 10000,
                            real_size = rs, f = mon, k = 5, delta = 0, z1 = 1, z2 = 0.25,
                            a1 = -0.3, a2 = 1.9, t1 = 4.7, t2 = 0.85, e1 = 1, e2 = 0.05,
                            c1 = 3, c2 = 3, Ker = K, return_gen = FALSE, n_dense = 801)
cv1 <- env_c1$cv_select_c(gen, m = m, cand_c = cand, n_folds = 5L, Ker = K, fold_seed = 99L)
fx$case1 <- list(seed = SEED, n = rs, Y = mask(gen$Y, ind), A = gen$A, v = v, m = m,
                 h1 = h1, h2 = 2.34 * h2, fhat_grid = out$results$fhat,
                 rhat = out$results$rhat[-1], idx = idx, fhat = out$dense$fhat[idx],
                 g = out$dense$g[idx], EY = out$dense$EY[idx], EZ0 = out$EZ0,
                 cv = cv1$cv, cv_folds = cv1$cv_folds, c_opt = cv1$c_opt)

## 2) Case 3 (unequal spacing), n = 400
rs <- 400
SEED <- 22
gen <- env_sim$generate_data_different_visit_interval(SEED = SEED, n = 10000, real_size = rs,
         f = ush, k = 5, v = 0.02, delta = 0, ch1 = 5, ch2 = 6, z1 = 1, z2 = 0.25, a1 = -0.3,
         a2 = 1.9, t1 = 4.7, t2 = 0.85, e1 = 1, e2 = 0.05, c1 = 3, c2 = 3)
ind <- gen$visit_age <= matrix(gen$T_tilde, rs, 6)
h1 <- env_sim$rot_bw_from_A(sd(gen$A), rs, -1 / 6)
h2 <- env_sim$rot_bw_from_A(sd(gen$A), rs, -1 / 5)
out <- env_sim$optimize_fhat_with_pairs_cvxr(m = 50, w1 = 0, threshold = 0.6, SEED = SEED,
         h1 = h1, h2 = 2.34 * h2, Ker = K, n = 10000, real_size = rs, f = ush, k = 5, delta = 0,
         z1 = 1, z2 = 0.25, a1 = -0.3, a2 = 1.9, t1 = 4.7, t2 = 0.85, e1 = 1, e2 = 0.05,
         c1 = 3, c2 = 3, n_dense = 801, ridge = 3e-5, ch1 = 5, ch2 = 6)
lg <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
cvl <- env_sim$cv_select_lambda(gen = gen, lambda_grid = lg, m = 50, ch1 = 5, ch2 = 6, h = h1,
         threshold = 0.6, Ker = K, n_folds = 5, seed_fold = 5, cstar_n = "train",
         dt = 1 / 200, tau = 1)
cv3 <- env_c3$cv_select_c(gen, m = 50, cand_c = cand, ch1 = 5, ch2 = 6, n_folds = 5L, Ker = K,
                          fold_seed = 77L, tau = 1)
fx$case3 <- list(seed = SEED, n = rs, Y = mask(gen$Y, ind), A = gen$A,
                 gaps = c(5, 5, 6, 5, 6) * 0.02, h1 = h1, h2 = 2.34 * h2, lambda = 3e-5,
                 counts = c(out$number_5_vec, out$number_6_vec), fhat_grid = out$results$fhat,
                 idx = idx, fhat = out$dense$fhat[idx], g = out$dense$g[idx],
                 EY = out$dense$EY[idx], lambda_grid = lg, cv_scores = cvl$fold_scores,
                 best_lambda = cvl$best_lambda, cv_c = cv3$cv, c_opt = cv3$c_opt)

## 3) Data-analysis pipeline on the age scale (informative truncation), n = 300
rs <- 300
SEED <- 23
v1 <- 1 / 22
gen <- env_dep$generate_data(SEED, 10000, rs, mon, 2, v1, 0.3, 1, 0.25, 3.2, 1.2, 4.7, 0.85,
                             1, 0.05, 3, 3)
ind <- outer(gen$A, (0:2) * v1, "+") <= matrix(gen$T_tilde, rs, 3)
tau0 <- 30
m <- 22
tau1 <- tau0 + m * 16 / 12
A <- tau0 + (tau1 - tau0) * gen$A
Y <- mask(gen$Y, ind)
h1 <- env_da$rot_bw_from_A(A, -1 / 6)
h2 <- 2.34 * env_da$rot_bw_from_A(A, -1 / 5)
res <- env_da$compute_fhat_data_analysis(Y, h1, h2, Ker = env_da$K_epan, m = m, A = A,
                                         tau0 = tau0, tau1 = tau1)
cvd <- env_da$cv_select_c_data_analysis(Y, A, m, tau0, tau1, cand, n_folds = 5,
                                        Ker = env_da$K_epan, fold_seed = 2026, n_repeats = 2)
bt <- env_da$bootstrap_fhat_data_analysis(Y, A, m = m, tau0 = tau0, tau1 = tau1, h1 = h1,
                                          h2 = h2, B = 30L, boot_seed = 7, verbose = FALSE)
sb <- res$survival_birth_cohort_result
nb <- res$result_no_birth_cohort
fx$data_analysis <- list(
  seed = SEED, Y = Y, A = A, tau0 = tau0, tau1 = tau1, m = m, h1 = h1, h2 = h2, idx = idx,
  fhat = sb$fhat[idx], bias_baseline = sb$survival_birth_cohort_bias[idx],
  mean_baseline = sb$EY_survival_birth_cohort[idx], nw_baseline = sb$NW_survival_birth_cohort[idx],
  bias_pooled = nb$surv_bias[idx], mean_pooled = nb$EY[idx], nw_pooled = nb$NW_no_birth_cohort[idx],
  EZ0 = c(baseline = res$EZ0_survival_birth_cohort, pooled = res$EZ0),
  cv = cvd$cv, cv_by_repeat = cvd$cv_by_repeat, c_opt = cvd$c_opt,
  boot = list(B = 30, seed = 7, se = bt$se[idx, ], lower = bt$ci_lower[idx, ],
              upper = bt$ci_upper[idx, ], n_incomplete = bt$n_incomplete))

## 4) Generators: informative truncation with unequal spacing
gen <- env_dep$generate_data_different_visit_interval(24, 10000, 300, ush, 5, 0.02, 0.3, 5, 6, 1,
         0.25, 3.2, 1.2, 4.7, 0.85, 1, 0.05, 3, 3)
ind <- gen$visit_age <= matrix(gen$T_tilde, 300, 6)
fx$gen_informative_unequal <- list(seed = 24, Y = mask(gen$Y, ind), A = gen$A, C = gen$C)

saveRDS(fx, "tests/testthat/fixtures/original-code.rds", version = 2)
cat("fixture size:", file.size("tests/testthat/fixtures/original-code.rds"), "bytes\n")
