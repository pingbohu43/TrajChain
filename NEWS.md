# TrajChain 0.1.0

First release, accompanying Hu, Wu, Zhao and Sun, *Nonparametric estimation of
marker trajectories in the presence of death and delayed entry*.

* `trajchain()` estimates the shape function by kernel ratio estimation and
  ratio chaining (equally spaced visits) or penalized ratio matching (widely
  or unequally spaced visits), the unconditional mean trajectory, the
  birth-cohort-survivor bias and the naive Nadaraya-Watson curve, with the
  baseline-only and pooled estimators of the latent mean.
* `cv_bandwidth()` and `cv_lambda()` choose the bandwidth constant and the
  penalty by cross-validation over subjects.
* `trajchain_boot()` computes subject-level bootstrap standard errors and
  pointwise confidence bands, optionally in parallel.
* `plot_trajectories()` and the `plot()` methods draw one or several fits with
  confidence bands.
* `simulate_cohort()`, `sim_truth()`, `shape_monotone()` and `shape_ushape()`
  reproduce the simulation designs of the paper; `simcohort` is a simulated
  example data set.
* Regression tests reproduce the authors' original research code (estimators,
  cross-validation, bootstrap and data generators) to numerical precision.
