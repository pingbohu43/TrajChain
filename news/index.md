# Changelog

## TrajChain 0.1.0

First release, accompanying Hu, Wu, Zhao and Sun, *Nonparametric
estimation of marker trajectories in the presence of death and delayed
entry*.

- [`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md)
  estimates the shape function by kernel ratio estimation and ratio
  chaining (equally spaced visits) or penalized ratio matching (widely
  or unequally spaced visits), the unconditional mean trajectory, the
  birth-cohort-survivor bias and the naive Nadaraya-Watson curve, with
  the baseline-only and pooled estimators of the latent mean.
- [`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)
  and
  [`cv_lambda()`](https://pingbohu43.github.io/TrajChain/reference/cv_lambda.md)
  choose the bandwidth constant and the penalty by cross-validation over
  subjects.
- [`trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot.md)
  computes subject-level bootstrap standard errors and pointwise
  confidence bands, optionally in parallel.
- [`plot_trajectories()`](https://pingbohu43.github.io/TrajChain/reference/plot_trajectories.md)
  and the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods draw one or several fits with confidence bands.
- [`simulate_cohort()`](https://pingbohu43.github.io/TrajChain/reference/simulate_cohort.md),
  [`sim_truth()`](https://pingbohu43.github.io/TrajChain/reference/sim_truth.md),
  [`shape_monotone()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)
  and
  [`shape_ushape()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)
  reproduce the simulation designs of the paper; `simcohort` is a
  simulated example data set.
- Regression tests reproduce the authors’ original research code
  (estimators, cross-validation, bootstrap and data generators) to
  numerical precision.
