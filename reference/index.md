# Package index

## Estimation

Fit the shape function, mean trajectories and birth-cohort-survivor
bias.

- [`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md)
  : Estimate marker trajectories in the presence of death and delayed
  entry
- [`print(`*`<trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain-methods.md)
  [`summary(`*`<trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain-methods.md)
  [`print(`*`<summary.trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain-methods.md)
  [`predict(`*`<trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain-methods.md)
  [`as.data.frame(`*`<trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain-methods.md)
  : Print, summarize and extract a trajchain fit
- [`ratio_estimate()`](https://pingbohu43.github.io/TrajChain/reference/ratio_estimate.md)
  : Kernel estimator of the ratio of the shape function at two ages

## Tuning

Bandwidths, kernels and cross-validation.

- [`bw_rot()`](https://pingbohu43.github.io/TrajChain/reference/bw_rot.md)
  : Rule-of-thumb bandwidths based on the spread of entry ages
- [`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)
  : Cross-validated choice of the bandwidth constant
- [`cv_lambda()`](https://pingbohu43.github.io/TrajChain/reference/cv_lambda.md)
  : Cross-validated choice of the penalty of the penalized estimator
- [`kernel_function()`](https://pingbohu43.github.io/TrajChain/reference/kernel_function.md)
  : Kernel functions with compact support

## Inference

- [`trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot.md)
  : Subject-level bootstrap standard errors and confidence bands
- [`print(`*`<trajchain_boot>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot-methods.md)
  [`summary(`*`<trajchain_boot>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot-methods.md)
  [`as.data.frame(`*`<trajchain_boot>`*`)`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot-methods.md)
  : Summaries of bootstrap results

## Graphics

- [`plot_trajectories()`](https://pingbohu43.github.io/TrajChain/reference/plot_trajectories.md)
  : Plot estimated trajectories for one or several groups
- [`plot(`*`<trajchain>`*`)`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md)
  [`plot(`*`<trajchain_boot>`*`)`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md)
  : Plot a trajchain fit or its bootstrap confidence bands
- [`plot(`*`<trajchain_cv>`*`)`](https://pingbohu43.github.io/TrajChain/reference/plot.trajchain_cv.md)
  : Plot a cross-validation criterion

## Data and simulation

- [`visit_matrix()`](https://pingbohu43.github.io/TrajChain/reference/visit_matrix.md)
  : Convert long-format longitudinal data to a visit matrix
- [`simcohort`](https://pingbohu43.github.io/TrajChain/reference/simcohort.md)
  : Simulated cohort with death and delayed entry on the age scale
- [`simulate_cohort()`](https://pingbohu43.github.io/TrajChain/reference/simulate_cohort.md)
  : Simulate longitudinal cohort data with death and delayed entry
- [`sim_truth()`](https://pingbohu43.github.io/TrajChain/reference/sim_truth.md)
  : True curves of the simulation designs
- [`shape_monotone()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)
  [`shape_ushape()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)
  : Shape functions used in the simulation studies
