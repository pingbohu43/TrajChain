# Kernel estimator of the ratio of the shape function at two ages

Computes the kernel ratio estimator of equation (4) of Hu et al.,
\$\$\hat r(t, t - v) = \frac{\sum\_{i=1}^n \int K_h(t-u)\\ Y_i(u)\\
dN_i(u)} {\sum\_{i=1}^n \int K_h(t-u)\\ Y_i(u - v)\\ dN_i(u)},\$\$ which
estimates \\r(t, t - v) = f(t)/f(t - v)\\, the relative change of the
shape function between two ages \\v\\ apart. The counting process
\\N_i\\ jumps at the follow-up visits of subject \\i\\ for which both
the current measurement and the measurement at the preceding visit
(scheduled \\v\\ earlier) are observed.

## Usage

``` r
ratio_estimate(
  Y,
  entry_age,
  visit_gaps,
  t,
  gap = NULL,
  bandwidth = "rot",
  kernel = "epanechnikov"
)
```

## Arguments

- Y:

  Numeric matrix of marker values with one row per subject and one
  column per scheduled visit; column 1 is the baseline measurement at
  study entry. Visits that were not observed (after death, dropout or
  end of follow-up) are `NA`.

- entry_age:

  Numeric vector of ages at study entry (\\A_i\\), one per row of `Y`.

- visit_gaps:

  Scheduled time between consecutive visits: a single number for equally
  spaced visits, or a vector of length `ncol(Y) - 1` giving \\v_1,
  \dots, v_k\\. Visit \\l\\ of subject \\i\\ takes place at age \\A_i +
  v_1 + \dots + v_l\\.

- t:

  Ages at which to evaluate the ratio.

- gap:

  The age difference \\v\\ of the ratio. It must equal one of the gaps
  in `visit_gaps`; it can be omitted when all gaps are equal.

- bandwidth:

  Bandwidth \\h\\ on the age scale, or `"rot"` for the rule-of-thumb
  \\\hat\sigma_A n^{-1/6}\\.

- kernel:

  Kernel; see
  [`kernel_function()`](https://pingbohu43.github.io/TrajChain/reference/kernel_function.md).

## Value

A data frame with one row per element of `t` and columns

- `t`:

  evaluation age;

- `ratio`:

  \\\hat r(t, t - v)\\ (`NA` when no observed pair falls in the kernel
  window);

- `numerator`, `denominator`:

  the two kernel sums;

- `n_pairs`:

  number of observed visit pairs whose current visit lies within one
  bandwidth of `t`.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## See also

[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md),
which chains these ratios into the shape function.

## Examples

``` r
sim <- simulate_cohort(n = 500, seed = 1)
ratio_estimate(sim$Y, sim$entry_age, sim$visit_gaps, t = c(0.3, 0.5, 0.7))
#>     t    ratio numerator denominator n_pairs
#> 1 0.3 1.067501   72.8492    68.24273     199
#> 2 0.5 1.115323  126.2191   113.16822     236
#> 3 0.7 1.069665  112.5800   105.24796     191

# true ratios f(t) / f(t - v)
v <- sim$visit_gaps[1]
shape_monotone(c(0.3, 0.5, 0.7)) / shape_monotone(c(0.3, 0.5, 0.7) - v)
#> [1] 1.054805 1.123692 1.049373
```
