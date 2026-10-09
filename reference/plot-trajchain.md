# Plot a trajchain fit or its bootstrap confidence bands

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md)
and
[`trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot.md)
results. They call
[`plot_trajectories()`](https://pingbohu43.github.io/TrajChain/reference/plot_trajectories.md)
for a single group; bootstrap results add pointwise confidence bands.

## Usage

``` r
# S3 method for class 'trajchain'
plot(x, which = c("shape", "mean", "nw", "bias"), mu_method = NULL, ...)

# S3 method for class 'trajchain_boot'
plot(
  x,
  which = c("shape", "mean", "nw", "bias"),
  mu_method = NULL,
  ci = c("percentile", "normal"),
  ...
)
```

## Arguments

- x:

  A `"trajchain"` or `"trajchain_boot"` object.

- which:

  Panels to draw, in order: any of `"shape"` (shape function), `"mean"`
  (unconditional mean), `"nw"` (simple kernel smoothing of the observed
  markers), `"bias"` (birth-cohort-survivor bias) and `"mu"`
  (\\\hat\mu\_{Z^0}\\).

- mu_method:

  `"baseline"` or `"pooled"`; defaults to the estimator chosen when
  fitting the first object.

- ...:

  Further arguments passed to
  [`plot_trajectories()`](https://pingbohu43.github.io/TrajChain/reference/plot_trajectories.md).

- ci:

  Type of band: pointwise `"percentile"` intervals (default) or
  `"normal"` intervals, estimate \\\pm z\_{1-\alpha/2}\\ times the
  bootstrap standard error.

## Value

`NULL`, invisibly.

## Examples

``` r
sim <- simulate_cohort(n = 500, seed = 8)
fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
                 mu_method = "pooled")
plot(fit)

plot(fit, which = c("shape", "bias"), truth = sim_truth(fit$curves$t))
```
