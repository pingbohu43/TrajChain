# Plot a cross-validation criterion

Plots the cross-validation criterion returned by
[`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)
or
[`cv_lambda()`](https://pingbohu43.github.io/TrajChain/reference/cv_lambda.md)
against the candidate values (logarithmic axis), together with the
criterion for individual partitions or folds (gray) and the selected
value (dashed line). A U-shaped curve with an interior minimum indicates
that the candidate set covers the optimum.

## Usage

``` r
# S3 method for class 'trajchain_cv'
plot(x, n_curves = 30, ...)
```

## Arguments

- x:

  A `"trajchain_cv"` object.

- n_curves:

  Maximum number of individual partitions or folds to draw.

- ...:

  Further graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`NULL`, invisibly.

## Examples

``` r
sim <- simulate_cohort(n = 500, seed = 9)
cvb <- cv_bandwidth(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
plot(cvb)
```
