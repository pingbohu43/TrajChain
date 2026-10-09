# Summaries of bootstrap results

Methods for objects returned by
[`trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot.md).
[`summary()`](https://rdrr.io/r/base/summary.html) reports the estimate,
bootstrap standard error and pointwise confidence limits of selected
curves at selected ages (values between grid ages are obtained by linear
interpolation of the pointwise quantities; ages outside \\\[\tau_0,
\tau_1\]\\ give `NA`).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
all curves in long format, convenient for plotting with other graphics
packages.

## Usage

``` r
# S3 method for class 'trajchain_boot'
print(x, digits = 4, ...)

# S3 method for class 'trajchain_boot'
summary(
  object,
  ages = NULL,
  which = c("shape", "mean", "bias"),
  mu_method = NULL,
  ...
)

# S3 method for class 'trajchain_boot'
as.data.frame(x, ...)
```

## Arguments

- x, object:

  A `"trajchain_boot"` object.

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- ages:

  Ages at which to report the results (default: six equally spaced ages
  spanning \\\[\tau_0, \tau_1\]\\).

- which:

  Curves to report: any of `"shape"`, `"mean"`, `"bias"`, `"mu"` and
  `"nw"`.

- mu_method:

  `"baseline"` or `"pooled"`; defaults to the estimator chosen when
  fitting.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with columns `quantity`, `age`, `estimate`, `se`, `lower` and `upper`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns a
data frame with columns `curve`, `t`, `estimate`, `se`, `lower`, `upper`
and `n_ok`. [`print()`](https://rdrr.io/r/base/print.html) returns `x`
invisibly.

## Examples

``` r
sim <- simulate_cohort(n = 400, seed = 6)
fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1,
                 n_grid = 101)
bt <- trajchain_boot(fit, B = 30, seed = 2)
summary(bt, ages = c(0.25, 0.5, 0.75), which = "bias")
#>   quantity  age  estimate         se     lower     upper
#> 1     bias 0.25 1.0175193 0.05946288 0.8872623 1.1260209
#> 2     bias 0.50 0.8556707 0.07318980 0.7293524 1.0134644
#> 3     bias 0.75 0.7082850 0.06930421 0.5917326 0.8446858
head(as.data.frame(bt))
#>   curve    t  estimate         se     lower     upper n_ok
#> 1 shape 0.00 0.6100592 0.04407184 0.5383895 0.6918409   30
#> 2 shape 0.01 0.6117020 0.04319075 0.5424971 0.6922013   30
#> 3 shape 0.02 0.6133448 0.04237264 0.5466047 0.6925618   30
#> 4 shape 0.03 0.6149877 0.04162121 0.5506374 0.6929222   30
#> 5 shape 0.04 0.6166305 0.04094015 0.5544305 0.6932826   30
#> 6 shape 0.05 0.6182733 0.04033302 0.5582236 0.6936431   30
```
