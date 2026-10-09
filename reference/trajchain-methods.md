# Print, summarize and extract a trajchain fit

Methods for objects returned by
[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md).
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) describe the design
and report the estimates at a few ages;
[`predict()`](https://rdrr.io/r/stats/predict.html) evaluates the
estimated curves at arbitrary ages;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the curves on the dense age grid.

## Usage

``` r
# S3 method for class 'trajchain'
print(x, digits = 4, ...)

# S3 method for class 'trajchain'
summary(object, ages = NULL, mu_method = NULL, ...)

# S3 method for class 'summary.trajchain'
print(x, digits = 4, ...)

# S3 method for class 'trajchain'
predict(
  object,
  newdata,
  type = c("shape", "mean", "bias", "mu", "nw"),
  mu_method = NULL,
  ...
)

# S3 method for class 'trajchain'
as.data.frame(x, ...)
```

## Arguments

- x, object:

  A `"trajchain"` object.

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- ages:

  Ages at which to report the estimates (default: six equally spaced
  ages spanning \\\[\tau_0, \tau_1\]\\).

- mu_method:

  `"baseline"` or `"pooled"`: which estimator of \\\mu\_{Z^0}\\ to
  report. Defaults to the one chosen in
  [`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md).

- newdata:

  Numeric vector of ages at which to evaluate the curve.

- type:

  The curve to evaluate: `"shape"` (\\\hat f\\), `"mean"` (unconditional
  mean), `"bias"` (birth-cohort-survivor bias), `"mu"`
  (\\\hat\mu\_{Z^0}\\) or `"nw"` (Nadaraya–Watson curve).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.trajchain"` whose `table` component is a data frame of
estimates at `ages`. [`predict()`](https://rdrr.io/r/stats/predict.html)
returns a numeric vector (`NA` outside \\\[\tau_0, \tau_1\]\\).
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
`x$curves`.

## Details

[`predict()`](https://rdrr.io/r/stats/predict.html) recomputes the
kernel sums at the requested ages when the fit contains the data
(`keep_data = TRUE`), so its values agree exactly with `x$curves` at the
grid ages; otherwise it interpolates `x$curves` linearly.

## Examples

``` r
sim <- simulate_cohort(n = 500, seed = 5)
fit <- trajchain(sim$Y, sim$entry_age, sim$visit_gaps, tau0 = 0, tau1 = 1)
summary(fit, ages = c(0.1, 0.5, 0.9))
#> <Summary of TrajChain fit>
#>   Shape estimator : ratio chaining (Section 3)
#>   Mean estimator  : baseline measurements (equation 5)
#>   Subjects        : 500, with 1452 measurements (2.90 per subject)
#>   Visit gaps      : 0.0625, equally spaced (5 follow-up visits)
#>   Age interval    : [0, 1], m = 16 grid points (step 0.0625)
#>   Bandwidths      : h = 0.1289 (ratio), h' = 0.2451 (mean), kernel: epanechnikov
#>   mu(tau0)        : 1.003 (baseline), 1.024 (pooled)
#> 
#> Estimates (shape = f(t); mean = unconditional mean; bias = birth-cohort-survivor bias;
#> nw = Nadaraya-Watson smoothing of the observed markers):
#>  age  shape   mean   bias     nw
#>  0.1 0.6351 0.6369 1.0808 0.6884
#>  0.5 1.0727 1.0759 0.8799 0.9467
#>  0.9 1.2620 1.2657 0.8696 1.1007
predict(fit, c(0.25, 0.75), type = "shape")
#> [1] 0.7024864 1.2677420
predict(fit, c(0.25, 0.75), type = "bias", mu_method = "pooled")
#> [1] 1.0573324 0.8438136
head(as.data.frame(fit))
#>         t     shape mean_baseline bias_baseline mu_baseline nw_baseline
#> 1 0.00000 0.6708630     0.6728228      1.000000    1.002921   0.6728228
#> 2 0.00125 0.6703743     0.6723326      1.000699    1.003623   0.6728029
#> 3 0.00250 0.6698856     0.6718425      1.001471    1.004397   0.6728311
#> 4 0.00375 0.6693969     0.6713524      1.002232    1.005160   0.6728507
#> 5 0.00500 0.6689082     0.6708623      1.003015    1.005946   0.6728852
#> 6 0.00625 0.6684195     0.6703721      1.003770    1.006702   0.6728993
#>   mean_pooled bias_pooled mu_pooled nw_pooled
#> 1   0.6868995    1.000000  1.023904 0.6868995
#> 2   0.6863992    1.000845  1.024769 0.6869789
#> 3   0.6858988    1.001702  1.025647 0.6870661
#> 4   0.6853984    1.002553  1.026519 0.6871485
#> 5   0.6848980    1.003411  1.027397 0.6872340
#> 6   0.6843976    1.004259  1.028265 0.6873123
```
