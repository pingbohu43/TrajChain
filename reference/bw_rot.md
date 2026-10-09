# Rule-of-thumb bandwidths based on the spread of entry ages

Computes \\h = c\\\hat\sigma_A\\ n^{\alpha}\\, where \\\hat\sigma_A\\ is
the sample standard deviation of the entry ages and \\n\\ is the number
of subjects. This is the form of all bandwidths used in the paper.

## Usage

``` r
bw_rot(entry_age, exponent = -1/6, constant = 1)
```

## Arguments

- entry_age:

  Numeric vector of entry ages \\A_i\\ (one per subject). Non-finite
  values are ignored.

- exponent:

  The exponent \\\alpha\\. The default \\-1/6\\ gives the "small"
  bandwidth for the ratio estimator (\\mh \to 0\\); the paper's "large"
  bandwidth uses \\-0.14\\ and the mean-trajectory bandwidth uses
  \\-1/5\\.

- constant:

  The multiplicative constant \\c\\ (\\\kappa\\ in the paper). Use
  `2.34` together with `exponent = -1/5` for the bandwidth \\h'\\ of the
  mean-trajectory estimators (the rule-of-thumb constant of the
  Epanechnikov kernel), and a value selected by
  [`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)
  for the ratio estimator.

## Value

A single positive number, on the same scale as `entry_age`.

## Details

[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md)
uses, by default, the ratio bandwidth \\h = \hat\sigma_A n^{-1/6}\\ and
the mean-trajectory bandwidth \\h' = 2.34\\\hat\sigma_A n^{-1/5}\\, as
in the authors' simulation and data-analysis code.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## See also

[`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)
to choose the constant by cross-validation.

## Examples

``` r
set.seed(1)
age <- runif(200, 25, 80)
bw_rot(age)                                   # ratio estimator
#> [1] 6.1187
bw_rot(age, exponent = -1/5, constant = 2.34) # mean trajectory
#> [1] 11.9998
```
