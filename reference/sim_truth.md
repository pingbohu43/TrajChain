# True curves of the simulation designs

Computes, by numerical integration over the log-normal latent variable,
the true shape function, \\\mu\_{Z^0}(t) = E(Z^0 \mid A^0 = t, T^0 \ge
t)\\, the birth-cohort-survivor bias \\\mu\_{Z^0}(t)/\mu\_{Z^0}(0)\\ and
the unconditional mean \\\mu\_{Z^0}(0) f(t)\\ for the designs of
[`simulate_cohort()`](https://pingbohu43.github.io/TrajChain/reference/simulate_cohort.md).
Under random truncation \\\mu\_{Z^0}(t) = E(Z^0 \mid T^0 \ge t)\\, the
pure survivor bias, and \\\mu\_{Z^0}(0) = E(Z^0)\\.

## Usage

``` r
sim_truth(
  t,
  shape = "monotone",
  truncation = c("random", "informative"),
  trunc_par = NULL,
  death_par = c(shape = 4.7, scale = 0.85),
  frailty_par = c(mean = 1, var = 0.25)
)
```

## Arguments

- t:

  Ages on the unit scale at which to evaluate the curves.

- shape:

  Shape function: `"monotone"`
  ([`shape_monotone()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)),
  `"ushape"`
  ([`shape_ushape()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)),
  or a vectorized function that is positive and integrates to one on
  \\\[0, 1\]\\.

- truncation:

  `"random"` or `"informative"`.

- trunc_par:

  Parameters of the entry-age distribution: for random truncation a
  named vector `c(min = , max = )` (default `c(-0.3, 1.9)`); for
  informative truncation `c(shape = , scale = , shift = )` (default
  `c(3.2, 1.2, 0.3)`).

- death_par:

  Weibull parameters `c(shape = , scale = )` of the age at death
  (default `c(4.7, 0.85)`).

- frailty_par:

  Mean and variance `c(mean = , var = )` of the log-normal latent
  variable \\Z^0\\ (default `c(1, 0.25)`).

## Value

A data frame with columns `t`, `shape`, `mu`, `bias` and `mean`.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## Examples

``` r
truth <- sim_truth(seq(0, 1, by = 0.1))
truth
#>      t     shape        mu      bias      mean
#> 1  0.0 0.6412312 1.0000000 1.0000000 0.6412312
#> 2  0.1 0.6494445 0.9999893 0.9999893 0.6494445
#> 3  0.2 0.6708552 0.9997218 0.9997218 0.6708552
#> 4  0.3 0.7230567 0.9981344 0.9981344 0.7230567
#> 5  0.4 0.8319574 0.9928502 0.9928502 0.8319574
#> 6  0.5 1.0000000 0.9800188 0.9800188 1.0000000
#> 7  0.6 1.1680426 0.9548596 0.9548596 1.1680426
#> 8  0.7 1.2769433 0.9132523 0.9132523 1.2769433
#> 9  0.8 1.3291448 0.8539093 0.8539093 1.3291448
#> 10 0.9 1.3505555 0.7795938 0.7795938 1.3505555
#> 11 1.0 1.3587688 0.6962223 0.6962223 1.3587688

# informative truncation: birth-cohort and survivor effects combined
sim_truth(c(0, 0.5, 1), truncation = "informative")
#>     t     shape        mu      bias      mean
#> 1 0.0 0.6412312 1.2454020 1.0000000 0.7985906
#> 2 0.5 1.0000000 1.1315579 0.9085885 1.2454020
#> 3 1.0 1.3587688 0.7014274 0.5632136 1.6922134
```
