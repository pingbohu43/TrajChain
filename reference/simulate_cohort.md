# Simulate longitudinal cohort data with death and delayed entry

Generates data from the simulation designs of Section 6 of Hu et al. on
the unit age scale \\\[0, 1\]\\. In the underlying population, the
latent variable \\Z^0\\ is log-normal, the age at death \\T^0\\ given
\\Z^0 = z\\ is Weibull with shape \\t_1\\ and scale \\t_2 z^{-1/t_1}\\,
and the marker is \\Y^0(t) = Z^0 f(t) e_t\\ with independent log-normal
noise \\e_t\\ of mean one. Entry ages \\A^0\\ follow either

- random truncation: \\A^0 \sim U(a_1, a_2)\\, independent of \\(Z^0,
  T^0)\\ (no birth-cohort effect), or

- informative truncation: \\A^0 + \delta\\ given \\Z^0 = z\\ is Weibull
  with shape \\a_1'\\ and scale \\a_2' z^{-1/a_1'}\\ (birth-cohort
  effect).

## Usage

``` r
simulate_cohort(
  n = 1000,
  shape = "monotone",
  visit_gaps = "equal",
  truncation = c("random", "informative"),
  trunc_par = NULL,
  death_par = c(shape = 4.7, scale = 0.85),
  dropout_par = c(shape = 3, scale = 3),
  frailty_par = c(mean = 1, var = 0.25),
  noise_var = 0.05,
  pop_size = max(10000, 4 * n),
  seed = NULL
)
```

## Arguments

- n:

  Number of enrolled subjects (the sample size).

- shape:

  Shape function: `"monotone"`
  ([`shape_monotone()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)),
  `"ushape"`
  ([`shape_ushape()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)),
  or a vectorized function that is positive and integrates to one on
  \\\[0, 1\]\\.

- visit_gaps:

  `"equal"` for the equally spaced design of the paper (\\k = 5\\
  follow-ups with gap \\v = 1/\lfloor 6 n^{0.16}\rfloor\\), `"unequal"`
  for the unequally spaced design \\(v_1, \dots, v_5) = (5, 5, 6, 5, 6)
  \times 0.02\\, or a numeric vector of gaps (its length is the number
  of follow-up visits \\k\\).

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

- dropout_par:

  Weibull parameters `c(shape = , scale = )` of \\W\\ (default
  `c(3, 3)`; the paper uses scale 5 in its small-sample design).

- frailty_par:

  Mean and variance `c(mean = , var = )` of the log-normal latent
  variable \\Z^0\\ (default `c(1, 0.25)`).

- noise_var:

  Variance of the log-normal noise \\e_t\\, whose mean is fixed at one
  so that \\E\\Y^0(t) \mid Z^0\\ = Z^0 f(t)\\ (default 0.05; use 0 for
  noise-free data).

- pop_size:

  Number of individuals drawn from the underlying population before
  truncation. The default `max(10000, 4 * n)` gives the population size
  used in the paper for \\n \le 2500\\.

- seed:

  Seed, or `NULL` to use the current random number stream. The caller's
  random number stream is left unchanged when a seed is supplied.

## Value

An object of class `"trajchain_sim"`, a list with components

- `Y`:

  \\n \times (k+1)\\ matrix of observed marker values, `NA` after death
  or dropout;

- `entry_age`, `visit_gaps`, `tau0`, `tau1`:

  entry ages, visit gaps and the age interval \\\[0, 1\]\\, ready to be
  passed to
  [`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md);

- `latent`:

  data frame with the latent variable `Z`, the age at death `death_age`
  and the dropout index `dropout` (\\C\\);

- `Y_complete`:

  the marker at every scheduled visit, observed or not;

- `shape_fun`:

  the true shape function;

- `prop_enrolled`, `mean_visits`, `settings`:

  the proportion of the population alive at entry, the mean number of
  observed measurements per subject, and all simulation parameters.

## Details

Individuals are enrolled only if alive at entry (\\A^0 \le T^0\\); the
first `n` enrolled individuals form the sample. Each subject has a
baseline visit at \\A\\ and up to \\k\\ follow-up visits at \\A + v_1 +
\dots + v_l\\; follow-up stops at death or after the \\C\\-th follow-up
visit, where \\C = \min\\k, \lfloor W \rfloor\\\\ and \\W\\ is Weibull,
so that \\P(C = l) = F_W(l + 1) - F_W(l)\\ for \\l \< k\\.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## See also

[`sim_truth()`](https://pingbohu43.github.io/TrajChain/reference/sim_truth.md)
for the true curves of these designs.

## Examples

``` r
sim <- simulate_cohort(n = 1000, seed = 1)
sim
#> <TrajChain simulated cohort>
#>   n = 1000 subjects enrolled (50.8% of the 10000 simulated individuals were alive at entry)
#>   truncation: random; visits: 0.05556, equally spaced (5 follow-up visits)
#>   mean number of observed measurements per subject: 2.91
#>   pass Y, entry_age and visit_gaps to trajchain() with tau0 = 0 and tau1 = 1
head(sim$Y)
#>        visit_0   visit_1   visit_2   visit_3 visit_4 visit_5
#> [1,] 0.4900233 0.3880605 0.6475243        NA      NA      NA
#> [2,] 0.6541136 0.6262784 0.5513313        NA      NA      NA
#> [3,] 0.3166980 0.3595320 0.2217208        NA      NA      NA
#> [4,] 1.2161601 1.6517830 2.2808034        NA      NA      NA
#> [5,] 0.5669925 0.6044165 0.7037160 0.5345847      NA      NA
#> [6,] 0.6852364 0.5777045 0.8769046        NA      NA      NA

# unequally spaced visits with informative truncation
sim2 <- simulate_cohort(n = 1000, shape = "ushape", visit_gaps = "unequal",
                        truncation = "informative", seed = 2)
sim2$visit_gaps
#> [1] 0.10 0.10 0.12 0.10 0.12
```
