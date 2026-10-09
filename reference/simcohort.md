# Simulated cohort with death and delayed entry on the age scale

A simulated longitudinal cohort in long format, used in the examples and
the vignette. Subjects enter the study at different ages, have a
baseline visit and up to two follow-up visits scheduled every 16 months,
and are lost to follow-up through death or dropout. The data were
generated with
[`simulate_cohort()`](https://pingbohu43.github.io/TrajChain/reference/simulate_cohort.md)
under the informative-truncation design of the paper (so that both
survivor and birth-cohort selection are present), with the unit age
scale mapped to ages 30 to 59.33, i.e. a grid of \\m = 22\\ points with
step 16 months, as in the data analysis of the paper. The marker in
group A follows the monotone shape function
[`shape_monotone()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md)
and in group B the U-shaped function
[`shape_ushape()`](https://pingbohu43.github.io/TrajChain/reference/shape_functions.md),
rescaled to the age interval. These are simulated data, not observations
on real individuals.

## Usage

``` r
simcohort
```

## Format

A data frame with 1502 rows (one per observed measurement) and 6
variables:

- `id`:

  subject identifier (1 to 600);

- `group`:

  factor with levels `A` and `B` (300 subjects each);

- `visit`:

  scheduled visit number: 0 (baseline), 1 or 2;

- `age`:

  age at the visit, in years;

- `entry_age`:

  age at study entry, in years;

- `marker`:

  the (positive) marker value.

## Source

Simulated; see `data-raw/simcohort.R` in the package sources.

## Examples

``` r
head(simcohort)
#>   id group visit     age entry_age  marker
#> 1  1     A     0 37.1463   37.1463 0.53748
#> 2  1     A     1 38.4796   37.1463 0.68528
#> 3  1     A     2 39.8129   37.1463 0.88437
#> 4  2     A     0 47.7841   47.7841 1.19542
#> 5  2     A     1 49.1174   47.7841 1.00355
#> 6  3     A     0 56.0942   56.0942 0.90858
groupA <- subset(simcohort, group == "A")
Y <- visit_matrix(groupA, id = "id", visit = "visit", value = "marker")
A <- groupA$entry_age[match(rownames(Y), groupA$id)]
fit <- trajchain(Y, A, visit_gaps = 16 / 12, tau0 = 30, m = 22)
fit
#> <TrajChain fit>
#>   Shape estimator : ratio chaining (Section 3)
#>   Mean estimator  : baseline measurements (equation 5)
#>   Subjects        : 300, with 762 measurements (2.54 per subject)
#>   Visit gaps      : 1.333, equally spaced (2 follow-up visits)
#>   Age interval    : [30, 59.33], m = 22 grid points (step 1.333)
#>   Bandwidths      : h = 2.832 (ratio), h' = 5.48 (mean)
#> 
#>    age   shape   mean   bias     nw
#>  30.00 0.02102 0.7559 1.0000 0.7559
#>  35.87 0.02016 0.7250 1.2547 0.9097
#>  41.73 0.02693 0.9685 1.0310 0.9985
#>  47.60 0.03929 1.4130 0.8901 1.2578
#>  53.47 0.04632 1.6659 0.8329 1.3874
#>  59.33 0.05049 1.8157 0.5123 0.9301
```
