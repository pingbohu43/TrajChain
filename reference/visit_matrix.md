# Convert long-format longitudinal data to a visit matrix

Reshapes long-format data (one row per measurement) into the
subject-by-visit matrix used by
[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md),
with one row per subject and one column per scheduled visit. Unobserved
visits (after death, dropout or a missed visit) are filled with `NA`.

## Usage

``` r
visit_matrix(data, id, visit, value)
```

## Arguments

- data:

  A data frame in long format.

- id:

  Name of the column identifying subjects.

- visit:

  Name of the column with the scheduled visit number. Visits must be
  numbered by consecutive integers; the smallest number in the data (for
  example 0 or 1) is taken to be the baseline visit at study entry.

- value:

  Name of the column with the marker values.

## Value

A numeric matrix with one row per subject (row names are the subject
identifiers, in order of first appearance in `data`) and one column per
scheduled visit (column names `"visit_<number>"`).

## See also

[`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md),
[simcohort](https://pingbohu43.github.io/TrajChain/reference/simcohort.md)
for an example data set in long format.

## Examples

``` r
long <- data.frame(
  id = c(1, 1, 1, 2, 2, 3),
  visit = c(0, 1, 2, 0, 1, 0),
  marker = c(1.2, 1.4, 1.5, 0.8, 0.9, 1.1)
)
visit_matrix(long, id = "id", visit = "visit", value = "marker")
#>   visit_0 visit_1 visit_2
#> 1     1.2     1.4     1.5
#> 2     0.8     0.9      NA
#> 3     1.1      NA      NA
```
