# Shape functions used in the simulation studies

The two shape functions of Section 6 of Hu et al., both probability
densities on \\\[0, 1\]\\:

- `shape_monotone()`: the monotonically increasing function \\f(t) =
  8/\\11 + 11\exp(5 - 10t)\\ + 7/11\\;

- `shape_ushape()`: the U-shaped function \\f(t) = 4/5 + (60/37)(t -
  3/10)^2\\.

## Usage

``` r
shape_monotone(t)

shape_ushape(t)
```

## Arguments

- t:

  Numeric vector of ages on the unit scale.

## Value

Numeric vector of \\f(t)\\.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## Examples

``` r
curve(shape_monotone(x), 0, 1, ylim = c(0.5, 1.6), ylab = "f(t)", xlab = "t")
curve(shape_ushape(x), 0, 1, add = TRUE, lty = 2)

integrate(shape_monotone, 0, 1)$value
#> [1] 1
integrate(shape_ushape, 0, 1)$value
#> [1] 1
```
