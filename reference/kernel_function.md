# Kernel functions with compact support

Returns one of the kernel functions supported by the package. All
kernels are symmetric probability densities supported on \\\[-1, 1\]\\,
as required by the theory in Hu et al. The default, the Epanechnikov
kernel \\K(x) = (3/4)(1 - x^2) I(\|x\| \le 1)\\, is the kernel used
throughout the paper.

## Usage

``` r
kernel_function(kernel = "epanechnikov")
```

## Arguments

- kernel:

  Either the name of a kernel (`"epanechnikov"`, `"biweight"`,
  `"triweight"`, `"triangular"` or `"uniform"`), or a user-supplied
  vectorized function. A user-supplied function must return a
  non-negative value for every element of its argument, preserve the
  dimensions of a matrix argument and vanish outside \\\[-1, 1\]\\.

## Value

A function of one argument `u` (a numeric vector or matrix) that returns
\\K(u)\\ with the same dimensions as `u`.

## Details

Kernel weights enter the estimators only through ratios, so the
normalizing factor \\1/h\\ of \\K_h(\cdot) = K(\cdot/h)/h\\ cancels and
is omitted. The rule-of-thumb constants used by
[`bw_rot()`](https://pingbohu43.github.io/TrajChain/reference/bw_rot.md)
(for example \\2.34\\ for the mean-trajectory bandwidth) were calibrated
for the Epanechnikov kernel.

## Examples

``` r
K <- kernel_function("epanechnikov")
K(c(-1.5, -0.5, 0, 0.5, 1.5))
#> [1] 0.0000 0.5625 0.7500 0.5625 0.0000
integrate(K, -1, 1)$value
#> [1] 1
```
