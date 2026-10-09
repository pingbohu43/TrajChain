# Plot estimated trajectories for one or several groups

Draws side-by-side panels of the estimated shape function, unconditional
mean, Nadaraya–Watson curve and birth-cohort-survivor bias (the layout
of Fig. 5 of Hu et al.), for one fit or for several groups (for example,
women and men analyzed separately). Bootstrap results add pointwise
confidence bands.

## Usage

``` r
plot_trajectories(
  x,
  which = c("shape", "mean", "nw", "bias"),
  mu_method = NULL,
  band = TRUE,
  ci = c("percentile", "normal"),
  truth = NULL,
  col = NULL,
  lty = NULL,
  lwd = 2,
  xlab = "Age",
  ylab = NULL,
  legend = NULL,
  include_zero = TRUE,
  mfrow = NULL,
  ...
)
```

## Arguments

- x:

  A `"trajchain"` or `"trajchain_boot"` object, or a named list of such
  objects (one per group; the names are used in the legend).

- which:

  Panels to draw, in order: any of `"shape"` (shape function), `"mean"`
  (unconditional mean), `"nw"` (simple kernel smoothing of the observed
  markers), `"bias"` (birth-cohort-survivor bias) and `"mu"`
  (\\\hat\mu\_{Z^0}\\).

- mu_method:

  `"baseline"` or `"pooled"`; defaults to the estimator chosen when
  fitting the first object.

- band:

  Draw confidence bands for bootstrap objects.

- ci:

  Type of band: pointwise `"percentile"` intervals (default) or
  `"normal"` intervals, estimate \\\pm z\_{1-\alpha/2}\\ times the
  bootstrap standard error.

- truth:

  Optional data frame of true curves, as returned by
  [`sim_truth()`](https://pingbohu43.github.io/TrajChain/reference/sim_truth.md),
  drawn as a dashed black line (useful for simulations; its `t` must be
  on the age scale of the fits).

- col, lty, lwd:

  Colors, line types and line width, one per group.

- xlab:

  Label of the horizontal axis.

- ylab:

  Labels of the vertical axes, one per panel (defaults to the names used
  in the paper).

- legend:

  Draw a legend above the panels. Defaults to `TRUE` when there is more
  than one group or a `truth` curve.

- include_zero:

  Include zero in the vertical range of every panel.

- mfrow:

  Panel layout `c(rows, columns)`; defaults to a single row.

- ...:

  Further graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`NULL`, invisibly. Called for its side effect.

## See also

[`plot.trajchain()`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md),
[`plot.trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md).

## Examples

``` r
simA <- simulate_cohort(n = 500, shape = "monotone", seed = 1)
simB <- simulate_cohort(n = 500, shape = "ushape", seed = 2)
fitA <- trajchain(simA$Y, simA$entry_age, simA$visit_gaps, tau0 = 0, tau1 = 1)
fitB <- trajchain(simB$Y, simB$entry_age, simB$visit_gaps, tau0 = 0, tau1 = 1)
plot_trajectories(list(Monotone = fitA, "U-shaped" = fitB),
                  which = c("shape", "mean"))

```
