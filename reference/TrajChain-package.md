# TrajChain: nonparametric marker trajectories under death and delayed entry

Population-level biomarker trajectories are usually summarized by a
single age-indexed curve estimated from a longitudinal cohort study that
enrolls individuals at different ages and follows them over time. Two
features of such studies distort naive curves: deaths during follow-up
(biomarkers are only observed in survivors, and death may be associated
with the biomarker through unobserved factors), and delayed entry (only
individuals alive at enrollment are sampled; when earlier-born and
later-born cohorts differ, this left truncation is informative).
`TrajChain` implements the nonparametric framework of Hu, Wu, Zhao and
Sun, in which a latent variable \\Z^0\\ captures individual
heterogeneity, may be associated with both death and birth cohort, and
multiplies a common shape function \\f(t)\\: \\E\\Y^0(t) \mid Z^0, T^0 =
s\\ = Z^0 f(t)\\ for \\t \le s\\.

## Workflow

1.  Arrange the data as a subject-by-visit matrix
    ([`visit_matrix()`](https://pingbohu43.github.io/TrajChain/reference/visit_matrix.md))
    with the entry ages and the scheduled gaps between visits.

2.  Optionally choose the bandwidth constant by cross-validation
    ([`cv_bandwidth()`](https://pingbohu43.github.io/TrajChain/reference/cv_bandwidth.md)).

3.  Fit the estimators with
    [`trajchain()`](https://pingbohu43.github.io/TrajChain/reference/trajchain.md):
    the shape function by kernel estimation of ratios and ratio chaining
    (equally spaced visits) or by penalized ratio matching (widely or
    unequally spaced visits, with the penalty chosen by
    [`cv_lambda()`](https://pingbohu43.github.io/TrajChain/reference/cv_lambda.md));
    then the unconditional mean trajectory and the birth-cohort-survivor
    bias.

4.  Quantify uncertainty with the subject-level bootstrap
    ([`trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/trajchain_boot.md))
    and display the results with
    [`plot.trajchain()`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md),
    [`plot.trajchain_boot()`](https://pingbohu43.github.io/TrajChain/reference/plot-trajchain.md)
    or
    [`plot_trajectories()`](https://pingbohu43.github.io/TrajChain/reference/plot_trajectories.md).

[`simulate_cohort()`](https://pingbohu43.github.io/TrajChain/reference/simulate_cohort.md)
and
[`sim_truth()`](https://pingbohu43.github.io/TrajChain/reference/sim_truth.md)
reproduce the simulation designs of the paper, and
[simcohort](https://pingbohu43.github.io/TrajChain/reference/simcohort.md)
is a simulated example data set on the age scale. See
[`vignette("TrajChain")`](https://pingbohu43.github.io/TrajChain/articles/TrajChain.md)
for a worked analysis.

## References

Hu, P., Wu, Y., Zhao, Y. and Sun, Y. Nonparametric estimation of marker
trajectories in the presence of death and delayed entry. Manuscript
under revision at *Biometrika*.

## See also

Useful links:

- <https://github.com/pingbohu43/TrajChain>

- <https://pingbohu43.github.io/TrajChain/>

- Report bugs at <https://github.com/pingbohu43/TrajChain/issues>

## Author

**Maintainer**: Pingbo Hu <pingbo.hu@gmail.com> \[copyright holder\]
