#' TrajChain: nonparametric marker trajectories under death and delayed entry
#'
#' Population-level biomarker trajectories are usually summarized by a single
#' age-indexed curve estimated from a longitudinal cohort study that enrolls
#' individuals at different ages and follows them over time. Two features of
#' such studies distort naive curves: deaths during follow-up (biomarkers are
#' only observed in survivors, and death may be associated with the biomarker
#' through unobserved factors), and delayed entry (only individuals alive at
#' enrollment are sampled; when earlier-born and later-born cohorts differ,
#' this left truncation is informative). `TrajChain` implements the
#' nonparametric framework of Hu, Wu, Zhao and Sun, in which a latent variable
#' \eqn{Z^0} captures individual heterogeneity, may be associated with both death
#' and birth cohort, and multiplies a common shape function \eqn{f(t)}:
#' \eqn{E\{Y^0(t) \mid Z^0, T^0 = s\} = Z^0 f(t)}{E{Y0(t) | Z0, T0 = s} = Z0 f(t)}
#' for \eqn{t \le s}.
#'
#' @section Workflow:
#' 1. Arrange the data as a subject-by-visit matrix ([visit_matrix()]) with
#'    the entry ages and the scheduled gaps between visits.
#' 2. Optionally choose the bandwidth constant by cross-validation
#'    ([cv_bandwidth()]).
#' 3. Fit the estimators with [trajchain()]: the shape function by kernel
#'    estimation of ratios and ratio chaining (equally spaced visits) or by
#'    penalized ratio matching (widely or unequally spaced visits, with the
#'    penalty chosen by [cv_lambda()]); then the unconditional mean trajectory
#'    and the birth-cohort-survivor bias.
#' 4. Quantify uncertainty with the subject-level bootstrap
#'    ([trajchain_boot()]) and display the results with [plot.trajchain()],
#'    [plot.trajchain_boot()] or [plot_trajectories()].
#'
#' [simulate_cohort()] and [sim_truth()] reproduce the simulation designs of the
#' paper, and [simcohort] is a simulated example data set on the age scale.
#' See `vignette("TrajChain")` for a worked analysis.
#'
#' @template ref
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
