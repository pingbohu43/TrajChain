#' Simulated cohort with death and delayed entry on the age scale
#'
#' A simulated longitudinal cohort in long format, used in the examples and the
#' vignette. Subjects enter the study at different ages, have a baseline visit
#' and up to two follow-up visits scheduled every 16 months, and are lost to
#' follow-up through death or dropout. The data were generated with
#' [simulate_cohort()] under the informative-truncation design of the paper
#' (so that both survivor and birth-cohort selection are present), with the
#' unit age scale mapped to ages 30 to 59.33, i.e. a grid of \eqn{m = 22}
#' points with step 16 months, as in the data analysis of the paper. The
#' marker in group A follows the monotone shape function [shape_monotone()] and
#' in group B the U-shaped function [shape_ushape()], rescaled to the age
#' interval. These are simulated data, not observations on real individuals.
#'
#' @format A data frame with 1502 rows (one per observed measurement) and 6
#'   variables:
#'   \describe{
#'     \item{`id`}{subject identifier (1 to 600);}
#'     \item{`group`}{factor with levels `A` and `B` (300 subjects each);}
#'     \item{`visit`}{scheduled visit number: 0 (baseline), 1 or 2;}
#'     \item{`age`}{age at the visit, in years;}
#'     \item{`entry_age`}{age at study entry, in years;}
#'     \item{`marker`}{the (positive) marker value.}
#'   }
#'
#' @source Simulated; see `data-raw/simcohort.R` in the package sources.
#'
#' @examples
#' head(simcohort)
#' groupA <- subset(simcohort, group == "A")
#' Y <- visit_matrix(groupA, id = "id", visit = "visit", value = "marker")
#' A <- groupA$entry_age[match(rownames(Y), groupA$id)]
#' fit <- trajchain(Y, A, visit_gaps = 16 / 12, tau0 = 30, m = 22)
#' fit
"simcohort"
