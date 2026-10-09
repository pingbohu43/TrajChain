# Code that creates the example data set `simcohort` (run from the package root).
#
# Two groups of 300 subjects simulated with simulate_cohort() under the
# informative-truncation design of the paper (birth-cohort effect), with a
# baseline visit and two follow-up visits scheduled 16 months apart. The unit
# age scale [0, 1] is mapped to ages [30, 30 + 22 * 16 / 12] = [30, 59.33], so
# that the grid of the analysis has m = 22 points, as in the ABC-DS analysis
# of the paper. Group A follows the monotone shape function and group B the
# U-shaped one.
devtools::load_all()

tau0 <- 30
gap_years <- 16 / 12
m <- 22
scale <- m * gap_years

make_group <- function(shape, seed, group, id_offset) {
  sim <- simulate_cohort(n = 300, shape = shape, visit_gaps = rep(1 / m, 2),
                         truncation = "informative", seed = seed)
  entry_age <- tau0 + scale * sim$entry_age
  Y <- sim$Y
  long <- data.frame(
    id = rep(id_offset + seq_len(nrow(Y)), times = ncol(Y)),
    group = group,
    visit = rep(seq_len(ncol(Y)) - 1L, each = nrow(Y)),
    entry_age = rep(entry_age, times = ncol(Y)),
    marker = as.vector(Y)
  )
  long$age <- long$entry_age + long$visit * gap_years
  long <- long[!is.na(long$marker), c("id", "group", "visit", "age", "entry_age", "marker")]
  long[order(long$id, long$visit), ]
}

simcohort <- rbind(make_group("monotone", 2026, "A", 0L),
                   make_group("ushape", 2027, "B", 300L))
simcohort$group <- factor(simcohort$group)
simcohort$age <- round(simcohort$age, 4)
simcohort$entry_age <- round(simcohort$entry_age, 4)
simcohort$marker <- round(simcohort$marker, 5)
rownames(simcohort) <- NULL

usethis::use_data(simcohort, overwrite = TRUE)
