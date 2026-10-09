# Input validation and the internal "design" object shared by all estimators.
#
# Internally every computation is carried out on the unit time scale
#   s = (t - tau0) / (tau1 - tau0),
# so that the grid is u_j = (j - 0.5) / m, j = 1, ..., m, exactly as in the
# simulation code of the paper (where tau0 = 0 and tau1 = 1). Results are
# mapped back to the original age scale before they are returned.

.check_Y <- function(Y) {
  if (is.data.frame(Y)) Y <- as.matrix(Y)
  if (!is.matrix(Y) || !is.numeric(Y)) {
    .stopf("`Y` must be a numeric matrix (one row per subject, one column per scheduled visit).")
  }
  if (ncol(Y) < 2L) {
    .stopf("`Y` must have at least two columns: the baseline visit and at least one follow-up visit.")
  }
  if (nrow(Y) < 2L) .stopf("`Y` must contain at least two subjects.")
  storage.mode(Y) <- "double"
  if (any(is.nan(Y) | is.infinite(Y))) {
    .stopf("`Y` contains NaN or infinite values; use NA for unobserved visits.")
  }
  if (anyNA(Y[, 1L])) {
    .stopf(paste0(
      "The first column of `Y` (the baseline measurement at study entry) must be ",
      "observed for every subject; %d subject(s) have a missing baseline value."),
      sum(is.na(Y[, 1L])))
  }
  if (any(Y < 0, na.rm = TRUE)) {
    .stopf("The marker must be non-negative (the model assumes a positive marker); `Y` has negative values.")
  }
  Y
}

.check_entry_age <- function(entry_age, n) {
  if (!is.numeric(entry_age) || length(entry_age) != n) {
    .stopf("`entry_age` must be a numeric vector with one value per row of `Y` (%d).", n)
  }
  if (any(!is.finite(entry_age))) .stopf("`entry_age` must not contain missing or infinite values.")
  as.numeric(entry_age)
}

.check_gaps <- function(visit_gaps, k) {
  if (!is.numeric(visit_gaps) || !length(visit_gaps) || any(!is.finite(visit_gaps)) ||
      any(visit_gaps <= 0)) {
    .stopf("`visit_gaps` must contain positive, finite numbers.")
  }
  if (length(visit_gaps) == 1L) visit_gaps <- rep(visit_gaps, k)
  if (length(visit_gaps) != k) {
    .stopf(paste0(
      "`visit_gaps` must have length 1 (equally spaced visits) or ncol(Y) - 1 = %d ",
      "(the scheduled time between consecutive visits); it has length %d."),
      k, length(visit_gaps))
  }
  as.numeric(visit_gaps)
}

# Build the design. `warn` controls the warning about intermittent missingness.
.prepare_design <- function(Y, entry_age, visit_gaps, tau0, tau1 = NULL, m = NULL,
                            grid_step = NULL, warn = TRUE) {
  Y <- .check_Y(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  k <- p - 1L
  A <- .check_entry_age(entry_age, n)
  gaps <- .check_gaps(visit_gaps, k)
  tau0 <- .check_number(tau0, "tau0")

  # ---- grid step: by default the greatest common divisor of the gaps ----
  if (is.null(grid_step)) {
    step <- .gcd_numeric(gaps)
  } else {
    step <- .check_number(grid_step, "grid_step", positive = TRUE)
  }
  units <- gaps / step
  if (any(abs(units - round(units)) > 1e-6 * pmax(1, units))) {
    .stopf(paste0(
      "Every visit gap must be an integer multiple of the grid step (%s). ",
      "Gaps in grid units: %s."), .format_num(step),
      paste(.format_num(units), collapse = ", "))
  }
  units <- as.integer(round(units))
  if (is.null(grid_step) && max(units) > 50L) {
    .stopf(paste0(
      "The greatest common divisor of the visit gaps (%s) is much smaller than the gaps ",
      "themselves (up to %d times), which would create a very fine grid. Supply the ",
      "scheduled (rounded) visit gaps rather than observed ones, or set `grid_step` ",
      "explicitly if this grid is intended."), .format_num(step), max(units))
  }

  # ---- age interval [tau0, tau1] and number of grid points m ----
  if (is.null(tau1) && is.null(m)) {
    .stopf("Supply the end of the age interval, `tau1`, or the number of grid points, `m`.")
  }
  if (!is.null(m)) {
    m <- .check_count(m, "m", min = 2L)
    tau1_m <- tau0 + m * step
    if (!is.null(tau1)) {
      tau1 <- .check_number(tau1, "tau1")
      if (abs(tau1 - tau1_m) > 1e-6 * max(1, abs(tau1))) {
        .stopf("`tau1` and `m` are inconsistent: with grid step %s, m = %d implies tau1 = %s.",
               .format_num(step), m, .format_num(tau1_m, 8))
      }
    } else {
      tau1 <- tau1_m
    }
  } else {
    tau1 <- .check_number(tau1, "tau1")
    if (tau1 <= tau0) .stopf("`tau1` must be larger than `tau0`.")
    m_real <- (tau1 - tau0) / step
    if (abs(m_real - round(m_real)) > 1e-6 * max(1, m_real)) {
      lo <- max(1, floor(m_real))
      .stopf(paste0(
        "(tau1 - tau0) must be an integer multiple of the grid step (%s). ",
        "For example, use tau1 = %s (m = %d) or tau1 = %s (m = %d)."),
        .format_num(step), .format_num(tau0 + lo * step, 8), as.integer(lo),
        .format_num(tau0 + (lo + 1) * step, 8), as.integer(lo + 1))
    }
    m <- as.integer(round(m_real))
  }
  if (m < 2L) .stopf("The grid must contain at least two points (m >= 2).")
  tau <- tau1 - tau0

  # ---- scheduled visit times (time since entry) ----
  equal <- all(abs(gaps - gaps[1L]) <= 1e-12 * gaps[1L])
  offsets <- if (equal) gaps[1L] * (0:k) else c(0, cumsum(gaps))
  age <- outer(A, offsets, "+")
  s_visit <- (age - tau0) / tau
  obs <- !is.na(Y)

  # ---- monotone missingness check: a missing visit followed by an observed one ----
  intermittent <- apply(obs, 1L, function(o) any(diff(o) > 0))
  if (warn && any(intermittent)) {
    .warnf(paste0(
      "%d subject(s) have intermittent missing visits (a missing visit followed by an ",
      "observed one). The method assumes monotone dropout; consecutive-visit pairs with ",
      "a missing visit are skipped."), sum(intermittent))
  }

  # ---- observed pairs of consecutive visits, grouped by gap (grid units) ----
  pairs <- list()
  for (g in sort(unique(units))) {
    cols <- which(units == g) + 1L
    u <- ycur <- ylag <- numeric(0)
    id <- integer(0)
    for (cc in cols) {
      ok <- obs[, cc] & obs[, cc - 1L]
      u <- c(u, s_visit[ok, cc])
      ycur <- c(ycur, Y[ok, cc])
      ylag <- c(ylag, Y[ok, cc - 1L])
      id <- c(id, which(ok))
    }
    pairs[[as.character(g)]] <- list(gap = g, u = u, ycur = ycur, ylag = ylag, id = id)
  }

  structure(list(
    Y = Y, entry_age = A, obs = obs,
    n = n, p = p, k = k,
    visit_gaps = gaps, offsets = offsets, equal_gaps = equal,
    grid_step = tau / m, gap_units = units,
    tau0 = tau0, tau1 = tau1, tau = tau, m = m,
    t_grid = tau0 + (seq_len(m) - 0.5) * (tau / m),
    s_grid = (seq_len(m) - 0.5) / m,
    s_visit = s_visit,
    pairs = pairs,
    pooled = list(u = s_visit[obs], y = Y[obs]),
    baseline = list(u = s_visit[, 1L], y = Y[, 1L]),
    n_obs = sum(obs),
    n_intermittent = sum(intermittent)
  ), class = "trajchain_design")
}

# Subset a design by subject (row) indices, e.g. for bootstrap or CV folds.
.subset_design <- function(design, rows) {
  .prepare_design(design$Y[rows, , drop = FALSE], design$entry_age[rows],
                  design$visit_gaps, design$tau0, tau1 = design$tau1,
                  grid_step = design$grid_step, warn = FALSE)
}

#' Convert long-format longitudinal data to a visit matrix
#'
#' Reshapes long-format data (one row per measurement) into the subject-by-visit
#' matrix used by [trajchain()], with one row per subject and one column per
#' scheduled visit. Unobserved visits (after death, dropout or a missed visit)
#' are filled with `NA`.
#'
#' @param data A data frame in long format.
#' @param id Name of the column identifying subjects.
#' @param visit Name of the column with the scheduled visit number. Visits must
#'   be numbered by consecutive integers; the smallest number in the data (for
#'   example 0 or 1) is taken to be the baseline visit at study entry.
#' @param value Name of the column with the marker values.
#'
#' @return A numeric matrix with one row per subject (row names are the subject
#'   identifiers, in order of first appearance in `data`) and one column per
#'   scheduled visit (column names `"visit_<number>"`).
#'
#' @seealso [trajchain()], [simcohort] for an example data set in long format.
#'
#' @examples
#' long <- data.frame(
#'   id = c(1, 1, 1, 2, 2, 3),
#'   visit = c(0, 1, 2, 0, 1, 0),
#'   marker = c(1.2, 1.4, 1.5, 0.8, 0.9, 1.1)
#' )
#' visit_matrix(long, id = "id", visit = "visit", value = "marker")
#'
#' @export
visit_matrix <- function(data, id, visit, value) {
  if (!is.data.frame(data)) .stopf("`data` must be a data frame.")
  for (nm in c(id, visit, value)) {
    if (!is.character(nm) || length(nm) != 1L || !nm %in% names(data)) {
      .stopf("Column `%s` was not found in `data`.", as.character(nm)[1L])
    }
  }
  vis <- data[[visit]]
  if (!is.numeric(vis) || anyNA(vis) || any(abs(vis - round(vis)) > 1e-8)) {
    .stopf("Column `%s` must contain whole visit numbers without missing values.", visit)
  }
  vis <- as.integer(round(vis))
  ids <- data[[id]]
  if (anyNA(ids)) .stopf("Column `%s` must not contain missing values.", id)
  key <- paste(ids, vis, sep = "\r")
  if (anyDuplicated(key)) {
    .stopf("Some subjects have more than one row for the same visit in `%s`.", visit)
  }
  val <- data[[value]]
  if (!is.numeric(val)) .stopf("Column `%s` must be numeric.", value)
  uid <- unique(ids)
  levels_v <- seq.int(min(vis), max(vis))
  out <- matrix(NA_real_, nrow = length(uid), ncol = length(levels_v),
                dimnames = list(as.character(uid), paste0("visit_", levels_v)))
  out[cbind(match(ids, uid), vis - min(vis) + 1L)] <- as.numeric(val)
  out
}
