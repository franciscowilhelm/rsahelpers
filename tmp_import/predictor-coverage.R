# predictor-coverage.R -- co-occurrence diagnostics for the sparse-region /
# bag-plot caveat, and the cell-proportions heatmap from report 07
# (07-predictor-complications).
#
# Exported from the rsasim project. predictor_coverage() is lib/simulate.R;
# the Likert predictor generator it is usually run on is carried over too;
# plot_predictor_coverage() is a NEW function that packages the inline ggplot
# from qmd/07-predictor-complications.qmd as a reusable plotter.
#
# Background (Nestler et al. 2019; Humberg et al. 2022): don't interpret the
# regression surface where there are no data points, and report how many cases
# populate the region a hypothesis is about. Sign convention: deficiency =
# P > E, excess = E > P.

# ---- Likert predictor generator ---------------------------------------------

#' Category-probability presets for 1-5 Likert margins
#'
#' Chosen to bracket the marginal distributions typical of P/E scales (floor
#' effects are rare, so no floor presets):
#' `normal` (symmetric, ~7% at each endpoint), `mild_ceiling`, `strong_ceiling`.
#' @export
likert_margin_presets <- list(
  normal         = c(.07, .24, .38, .24, .07),
  mild_ceiling   = c(.03, .10, .25, .38, .24),
  strong_ceiling = c(.01, .05, .14, .35, .45)
)

#' Latent normal thresholds for a Likert margin
#'
#' @param scale Integer scale points (default `1:5`).
#' @param thresholds A preset name (see [likert_margin_presets]) or numeric
#'   cutpoints of length `length(scale) - 1`. Uniform margins are rejected.
#' @return Numeric cutpoints on the latent standard-normal scale.
#' @export
likert_thresholds <- function(scale = 1:5, thresholds = "normal") {
  k <- length(scale)
  if (is.character(thresholds)) {
    if (identical(thresholds, "uniform")) {
      stop(
        "Uniform Likert margins are not supported: equal category ",
        "probabilities do not occur in empirical P/E scales. Use one of the ",
        "presets ('", paste(names(likert_margin_presets), collapse = "', '"),
        "') or pass numeric cutpoints."
      )
    }
    probs <- likert_margin_presets[[thresholds]]
    if (is.null(probs)) {
      stop("Unknown margin preset '", thresholds, "'. Available: '",
           paste(names(likert_margin_presets), collapse = "', '"), "'.")
    }
    if (k != length(probs)) {
      stop("Margin presets are defined for a ", length(probs),
           "-point scale; pass numeric cutpoints for other scales.")
    }
    return(stats::qnorm(cumsum(probs)[seq_len(k - 1)]))
  }
  stopifnot(is.numeric(thresholds), length(thresholds) == k - 1, !is.unsorted(thresholds))
  thresholds
}

#' Simulate correlated person-environment predictors
#'
#' Predictors come from a latent bivariate standard normal with correlation
#' `rho`, optionally discretized to a Likert scale via [likert_thresholds].
#' `rho` is the *latent* correlation; discretization attenuates the manifest
#' Likert correlation by roughly 10%. Predictors are pooled-centered at the
#' scale midpoint.
#'
#' @param n Sample size.
#' @param rho Latent correlation between P and E (`abs(rho) < 1`).
#' @param type `"likert"` (discretized) or `"continuous"` (latent, shifted to
#'   the scale midpoint; not clamped).
#' @param thresholds Preset name or numeric cutpoints (see [likert_thresholds]).
#' @param scale Integer scale points (default `1:5`).
#' @return A [tibble::tibble] with `P_raw`, `E_raw` (raw scale) and centered
#'   `P`, `E`.
#' @export
simulate_pe_predictors <- function(n, rho = 0,
                                   type = c("likert", "continuous"),
                                   thresholds = "normal",
                                   scale = 1:5) {
  type <- match.arg(type)
  stopifnot(abs(rho) < 1)
  sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  latent <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma)
  center_value <- mean(scale)

  if (type == "likert") {
    cuts <- likert_thresholds(scale, thresholds)
    p_raw <- scale[findInterval(latent[, 1], cuts) + 1]
    e_raw <- scale[findInterval(latent[, 2], cuts) + 1]
  } else {
    p_raw <- center_value + latent[, 1]
    e_raw <- center_value + latent[, 2]
  }

  tibble::tibble(
    P_raw = p_raw,
    E_raw = e_raw,
    P = p_raw - center_value,
    E = e_raw - center_value
  )
}

# ---- Coverage diagnostics ----------------------------------------------------

#' Co-occurrence coverage of the raw P x E grid
#'
#' Returns the raw P x E cell proportions plus the shares of the regions RSA
#' hypotheses live in. Sign convention: deficiency = P > E, excess = E > P.
#'
#' @param data A data frame with raw predictors `P_raw`, `E_raw`.
#' @param scale Integer scale points (default `1:5`).
#' @return A list with `cell_props` (a proportions [table] over P x E) and
#'   `summary` (a one-row tibble: manifest correlation, discrepancy-region
#'   shares `deficiency_ge1/ge2`, `excess_ge1/ge2`, `low_low`/`high_high` LOC
#'   ends, `cells_below_1pct`, `min_cell`).
#' @export
predictor_coverage <- function(data, scale = 1:5) {
  n <- nrow(data)
  cells <- table(
    P = factor(data$P_raw, levels = scale),
    E = factor(data$E_raw, levels = scale)
  ) / n
  disc <- data$P_raw - data$E_raw
  lo <- scale[2]
  hi <- scale[length(scale) - 1]
  list(
    cell_props = cells,
    summary = tibble::tibble(
      r_manifest = stats::cor(data$P_raw, data$E_raw),
      deficiency_ge1 = mean(disc >= 1),
      deficiency_ge2 = mean(disc >= 2),
      excess_ge1 = mean(disc <= -1),
      excess_ge2 = mean(disc <= -2),
      low_low = mean(data$P_raw <= lo & data$E_raw <= lo),
      high_high = mean(data$P_raw >= hi & data$E_raw >= hi),
      cells_below_1pct = sum(cells < 0.01),
      min_cell = min(cells)
    )
  )
}

#' Long-format cell proportions from one or more coverage results
#'
#' Reshapes the `cell_props` table(s) into a tidy data frame suitable for
#' ggplot. Pass a single [predictor_coverage] result, its `cell_props` table,
#' or a named list of coverage results (the names become a `group` column for
#' faceting).
#'
#' @param coverage A coverage result, a `cell_props` table, or a named list of
#'   coverage results.
#' @return A [tibble::tibble] with integer `P`, `E`, proportion `n`, and (for a
#'   named list) a `group` column.
#' @export
coverage_cells_long <- function(coverage) {
  one <- function(cov) {
    props <- if (is.list(cov) && !is.null(cov$cell_props)) cov$cell_props else cov
    df <- tibble::as_tibble(props)
    df$P <- as.integer(as.character(df$P))
    df$E <- as.integer(as.character(df$E))
    df
  }
  if (is.list(coverage) && is.null(coverage$cell_props) &&
      !inherits(coverage, "table")) {
    purrr::imap_dfr(coverage, function(cov, nm) {
      out <- one(cov)
      out$group <- nm
      out
    })
  } else {
    one(coverage)
  }
}

#' Cell-proportions heatmap of the raw P x E grid
#'
#' The reusable form of the report-07 figure. Tiles the raw grid by proportion
#' of cases; cells below `min_prop` (the regions a bag plot would exclude) are
#' marked with a red x.
#'
#' @param coverage A [predictor_coverage] result, a `cell_props` table, or a
#'   named list of coverage results (each panel faceted, keyed by name).
#' @param min_prop Sparse-cell threshold for the x marks (default `0.01`,
#'   ~5 expected cases at n = 500).
#' @param mark_sparse Logical; draw the x marks (default `TRUE`).
#' @return A [ggplot2::ggplot] object.
#' @export
plot_predictor_coverage <- function(coverage, min_prop = 0.01,
                                     mark_sparse = TRUE) {
  cells <- coverage_cells_long(coverage)
  p <- ggplot2::ggplot(cells, ggplot2::aes(.data$E, .data$P, fill = .data$n)) +
    ggplot2::geom_tile()
  if (mark_sparse) {
    p <- p + ggplot2::geom_point(
      data = cells[cells$n < min_prop, ],
      shape = 4, size = 2, color = "red"
    )
  }
  if ("group" %in% names(cells)) {
    p <- p + ggplot2::facet_wrap(~ .data$group)
  }
  p +
    ggplot2::scale_fill_viridis_c(name = "Proportion", trans = "sqrt") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "E (raw)", y = "P (raw)") +
    ggplot2::theme_minimal()
}
