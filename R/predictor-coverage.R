#' Diagnose joint predictor coverage for response-surface analysis
#'
#' @param data A data frame containing the two commensurate predictors.
#' @param x,y Predictor columns, supplied as unquoted names or strings.
#' @param breaks Either the number of shared equal-width bins (default `10`) or
#'   a numeric vector of common breakpoints for both predictors.
#' @param discrepancy Optional positive cutoffs at which to summarize the
#'   proportions of observations with `X - Y` or `Y - X` at least that large.
#' @param min_count Minimum cell count considered adequately covered.
#'
#' @return An object of class `predictor_coverage`. It contains `cells`, a
#'   one-row `summary`, an optional `discrepancy` table, the shared `breaks`,
#'   and plotting metadata.
#' @details
#' Shared breaks preserve the geometry of the line `X = Y`. This is especially
#' useful for scale scores derived from Likert items: using every observed
#' decimal value creates an unnecessarily sparse grid, while separate quantile
#' bins distort congruence and discrepancy regions.
#' @export
#'
#' @examples
#' scores <- data.frame(
#'   self = c(1.0, 1.5, 2.2, 3.1, 3.8, 4.4, 4.8),
#'   other = c(1.2, 1.6, 2.7, 2.9, 3.5, 4.1, 4.7)
#' )
#' coverage <- predictor_coverage(
#'   scores,
#'   self,
#'   other,
#'   breaks = 4,
#'   discrepancy = c(0.5, 1)
#' )
#' coverage$summary
predictor_coverage <- function(
  data,
  x,
  y,
  breaks = 10L,
  discrepancy = NULL,
  min_count = 5L
) {
  x_name <- resolve_public_variable(x, substitute(x), "x")
  y_name <- resolve_public_variable(y, substitute(y), "y")
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  missing_columns <- setdiff(c(x_name, y_name), names(data))
  if (length(missing_columns) > 0L) {
    rlang::abort(sprintf(
      "Missing predictor column%s: %s.",
      if (length(missing_columns) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_columns), collapse = ", ")
    ))
  }
  x_values <- data[[x_name]]
  y_values <- data[[y_name]]
  if (!is.numeric(x_values) || !is.numeric(y_values)) {
    rlang::abort("`x` and `y` must select numeric columns.")
  }
  validate_coverage_min_count(min_count)
  discrepancy <- validate_discrepancy_cutoffs(discrepancy)

  finite <- is.finite(x_values) & is.finite(y_values)
  n_total <- length(x_values)
  n_complete <- sum(finite)
  if (n_complete < 2L) {
    rlang::abort("At least two finite X/Y pairs are required.")
  }
  x_complete <- x_values[finite]
  y_complete <- y_values[finite]
  pooled_range <- range(c(x_complete, y_complete))
  if (diff(pooled_range) <= 0) {
    rlang::abort("The pooled X/Y values must contain variation.")
  }
  shared_breaks <- coverage_breaks(breaks, pooled_range)
  n_bins <- length(shared_breaks) - 1L

  x_bin <- cut(
    x_complete,
    breaks = shared_breaks,
    labels = FALSE,
    include.lowest = TRUE
  )
  y_bin <- cut(
    y_complete,
    breaks = shared_breaks,
    labels = FALSE,
    include.lowest = TRUE
  )
  cell_id <- x_bin + (y_bin - 1L) * n_bins
  counts <- tabulate(cell_id, nbins = n_bins^2)
  cells <- expand.grid(
    x_bin = seq_len(n_bins),
    y_bin = seq_len(n_bins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  cells$count <- counts
  cells$proportion <- counts / n_complete
  cells$x_low <- shared_breaks[cells$x_bin]
  cells$x_high <- shared_breaks[cells$x_bin + 1L]
  cells$y_low <- shared_breaks[cells$y_bin]
  cells$y_high <- shared_breaks[cells$y_bin + 1L]
  cells$x_mid <- (cells$x_low + cells$x_high) / 2
  cells$y_mid <- (cells$y_low + cells$y_high) / 2
  cells$sparse <- cells$count < min_count

  difference <- x_complete - y_complete
  correlation <- if (stats::sd(x_complete) == 0 || stats::sd(y_complete) == 0) {
    NA_real_
  } else {
    stats::cor(x_complete, y_complete)
  }
  summary <- data.frame(
    x = x_name,
    y = y_name,
    n = n_complete,
    n_missing = n_total - n_complete,
    correlation = correlation,
    x_greater_than_y = mean(difference > 0),
    y_greater_than_x = mean(difference < 0),
    ties = mean(difference == 0),
    empty_cells = sum(cells$count == 0L),
    sparse_cells = sum(cells$sparse),
    min_cell_count = min(cells$count),
    min_count = as.integer(min_count),
    stringsAsFactors = FALSE
  )

  discrepancy_table <- coverage_discrepancy_table(
    difference,
    discrepancy,
    n_complete
  )
  structure(
    list(
      cells = cells,
      summary = summary,
      discrepancy = discrepancy_table,
      breaks = shared_breaks,
      metadata = list(
        x = x_name,
        y = y_name,
        min_count = as.integer(min_count)
      )
    ),
    class = "predictor_coverage"
  )
}

#' Plot joint predictor coverage
#'
#' @param coverage A [`predictor_coverage()`] result or a named list of results.
#' @param mark_sparse If `TRUE`, mark cells below the result's minimum count.
#' @param show_diagonal If `TRUE`, draw the line `X = Y`.
#'
#' @return A [`ggplot2::ggplot`] object.
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   scores <- data.frame(
#'     x = c(1, 1.5, 2, 3, 4, 4.5),
#'     y = c(1.2, 1.4, 2.5, 2.8, 3.7, 4.6)
#'   )
#'   plot_predictor_coverage(predictor_coverage(scores, x, y, breaks = 4))
#' }
plot_predictor_coverage <- function(
  coverage,
  mark_sparse = TRUE,
  show_diagonal = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    rlang::abort("Package `ggplot2` is required to plot predictor coverage.")
  }
  prepared <- coverage_plot_data(coverage)
  cells <- prepared$cells

  plot <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = x_mid, y = y_mid, fill = proportion)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Proportion", trans = "sqrt") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = prepared$x_label,
      y = prepared$y_label
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(mark_sparse)) {
    plot <- plot +
      ggplot2::geom_point(
        data = cells[cells$sparse, , drop = FALSE],
        shape = 4,
        size = 2,
        colour = "red"
      )
  }
  if (isTRUE(show_diagonal)) {
    plot <- plot +
      ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed",
        colour = "grey30"
      )
  }
  if ("group" %in% names(cells)) {
    plot <- plot + ggplot2::facet_wrap(~group)
  }
  plot
}

coverage_breaks <- function(breaks, pooled_range) {
  if (
    is.numeric(breaks) &&
      length(breaks) == 1L &&
      !is.na(breaks) &&
      is.finite(breaks) &&
      breaks >= 1 &&
      breaks == as.integer(breaks)
  ) {
    return(seq(
      pooled_range[[1]],
      pooled_range[[2]],
      length.out = as.integer(breaks) + 1L
    ))
  }
  if (
    !is.numeric(breaks) ||
      length(breaks) < 2L ||
      anyNA(breaks) ||
      any(!is.finite(breaks)) ||
      is.unsorted(breaks, strictly = TRUE)
  ) {
    rlang::abort(
      "`breaks` must be a positive integer or strictly increasing numeric vector."
    )
  }
  if (
    pooled_range[[1]] < breaks[[1]] ||
      pooled_range[[2]] > breaks[[length(breaks)]]
  ) {
    rlang::abort("Explicit `breaks` must span all finite X and Y values.")
  }
  breaks
}

validate_coverage_min_count <- function(min_count) {
  if (
    !is.numeric(min_count) ||
      length(min_count) != 1L ||
      is.na(min_count) ||
      !is.finite(min_count) ||
      min_count < 1 ||
      min_count != as.integer(min_count)
  ) {
    rlang::abort("`min_count` must be one positive integer.")
  }
  invisible(min_count)
}

validate_discrepancy_cutoffs <- function(discrepancy) {
  if (is.null(discrepancy)) {
    return(NULL)
  }
  if (
    !is.numeric(discrepancy) ||
      anyNA(discrepancy) ||
      any(!is.finite(discrepancy)) ||
      any(discrepancy <= 0)
  ) {
    rlang::abort("`discrepancy` must contain finite positive numbers.")
  }
  sort(unique(discrepancy))
}

coverage_discrepancy_table <- function(difference, cutoffs, n) {
  if (is.null(cutoffs)) {
    return(NULL)
  }
  rows <- lapply(cutoffs, function(cutoff) {
    x_count <- sum(difference >= cutoff)
    y_count <- sum(difference <= -cutoff)
    data.frame(
      cutoff = cutoff,
      direction = c("x_greater_than_y", "y_greater_than_x"),
      count = c(x_count, y_count),
      proportion = c(x_count, y_count) / n,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

coverage_plot_data <- function(coverage) {
  if (inherits(coverage, "predictor_coverage")) {
    return(list(
      cells = coverage$cells,
      x_label = coverage$metadata$x,
      y_label = coverage$metadata$y
    ))
  }
  if (
    !is.list(coverage) ||
      length(coverage) == 0L ||
      is.null(names(coverage)) ||
      any(!nzchar(names(coverage))) ||
      !all(vapply(coverage, inherits, logical(1), "predictor_coverage"))
  ) {
    rlang::abort(
      "`coverage` must be a predictor-coverage result or a named list of them."
    )
  }
  cells <- do.call(
    rbind,
    lapply(seq_along(coverage), function(i) {
      result <- coverage[[i]]$cells
      result$group <- names(coverage)[[i]]
      result
    })
  )
  rownames(cells) <- NULL
  x_labels <- unique(vapply(coverage, function(x) x$metadata$x, character(1)))
  y_labels <- unique(vapply(coverage, function(x) x$metadata$y, character(1)))
  list(
    cells = cells,
    x_label = if (length(x_labels) == 1L) x_labels else "X",
    y_label = if (length(y_labels) == 1L) y_labels else "Y"
  )
}
