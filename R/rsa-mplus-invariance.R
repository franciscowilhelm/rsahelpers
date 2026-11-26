#' Test measurement invariance between the X and Y factors
#'
#' Writes, runs, and reads a sequence of Mplus confirmatory factor models that
#' test whether matched X and Y indicators measure their factors in the same
#' way: configural, metric (equal loadings), and scalar (equal loadings and
#' intercepts). Scalar invariance is what `identification = "invariance"` in
#' [create_rsa_mplus_model()] assumes when it reads the latent line X = Y as
#' equality on the response scale.
#'
#' @param data A data frame containing the indicators.
#' @param measurement A named list mapping the factors to their indicators.
#'   The two factors of a pair need the same number of indicators, in
#'   matching order.
#' @param pairs A list of character vectors of length two, each naming an X
#'   and a Y factor of `measurement`. With several pairs, all factors are
#'   estimated in one model with free covariances. The default pairs the two
#'   factors of a two-factor `measurement`.
#' @param free_intercepts Names of indicators whose intercept is not shared
#'   within its pair. Naming either member frees the pair. When given, a
#'   `"partial"` level is fitted between `"metric"` and `"scalar"`. At least
#'   one matched pair of each X/Y factor pair must keep its shared intercept.
#' @param blocks Named list of Mplus input blocks, as in
#'   [create_rsa_mplus_model()]. The defaults are `ESTIMATOR = MLR;` and
#'   `OUTPUT: STDYX;`.
#' @param usevariables Columns of `data` to export. Defaults to the
#'   indicators in `measurement`.
#' @param modelout Path of an Mplus `.inp` file that names the set. Each level
#'   is written next to it with its name appended, for example
#'   `inv_metric.inp` for `modelout = "inv.inp"`, and all levels read the
#'   data file `inv.dat`. On a case-sensitive file system the file name must
#'   be lower case (see [write_rsa_mplus_model()]).
#' @param overwrite If `TRUE`, replace existing input and data files.
#' @param run If `TRUE` (the default), run the models. If `FALSE`, only write
#'   the input and data files.
#' @param Mplus_command Optional path or command used to invoke Mplus.
#' @param quiet If `TRUE`, suppress status messages from MplusAutomation.
#'
#' @details
#' The levels are:
#' * `"configural"`: the same factor structure, all loadings and intercepts
#'   free (marker loadings fixed at 1, factor means at 0).
#' * `"metric"`: matched indicators share their loading.
#' * `"partial"`: metric, and matched indicators share their intercept except
#'   those in `free_intercepts`. The Y factor mean is freed.
#' * `"scalar"`: metric, and all matched indicators share their intercept.
#'
#' `comparison` tests each level against a less restricted one with the
#' Satorra-Bentler scaled chi-square difference: metric against configural,
#' and partial and scalar each against metric. `delta_CFI` is the change in
#' CFI relative to the same model.
#'
#' @return An object of class `rsa_mplus_invariance`, a list with
#'   * `models`: the named list of workflows, one per level.
#'   * `comparison`: a data frame with one row per level: `level`, `k`,
#'     `chisq`, `df`, `scaling`, `CFI`, `RMSEA`, `SRMR`, the level it is
#'     tested `against`, `delta_CFI`, and the scaled difference test
#'     (`diff_chisq`, `diff_df`, `diff_p`). It is `NULL` when `run = FALSE`.
#'   * `files`: the `input` paths, named by level, and the shared `data` path.
#' @seealso [create_rsa_mplus_model()]
#' @export
#'
#' @examples
#' \dontrun{
#' invariance <- fit_rsa_mplus_invariance(
#'   data = dat,
#'   measurement = list(
#'     X = c("x1", "x2", "x3"),
#'     Y = c("y1", "y2", "y3")
#'   ),
#'   free_intercepts = "x2",
#'   modelout = "invariance/inv.inp"
#' )
#' invariance$comparison
#' }
fit_rsa_mplus_invariance <- function(
  data,
  measurement,
  pairs = NULL,
  free_intercepts = NULL,
  blocks = list(),
  usevariables = NULL,
  modelout,
  overwrite = FALSE,
  run = TRUE,
  Mplus_command = NULL,
  quiet = TRUE
) {
  require_rsa_mplus_namespace("MplusAutomation")
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }

  design <- rsa_mplus_invariance_design(measurement, pairs, free_intercepts)
  if (is.null(usevariables)) {
    usevariables <- unlist(design$measurement, use.names = FALSE)
  }
  validate_rsa_mplus_variables(data, usevariables)
  missing_indicators <- setdiff(
    unlist(design$measurement, use.names = FALSE),
    usevariables
  )
  if (length(missing_indicators) > 0L) {
    rlang::abort(sprintf(
      "Measurement indicator%s %s not found in exported data.",
      if (length(missing_indicators) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_indicators), collapse = ", ")
    ))
  }

  levels <- c(
    "configural",
    "metric",
    if (length(unlist(design$free)) > 0L) "partial",
    "scalar"
  )
  files <- rsa_mplus_family_files(modelout, levels, overwrite)
  if (isTRUE(run)) {
    Mplus_command <- resolve_rsa_mplus_command(Mplus_command)
  }

  block_args <- validate_rsa_mplus_blocks(blocks)
  defaults <- list(ANALYSIS = "ESTIMATOR = MLR;", OUTPUT = "STDYX;")
  defaults[names(block_args)] <- block_args
  block_args <- defaults
  block_args$VARIABLE <- rsa_mplus_variable_block(
    usevariables = usevariables,
    extra = block_args$VARIABLE
  )

  models <- lapply(stats::setNames(levels, levels), function(level) {
    title <- block_args$TITLE
    if (is.null(title)) {
      title <- sprintf("X/Y measurement invariance: %s;", level)
    }
    args <- c(
      list(TITLE = title),
      block_args[setdiff(names(block_args), "TITLE")],
      list(
        MODEL = paste(
          rsa_mplus_invariance_syntax(design, level),
          collapse = "\n"
        ),
        usevariables = usevariables,
        rdata = data[, usevariables, drop = FALSE],
        autov = FALSE,
        quiet = quiet
      )
    )
    workflow <- new_rsa_mplus_workflow(
      spec = list(
        level = level,
        usevariables = usevariables,
        invariance = design[c("measurement", "pairs", "free")]
      ),
      mplus = do.call(MplusAutomation::mplusObject, args),
      status = "created"
    )
    # Existing files were checked for the whole set beforehand.
    workflow <- write_rsa_mplus_model(
      workflow,
      modelout = files$input[[level]],
      dataout = files$data,
      overwrite = TRUE,
      quiet = quiet
    )
    if (isTRUE(run)) {
      workflow <- run_rsa_mplus_model(
        workflow,
        Mplus_command = Mplus_command,
        quiet = quiet
      )
    }
    workflow
  })

  structure(
    list(
      models = models,
      comparison = if (isTRUE(run)) rsa_mplus_invariance_comparison(models),
      files = files
    ),
    class = "rsa_mplus_invariance"
  )
}

#' @export
print.rsa_mplus_invariance <- function(x, digits = 3, ...) {
  cat("<rsa_mplus_invariance>\n")
  if (is.null(x$comparison)) {
    cat("Levels written:", paste(names(x$models), collapse = ", "), "\n")
    return(invisible(x))
  }
  shown <- x$comparison
  for (column in c("chisq", "diff_chisq")) {
    shown[[column]] <- round(shown[[column]], 2L)
  }
  for (column in c("scaling", "CFI", "RMSEA", "SRMR", "delta_CFI")) {
    shown[[column]] <- round(shown[[column]], digits)
  }
  shown$diff_p <- ifelse(
    is.na(shown$diff_p),
    NA_character_,
    format.pval(shown$diff_p, digits = digits, eps = 0.001)
  )
  print(shown, row.names = FALSE, na.print = "")
  invisible(x)
}

# Validates the factors, pairs, and freed intercepts. `free` holds, per pair,
# the positions of the indicator pairs whose intercepts are not shared.
rsa_mplus_invariance_design <- function(
  measurement,
  pairs,
  free_intercepts,
  call = rlang::caller_env()
) {
  valid_measurement <- is.list(measurement) &&
    !is.data.frame(measurement) &&
    !is.null(names(measurement)) &&
    all(nzchar(names(measurement))) &&
    !anyDuplicated(toupper(names(measurement))) &&
    all(vapply(
      measurement,
      function(x) is.character(x) && length(x) > 0L && !anyNA(x),
      logical(1)
    ))
  if (!valid_measurement) {
    rlang::abort(
      "`measurement` must be a named list mapping factors to indicator names.",
      call = call
    )
  }
  if (anyDuplicated(unlist(measurement, use.names = FALSE))) {
    rlang::abort(
      "Each observed indicator may belong to only one measurement factor.",
      call = call
    )
  }

  if (is.null(pairs)) {
    if (length(measurement) != 2L) {
      rlang::abort(
        "`pairs` is needed when `measurement` does not have exactly two factors.",
        call = call
      )
    }
    pairs <- list(names(measurement))
  }
  if (is.character(pairs)) {
    pairs <- list(pairs)
  }
  valid_pairs <- is.list(pairs) &&
    length(pairs) > 0L &&
    all(vapply(
      pairs,
      function(x) is.character(x) && length(x) == 2L && !anyNA(x),
      logical(1)
    ))
  if (!valid_pairs) {
    rlang::abort(
      "`pairs` must be a list of character vectors naming an X and a Y factor.",
      call = call
    )
  }
  paired <- unlist(pairs, use.names = FALSE)
  if (anyDuplicated(paired) || !setequal(paired, names(measurement))) {
    rlang::abort(
      "Every factor in `measurement` must belong to exactly one of the `pairs`.",
      call = call
    )
  }
  for (pair in pairs) {
    sizes <- lengths(measurement[pair])
    if (sizes[[1]] != sizes[[2]] || sizes[[1]] < 2L) {
      rlang::abort(
        sprintf(
          "`%s` and `%s` need the same number of indicators (at least two), in matching order.",
          pair[[1]],
          pair[[2]]
        ),
        call = call
      )
    }
  }

  if (
    !is.null(free_intercepts) &&
      (!is.character(free_intercepts) || anyNA(free_intercepts))
  ) {
    rlang::abort(
      "`free_intercepts` must be a character vector of indicator names.",
      call = call
    )
  }
  unknown <- setdiff(free_intercepts, unlist(measurement, use.names = FALSE))
  if (length(unknown) > 0L) {
    rlang::abort(
      sprintf(
        "`free_intercepts` must name indicators in `measurement`; not found: %s.",
        paste(sprintf("`%s`", unknown), collapse = ", ")
      ),
      call = call
    )
  }
  free <- lapply(pairs, function(pair) {
    positions <- c(
      match(free_intercepts, measurement[[pair[[1]]]]),
      match(free_intercepts, measurement[[pair[[2]]]])
    )
    sort(unique(positions[!is.na(positions)]))
  })
  # The freed Y factor mean is identified by the intercepts that stay shared.
  for (index in seq_along(pairs)) {
    pair <- pairs[[index]]
    if (length(free[[index]]) == length(measurement[[pair[[1]]]])) {
      rlang::abort(
        sprintf(
          "`free_intercepts` frees every intercept of `%s` and `%s`; at least one matched pair must share its intercept.",
          pair[[1]],
          pair[[2]]
        ),
        call = call
      )
    }
  }

  list(measurement = measurement, pairs = pairs, free = free)
}

# MODEL statements of one invariance level. Labels carry the pair number, so
# several pairs can be estimated in one model.
rsa_mplus_invariance_syntax <- function(design, level) {
  one_pair <- function(index) {
    pair <- design$pairs[[index]]
    x_items <- design$measurement[[pair[[1]]]]
    y_items <- design$measurement[[pair[[2]]]]
    k <- seq_along(x_items)

    loading_labels <- rep("", length(k))
    if (level != "configural") {
      loading_labels[-1] <- sprintf(" (l%d_%d)", index, k[-1])
    }
    loadings <- c(
      sprintf("%s BY %s%s;", pair[[1]], x_items, loading_labels),
      sprintf("%s BY %s%s;", pair[[2]], y_items, loading_labels)
    )
    if (!level %in% c("partial", "scalar")) {
      return(loadings)
    }

    free <- if (level == "partial") design$free[[index]] else integer()
    intercept_labels <- ifelse(k %in% free, "", sprintf(" (i%d_%d)", index, k))
    c(
      loadings,
      sprintf("[%s]%s;", x_items, intercept_labels),
      sprintf("[%s]%s;", y_items, intercept_labels),
      sprintf("[%s@0];", pair[[1]]),
      sprintf("[%s];", pair[[2]])
    )
  }
  unlist(lapply(seq_along(design$pairs), one_pair))
}

# Fit of each level and the scaled chi-square difference test (Satorra &
# Bentler, 2001) against the less restricted level.
rsa_mplus_invariance_comparison <- function(models) {
  value <- function(workflow, name) {
    x <- workflow$results$summaries[[name]]
    if (length(x) == 0L) NA_real_ else as.numeric(x[[1]])
  }
  fits <- data.frame(
    level = names(models),
    k = vapply(models, value, numeric(1), "Parameters"),
    chisq = vapply(models, value, numeric(1), "ChiSqM_Value"),
    df = vapply(models, value, numeric(1), "ChiSqM_DF"),
    scaling = vapply(models, value, numeric(1), "ChiSqM_ScalingCorrection"),
    CFI = vapply(models, value, numeric(1), "CFI"),
    RMSEA = vapply(models, value, numeric(1), "RMSEA_Estimate"),
    SRMR = vapply(models, value, numeric(1), "SRMR"),
    against = NA_character_,
    delta_CFI = NA_real_,
    diff_chisq = NA_real_,
    diff_df = NA_real_,
    diff_p = NA_real_
  )
  rownames(fits) <- NULL

  parents <- c(metric = "configural", partial = "metric", scalar = "metric")
  for (level in intersect(names(parents), fits$level)) {
    i <- match(level, fits$level)
    j <- match(parents[[level]], fits$level)
    fits$against[[i]] <- parents[[level]]
    fits$delta_CFI[[i]] <- fits$CFI[[i]] - fits$CFI[[j]]
    test <- rsa_mplus_scaled_chisq_difference(
      chisq0 = fits$chisq[[i]],
      df0 = fits$df[[i]],
      c0 = fits$scaling[[i]],
      chisq1 = fits$chisq[[j]],
      df1 = fits$df[[j]],
      c1 = fits$scaling[[j]]
    )
    fits$diff_chisq[[i]] <- test$chisq
    fits$diff_df[[i]] <- test$df
    fits$diff_p[[i]] <- test$p
  }
  fits
}

# `chisq0` and `chisq1` are the scaled statistics Mplus prints for the
# restricted and the less restricted model; without scaling factors (ML) the
# test is the ordinary difference test.
rsa_mplus_scaled_chisq_difference <- function(
  chisq0,
  df0,
  c0,
  chisq1,
  df1,
  c1
) {
  if (is.na(c0)) c0 <- 1
  if (is.na(c1)) c1 <- 1
  df <- df0 - df1
  out <- list(chisq = NA_real_, df = df, p = NA_real_)
  if (anyNA(c(chisq0, chisq1, df)) || df <= 0) {
    return(out)
  }
  cd <- (df0 * c0 - df1 * c1) / df
  if (cd <= 0) {
    return(out)
  }
  out$chisq <- (chisq0 * c0 - chisq1 * c1) / cd
  out$p <- stats::pchisq(out$chisq, df = df, lower.tail = FALSE)
  out
}
