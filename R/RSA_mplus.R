#' Extract and plot RSA surfaces from Mplus output
#'
#' @param model Path to an Mplus `.out` file, an object returned by
#'   [MplusAutomation::readModels()], an `mplusObject` containing results, or
#'   an [rsa_mplus_workflow()] object.
#' @param outcome Dependent variable label from the Mplus regression table.
#'   Inferred for generated RSA-Mplus workflows.
#' @param pred_x Label of the linear X predictor in the Mplus output. Inferred
#'   for generated RSA-Mplus workflows.
#' @param pred_y Label of the linear Y predictor in the Mplus output. Inferred
#'   for generated RSA-Mplus workflows.
#' @param pred_x2 Label of the squared X term in the Mplus output. Inferred for
#'   generated RSA-Mplus workflows.
#' @param pred_xy Label of the XY interaction term in the Mplus output. Inferred
#'   for generated RSA-Mplus workflows.
#' @param pred_y2 Label of the squared Y term in the Mplus output. Inferred for
#'   generated RSA-Mplus workflows.
#' @param rsa_model Name of the fitted model in the `RSA` package's model
#'   family. Terms that this model sets to zero may be missing from the Mplus
#'   output, and their labels are optional. One of `"full"`, `"additive"`,
#'   `"IA"`, `"SQD"`, `"SRRR"`, `"SRR"`, `"RR"`, `"SSQD"`, `"SRSQD"`, `"null"`,
#'   `"onlyx2"`, `"onlyy2"`, `"onlyx"`, `"onlyy"`, `"mean"`, `"diff"`,
#'   `"weak"`, or `"strong"`. If `NULL` (the default), the model that
#'   [create_rsa_mplus_model()] generated for a workflow, and otherwise
#'   `"full"`. Also passed to [RSA::plotRSA()] as `model`. See Details.
#' @param b0 Optional intercept passed to [RSA::plotRSA()]. If `NULL`, the
#'   function looks for `<outcome><-Intercepts>` in the Mplus expectation
#'   parameters and otherwise falls back to `0`.
#' @param coef_type Which Mplus coefficient table to use. Passed to
#'   [stats::coef()] for `mplus.model` objects.
#' @param new_labels Optional character vector of `NEW` parameter labels to
#'   return. Matching is case-insensitive.
#' @param include_new If `TRUE`, include `NEW` parameters from `MODEL
#'   CONSTRAINT` in the returned object.
#' @param plot If `TRUE`, call [RSA::plotRSA()]. If `FALSE`, only return the
#'   extracted coefficients and metadata.
#' @param xlab Optional x-axis label passed to [RSA::plotRSA()]. Defaults to
#'   `pred_x`.
#' @param ylab Optional y-axis label passed to [RSA::plotRSA()]. Defaults to
#'   `pred_y`.
#' @param zlab Optional z-axis label passed to [RSA::plotRSA()]. Defaults to
#'   `outcome`.
#' @param ... Additional arguments passed to [RSA::plotRSA()].
#'
#' @return A list with the extracted RSA coefficients, selected Mplus
#'   parameters, optional `NEW` parameters, and the plot object when
#'   `plot = TRUE`. `coefficients` always contains all five polynomial terms
#'   and `b0`. `regression_parameters` contains only the terms found in the
#'   Mplus output, and `omitted_terms` names the terms that were absent and
#'   set to zero.
#' @details `RSA::plotRSA()` accepts the polynomial coefficients directly
#'   (`x`, `y`, `x2`, `xy`, `y2`, and optionally `b0`).
#'
#'   The `RSA` package compares a family of restricted second-order surfaces
#'   against the full polynomial. In Mplus a restriction is written in one of
#'   three ways: a term fixed at zero (`Z ON XS@0;`), shared labels with
#'   `MODEL CONSTRAINT` (as in the ridge models), or a term dropped from the
#'   model, which leaves no row in the output. Terms fixed at zero and
#'   constrained terms are read as usual. Dropped terms need `rsa_model`, which
#'   names the model and thereby the terms it sets to zero:
#'
#'   | `rsa_model` | Terms set to zero |
#'   |---|---|
#'   | `"full"`, `"SRRR"`, `"SRR"`, `"RR"`, `"SSQD"`, `"SRSQD"`, `"weak"`, `"strong"`, `"a5zero"` | none |
#'   | `"noIA"` | `xy` |
#'   | `"addx2"` | `xy`, `y2` |
#'   | `"addy2"` | `x2`, `xy` |
#'   | `"SQD"` | `x`, `y` |
#'   | `"IA"` | `x2`, `y2` |
#'   | `"additive"`, `"mean"`, `"diff"` | `x2`, `xy`, `y2` |
#'   | `"onlyx2"` | `y`, `xy`, `y2` |
#'   | `"onlyy2"` | `x`, `x2`, `xy` |
#'   | `"onlyx"` | `y`, `x2`, `xy`, `y2` |
#'   | `"onlyy"` | `x`, `x2`, `xy`, `y2` |
#'   | `"null"` | all five |
#'
#'   A term set to zero is taken as 0 when its label is `NULL` or its row is
#'   absent from the output. If the row is present, its estimate must be 0;
#'   a nonzero estimate means the output is not the named model and is an
#'   error. The remaining terms must be found, so a mistyped label still
#'   fails. One set of labels can therefore be reused across all variants.
#'   Equality constraints such as `b1 = b2` are not checked. The cubic and
#'   absolute-difference models of `RSA` are not polynomial surfaces of second
#'   order and are not supported.
#'
#'   `NEW` parameters from
#'   `MODEL CONSTRAINT` are returned for inspection. When the generated
#'   surface parameters (`CS`, `CC`, `IS`, `IC`, and `A5`) are available, the
#'   default three-dimensional plot annotation uses these Mplus estimates and
#'   adds significance markers. Surface parameters that `rsa_model` fixes at
#'   zero, such as `CC` in the ridge models, may be missing and are shown as 0
#'   without a marker. For frequentist estimators the markers are `*`,
#'   `**`, and `***` at p-values no greater than .05, .01, and .001,
#'   respectively. For `ESTIMATOR = BAYES` Mplus reports a 95% credibility
#'   interval next to a one-tailed posterior p-value; the credibility interval
#'   is the preferred decision rule, so a single `*` marks parameters whose
#'   interval excludes zero and the posterior p-value is ignored. The
#'   credibility interval bounds are also added to `new_parameters` as
#'   `lower_2.5ci` and `upper_2.5ci`. Other plot types and incomplete output
#'   retain the annotation behavior of `RSA::plotRSA()`.
#' @export
#'
#' @examples
#' \dontrun{
#' rsa_pos <- RSA_mplus(
#'   system.file("extdata", "congruence_sim.out", package = "rsahelpers"),
#'   outcome = "Z",
#'   pred_x = "X",
#'   pred_y = "Y",
#'   pred_x2 = "XS",
#'   pred_xy = "XY",
#'   pred_y2 = "YS",
#'   plot = TRUE
#' )
#'
#' # An additive model without quadratic terms in the Mplus syntax
#' additive <- RSA_mplus(
#'   "additive.out",
#'   outcome = "Z",
#'   pred_x = "X",
#'   pred_y = "Y",
#'   rsa_model = "additive"
#' )
#' }
RSA_mplus <- function(
  model,
  outcome = NULL,
  pred_x = NULL,
  pred_y = NULL,
  pred_x2 = NULL,
  pred_xy = NULL,
  pred_y2 = NULL,
  rsa_model = NULL,
  b0 = NULL,
  coef_type = c("un", "std", "stdy", "stdyx"),
  new_labels = NULL,
  include_new = TRUE,
  plot = TRUE,
  xlab = NULL,
  ylab = NULL,
  zlab = NULL,
  ...
) {
  workflow <- NULL
  metadata <- NULL

  if (inherits(model, "rsa_mplus_workflow")) {
    workflow <- model
    metadata <- workflow$spec$metadata
    if (!is.null(workflow$results)) {
      model <- workflow$results
    } else if (!is.null(workflow$mplus) && !is.null(workflow$mplus$results)) {
      model <- workflow$mplus
    } else {
      rlang::abort("The RSA-Mplus workflow does not contain fitted results.")
    }
  }

  inferred <- list(
    outcome = outcome,
    pred_x = pred_x,
    pred_y = pred_y,
    pred_x2 = pred_x2,
    pred_xy = pred_xy,
    pred_y2 = pred_y2
  )
  for (argument in names(inferred)) {
    if (is.null(inferred[[argument]]) && !is.null(metadata[[argument]])) {
      inferred[[argument]] <- metadata[[argument]]
    }
  }

  if (is.null(rsa_model)) {
    rsa_model <- if (is.null(metadata$rsa_model)) "full" else metadata$rsa_model
  }
  rsa_model <- rlang::arg_match(rsa_model, names(rsa_mplus_zero_terms))
  zero_terms <- rsa_mplus_zero_terms[[rsa_model]]
  term_arguments <- c(
    x = "pred_x",
    y = "pred_y",
    x2 = "pred_x2",
    xy = "pred_xy",
    y2 = "pred_y2"
  )

  required_arguments <- c(
    "outcome",
    unname(term_arguments[setdiff(names(term_arguments), zero_terms)])
  )
  missing_arguments <- required_arguments[
    vapply(inferred[required_arguments], is.null, logical(1))
  ]
  if (length(missing_arguments) > 0L) {
    error_message <- sprintf(
      "Missing required argument%s: %s.",
      if (length(missing_arguments) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_arguments), collapse = ", ")
    )
    if (rsa_model == "full" && !identical(missing_arguments, "outcome")) {
      error_message <- c(
        error_message,
        i = "For a surface without some terms, name it with `rsa_model`, e.g. `rsa_model = \"additive\"`."
      )
    }
    rlang::abort(error_message)
  }
  outcome <- inferred$outcome
  pred_x <- inferred$pred_x
  pred_y <- inferred$pred_y
  pred_x2 <- inferred$pred_x2
  pred_xy <- inferred$pred_xy
  pred_y2 <- inferred$pred_y2

  if (is.null(b0) && !is.null(metadata$b0)) {
    b0 <- metadata$b0
  }

  normalize_label <- function(x) {
    toupper(trimws(x))
  }

  match_rows <- function(parameters, labels, context) {
    parameters <- parameters
    parameters$Label <- normalize_label(parameters$Label)

    if (length(labels) == 0L) {
      matched_rows <- parameters[0L, , drop = FALSE]
      matched_rows$term <- character()
      return(matched_rows)
    }

    matched_rows <- lapply(labels, function(label) {
      hits <- which(parameters$Label == normalize_label(label))

      if (length(hits) == 0L) {
        rlang::abort(sprintf(
          "Could not find `%s` in the Mplus %s.",
          label,
          context
        ))
      }

      if (length(hits) > 1L) {
        rlang::abort(sprintf(
          "Found multiple matches for `%s` in the Mplus %s.",
          label,
          context
        ))
      }

      parameters[hits, , drop = FALSE]
    })

    matched_rows <- do.call(rbind, matched_rows)
    matched_rows$term <- names(labels)
    rownames(matched_rows) <- NULL
    matched_rows
  }

  extract_intercept <- function(model, outcome, coef_type) {
    expectation_parameters <- tryCatch(
      as.data.frame(
        stats::coef(model, type = coef_type, params = "expectation"),
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(expectation_parameters) || nrow(expectation_parameters) == 0L) {
      return(NA_real_)
    }

    expectation_parameters$Label <- normalize_label(
      expectation_parameters$Label
    )
    intercept_label <- paste0(outcome, "<-Intercepts")
    hits <- which(
      expectation_parameters$Label == normalize_label(intercept_label)
    )

    if (length(hits) == 0L) {
      return(NA_real_)
    }

    expectation_parameters$est[hits[[1]]]
  }

  extract_new <- function(model, coef_type, new_labels = NULL) {
    new_parameters <- tryCatch(
      as.data.frame(
        stats::coef(model, type = coef_type, params = "new"),
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(new_parameters) || nrow(new_parameters) == 0L) {
      return(NULL)
    }

    new_parameters$Label <- normalize_label(new_parameters$Label)
    new_parameters <- rsa_mplus_add_credibility_intervals(
      parameters = new_parameters,
      model = model,
      coef_type = coef_type
    )

    if (is.null(new_labels)) {
      rownames(new_parameters) <- NULL
      return(new_parameters)
    }

    match_rows(
      parameters = new_parameters,
      labels = stats::setNames(new_labels, new_labels),
      context = "NEW parameters"
    )
  }

  if (inherits(model, "mplus.model") || inherits(model, "mplusObject")) {
    model <- model
  } else {
    if (!is.character(model) || length(model) != 1L || is.na(model)) {
      rlang::abort(
        "`model` must be a file path, Mplus model object, or RSA-Mplus workflow."
      )
    }

    if (!requireNamespace("MplusAutomation", quietly = TRUE)) {
      rlang::abort(
        "Package `MplusAutomation` must be installed to read Mplus output files."
      )
    }

    model <- MplusAutomation::readModels(target = model, quiet = TRUE)
  }

  coef_type <- match.arg(coef_type)

  # MplusAutomation errors instead of returning an empty table when a model
  # has no regressions, as in the `null` model.
  regression_parameters <- tryCatch(
    as.data.frame(
      stats::coef(model, type = coef_type, params = "regression"),
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(regression_parameters) || nrow(regression_parameters) == 0L) {
    if (length(zero_terms) < length(term_arguments)) {
      rlang::abort("No regression parameters were found in the Mplus output.")
    }
    regression_parameters <- data.frame(
      Label = character(),
      est = numeric(),
      se = numeric(),
      pval = numeric()
    )
  }

  term_labels <- c(
    x = pred_x,
    y = pred_y,
    x2 = pred_x2,
    xy = pred_xy,
    y2 = pred_y2
  )
  requested_labels <- stats::setNames(
    paste0(outcome, "<-", term_labels, recycle0 = TRUE),
    names(term_labels)
  )

  # Terms the model sets to zero are matched only when their row exists, i.e.
  # when Mplus fixed them (`@0`) rather than dropping them from the model.
  found <- normalize_label(requested_labels) %in%
    normalize_label(regression_parameters$Label)
  matched_terms <- names(requested_labels)[
    found | !names(requested_labels) %in% zero_terms
  ]

  matched_regressions <- match_rows(
    regression_parameters,
    requested_labels[matched_terms],
    context = "regression parameters"
  )

  estimated_zero <- matched_regressions$term %in% zero_terms &
    !matched_regressions$est %in% 0
  if (any(estimated_zero)) {
    rlang::abort(c(
      sprintf(
        "`rsa_model = \"%s\"` sets %s to zero, but the Mplus output has nonzero estimates:",
        rsa_model,
        paste(sprintf("`%s`", zero_terms), collapse = ", ")
      ),
      stats::setNames(
        sprintf(
          "`%s` = %s",
          matched_regressions$Label[estimated_zero],
          format(matched_regressions$est[estimated_zero])
        ),
        rep("x", sum(estimated_zero))
      )
    ))
  }

  rsa_terms <- names(term_arguments)
  omitted_terms <- setdiff(rsa_terms, matched_regressions$term)
  coefficients <- stats::setNames(numeric(length(rsa_terms)), rsa_terms)
  coefficients[matched_regressions$term] <- matched_regressions$est

  if (is.null(b0)) {
    b0 <- extract_intercept(model, outcome = outcome, coef_type = coef_type)

    if (is.na(b0)) {
      b0 <- 0
      warning(
        sprintf(
          "No intercept for outcome `%s` was found in the Mplus expectation parameters; using `b0 = 0`.",
          outcome
        ),
        call. = FALSE
      )
    }
  }

  coefficients <- c(coefficients, b0 = b0)

  new_parameters <- NULL
  if (isTRUE(include_new)) {
    new_parameters <- extract_new(
      model = model,
      coef_type = coef_type,
      new_labels = new_labels
    )
  }

  plot_object <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("RSA", quietly = TRUE)) {
      rlang::abort("Package `RSA` must be installed when `plot = TRUE`.")
    }

    if (is.null(xlab)) {
      xlab <- pred_x
    }
    if (is.null(ylab)) {
      ylab <- pred_y
    }
    if (is.null(zlab)) {
      zlab <- outcome
    }

    plot_args <- c(
      list(
        x = unname(coefficients["x"]),
        y = unname(coefficients["y"]),
        x2 = unname(coefficients["x2"]),
        y2 = unname(coefficients["y2"]),
        xy = unname(coefficients["xy"]),
        b0 = unname(coefficients["b0"]),
        model = rsa_model,
        xlab = xlab,
        ylab = ylab,
        zlab = zlab
      ),
      list(...)
    )

    plot_object <- do.call(RSA::plotRSA, plot_args)

    parameter_annotation <- rsa_mplus_parameter_annotation(
      new_parameters,
      zero = rsa_mplus_zero_surface(rsa_model)
    )
    if (
      !is.null(parameter_annotation) &&
        inherits(plot_object, "trellis") &&
        !is.null(plot_object$panel.args.common$SPs)
    ) {
      plot_object$panel.args.common$SPs <- parameter_annotation
    }
  }

  out <- list(
    plot = plot_object,
    coefficients = coefficients,
    regression_parameters = matched_regressions,
    rsa_model = rsa_model,
    omitted_terms = omitted_terms,
    new_parameters = new_parameters,
    model = model,
    workflow = workflow,
    outcome = outcome,
    coefficient_type = coef_type
  )

  class(out) <- "rsa_mplus"
  out
}

rsa_mplus_ci_columns <- c(lower = "lower_2.5ci", upper = "upper_2.5ci")

# Polynomial terms that each second-order model of `RSA::RSA()` fixes at zero.
# Models that only add equality constraints keep all five terms.
rsa_mplus_zero_terms <- list(
  full = character(),
  additive = c("x2", "xy", "y2"),
  IA = c("x2", "y2"),
  SQD = c("x", "y"),
  SRRR = character(),
  SRR = character(),
  RR = character(),
  SSQD = character(),
  SRSQD = character(),
  null = c("x", "y", "x2", "xy", "y2"),
  onlyx2 = c("y", "xy", "y2"),
  onlyy2 = c("x", "x2", "xy"),
  onlyx = c("y", "x2", "xy", "y2"),
  onlyy = c("x", "x2", "xy", "y2"),
  mean = c("x2", "xy", "y2"),
  diff = c("x2", "xy", "y2"),
  weak = character(),
  strong = character(),
  addx2 = c("xy", "y2"),
  addy2 = c("x2", "xy"),
  noIA = "xy",
  a5zero = character()
)

# Surface parameters that a generated model fixes at zero, so that
# create_rsa_mplus_model() does not define them. Models it does not generate,
# such as "weak", have all five.
rsa_mplus_zero_surface <- function(rsa_model) {
  if (!rsa_model %in% names(rsa_mplus_model_restrictions)) {
    return(character())
  }
  defined <- rsa_mplus_surface_definitions(rsa_mplus_model_spec(rsa_model))
  toupper(setdiff(names(rsa_mplus_surface_weights), names(defined)))
}

rsa_mplus_parameter_annotation <- function(parameters, zero = character()) {
  required_columns <- c("Label", "est")
  if (
    is.null(parameters) ||
      !all(required_columns %in% names(parameters)) ||
      nrow(parameters) == 0L
  ) {
    return(NULL)
  }

  labels <- c("CS", "CC", "IS", "IC", "A5")
  estimated <- setdiff(labels, zero)
  normalized <- toupper(trimws(parameters$Label))
  rows <- match(estimated, normalized)
  if (
    length(estimated) == 0L || anyNA(rows) || anyDuplicated(normalized[rows])
  ) {
    return(NULL)
  }

  surface_parameters <- parameters[rows, , drop = FALSE]
  if (anyNA(surface_parameters$est)) {
    return(NULL)
  }

  stars <- rsa_mplus_significance_stars(surface_parameters)
  if (is.null(stars)) {
    return(NULL)
  }

  # Parameters fixed at zero by the model are shown as 0 without a test.
  estimates <- stats::setNames(numeric(length(labels)), labels)
  estimates[estimated] <- surface_parameters$est
  markers <- stats::setNames(character(length(labels)), labels)
  markers[estimated] <- stars

  paste0(
    "a",
    seq_along(labels),
    ": ",
    formatC(unname(estimates), format = "f", digits = 2L),
    markers,
    collapse = "    "
  )
}

#' Significance markers for Mplus surface parameters
#'
#' Bayesian Mplus output reports a 95% credibility interval alongside a
#' one-tailed posterior p-value. The credibility interval is the preferred
#' decision rule, so a single `*` marks parameters whose interval excludes
#' zero and the posterior p-value is ignored. Frequentist output has no
#' interval in its estimates table and falls back to p-value thresholds.
#'
#' @param parameters Data frame of Mplus parameters.
#' @return Character vector of significance markers, or `NULL` when neither
#'   decision rule is available.
#' @noRd
rsa_mplus_significance_stars <- function(parameters) {
  lower <- parameters[[rsa_mplus_ci_columns[["lower"]]]]
  upper <- parameters[[rsa_mplus_ci_columns[["upper"]]]]

  if (!is.null(lower) && !is.null(upper) && !anyNA(lower) && !anyNA(upper)) {
    return(ifelse(lower > 0 | upper < 0, "*", ""))
  }

  p_values <- parameters$pval
  if (is.null(p_values) || anyNA(p_values)) {
    return(NULL)
  }

  stars <- rep("", length(p_values))
  stars[p_values <= 0.05] <- "*"
  stars[p_values <= 0.01] <- "**"
  stars[p_values <= 0.001] <- "***"
  stars
}

#' Raw Mplus parameter table for a coefficient type
#'
#' @param model An `mplus.model` or `mplusObject`.
#' @param coef_type One of `"un"`, `"std"`, `"stdy"`, or `"stdyx"`.
#' @return Data frame of parameters, or `NULL` when unavailable.
#' @noRd
rsa_mplus_raw_parameters <- function(model, coef_type) {
  table_name <- switch(
    coef_type,
    un = "unstandardized",
    std = "std.standardized",
    stdy = "stdy.standardized",
    stdyx = "stdyx.standardized",
    NULL
  )
  if (is.null(table_name)) {
    return(NULL)
  }

  parameters <- tryCatch(
    model$parameters[[table_name]],
    error = function(e) NULL
  )
  if (is.null(parameters) || nrow(parameters) == 0L) {
    return(NULL)
  }

  as.data.frame(parameters, stringsAsFactors = FALSE)
}

#' Attach Bayesian credibility intervals to `NEW` parameters
#'
#' [MplusAutomation::readModels()] keeps the 95% credibility interval bounds of
#' a Bayesian model in its estimates table, but [stats::coef()] drops them.
#' This helper joins them back on by parameter label. Frequentist output has no
#' such columns and is returned unchanged.
#'
#' @param parameters Data frame of `NEW` parameters with a `Label` column.
#' @param model An `mplus.model` or `mplusObject`.
#' @param coef_type One of `"un"`, `"std"`, `"stdy"`, or `"stdyx"`.
#' @return `parameters`, with credibility interval columns when available.
#' @noRd
rsa_mplus_add_credibility_intervals <- function(parameters, model, coef_type) {
  raw_parameters <- rsa_mplus_raw_parameters(model, coef_type)
  if (
    is.null(raw_parameters) ||
      !all(
        c("paramHeader", "param", rsa_mplus_ci_columns) %in%
          names(raw_parameters)
      )
  ) {
    return(parameters)
  }

  raw_parameters <- raw_parameters[
    raw_parameters$paramHeader == "New.Additional.Parameters",
    ,
    drop = FALSE
  ]
  if (nrow(raw_parameters) == 0L) {
    return(parameters)
  }

  rows <- match(
    toupper(trimws(parameters$Label)),
    toupper(trimws(raw_parameters$param))
  )
  for (column in rsa_mplus_ci_columns) {
    parameters[[column]] <- raw_parameters[[column]][rows]
  }

  parameters
}
