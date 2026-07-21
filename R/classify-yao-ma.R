#' Classify a response surface using the Yao-Ma typology
#'
#' @param x An [`RSA_mplus()`] result, an `rsa_mplus_workflow`, an
#'   `mplus.model`, an `mplusObject`, or a data frame containing response-surface
#'   parameters. A data frame must have columns `parameter`, `estimate`,
#'   `conf.low`, and `conf.high`.
#' @param valence Whether higher (`"positive"`) or lower (`"negative"`) values
#'   of the outcome are desirable. This argument is required because it selects
#'   the relevant principal axis.
#' @param equivalence Optional named numeric vector of practical-equivalence
#'   margins. Supply one nonnegative margin for each of `CS`, `CC`, `IS`, `IC`,
#'   `axis_intercept`, and `axis_slope`. When omitted, confidence intervals that
#'   include the null are treated as absence of an effect. When supplied, each
#'   condition is classified as equivalent, different, or indeterminate.
#' @param level Confidence level used when extracting intervals from Mplus.
#'   Intervals already present in a parameter data frame are used unchanged.
#'
#' @return An object of class `yao_ma_classification`. It contains a one-row
#'   `summary`, the normalized `parameters`, and a condition-level `decisions`
#'   data frame.
#' @details
#' The function uses `CS`, `CC`, `IS`, and `IC`, corresponding to the usual
#' RSA parameters `a1` through `a4`. For positively valenced outcomes it uses
#' `P10` and `P11`; for negatively valenced outcomes it uses `P20` and `P21`.
#'
#' Principal-axis parameters are nonlinear functions of the polynomial
#' coefficients. Mplus confidence intervals for these `MODEL CONSTRAINT`
#' parameters use its delta-method standard errors. A missing or degenerate
#' principal axis therefore produces a partial result rather than an invented
#' Yao-Ma subtype.
#' @export
#'
#' @examples
#' params <- data.frame(
#'   parameter = c("CS", "CC", "IS", "IC", "P10", "P11"),
#'   estimate = c(0.4, 0.01, 0.02, -0.5, 0.03, 1.02),
#'   conf.low = c(0.2, -0.1, -0.1, -0.7, -0.1, 0.9),
#'   conf.high = c(0.6, 0.1, 0.1, -0.3, 0.1, 1.1)
#' )
#' classify_yao_ma(params, valence = "positive")
classify_yao_ma <- function(
  x,
  valence,
  equivalence = NULL,
  level = 0.95
) {
  valence <- match.arg(valence, c("positive", "negative"))
  validate_yao_ma_level(level)
  margins <- validate_yao_ma_margins(equivalence)
  parameters <- yao_ma_parameters(x, level = level)

  surface <- c("CS", "CC", "IS", "IC")
  missing_surface <- setdiff(surface, parameters$parameter)
  if (length(missing_surface) > 0L) {
    rlang::abort(sprintf(
      "Missing required surface parameter%s: %s.",
      if (length(missing_surface) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_surface), collapse = ", ")
    ))
  }
  validate_yao_ma_intervals(parameters, surface)

  axis_parameters <- if (valence == "positive") {
    c(axis_intercept = "P10", axis_slope = "P11")
  } else {
    c(axis_intercept = "P20", axis_slope = "P21")
  }

  method <- if (is.null(margins)) "significance" else "equivalence"
  surface_margins <- if (is.null(margins)) {
    stats::setNames(rep(NA_real_, length(surface)), surface)
  } else {
    margins[surface]
  }

  decisions <- do.call(
    rbind,
    lapply(surface, function(parameter) {
      yao_ma_decision(
        parameters,
        condition = parameter,
        parameter = parameter,
        null = 0,
        margin = surface_margins[[parameter]],
        method = method
      )
    })
  )

  axis_available <- all(axis_parameters %in% parameters$parameter)
  if (axis_available) {
    axis_rows <- parameters$parameter %in% unname(axis_parameters)
    axis_values <- parameters[
      axis_rows,
      c("estimate", "conf.low", "conf.high")
    ]
    axis_available <- all(stats::complete.cases(axis_values)) &&
      all(is.finite(as.matrix(axis_values))) &&
      all(axis_values$conf.low <= axis_values$conf.high)
  }

  if (axis_available) {
    axis_margins <- if (is.null(margins)) {
      c(axis_intercept = NA_real_, axis_slope = NA_real_)
    } else {
      margins[c("axis_intercept", "axis_slope")]
    }
    decisions <- rbind(
      decisions,
      yao_ma_decision(
        parameters,
        condition = "axis_intercept",
        parameter = axis_parameters[["axis_intercept"]],
        null = 0,
        margin = axis_margins[["axis_intercept"]],
        method = method
      ),
      yao_ma_decision(
        parameters,
        condition = "axis_slope",
        parameter = axis_parameters[["axis_slope"]],
        null = 1,
        margin = axis_margins[["axis_slope"]],
        method = method
      )
    )
  }

  result <- yao_ma_classification_summary(
    parameters = parameters,
    decisions = decisions,
    valence = valence,
    method = method,
    axis_available = axis_available
  )

  structure(
    list(
      summary = result,
      parameters = parameters,
      decisions = decisions
    ),
    class = "yao_ma_classification"
  )
}

#' @export
print.yao_ma_classification <- function(x, ...) {
  summary <- x$summary
  cat("<yao_ma_classification>\n")
  cat("Method:", summary$method, "\n")
  cat("Valence:", summary$valence, "\n")
  cat("Status:", summary$status, "\n")
  if (!is.na(summary$type)) {
    cat("Type:", summary$type, "-", summary$label, "\n")
  } else {
    cat("Result:", summary$label, "\n")
  }
  invisible(x)
}

validate_yao_ma_level <- function(level) {
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1
  ) {
    rlang::abort("`level` must be one number strictly between 0 and 1.")
  }
  invisible(level)
}

validate_yao_ma_margins <- function(equivalence) {
  if (is.null(equivalence)) {
    return(NULL)
  }
  required <- c("CS", "CC", "IS", "IC", "axis_intercept", "axis_slope")
  if (
    !is.numeric(equivalence) ||
      is.null(names(equivalence)) ||
      anyNA(equivalence) ||
      any(!is.finite(equivalence)) ||
      any(equivalence < 0)
  ) {
    rlang::abort(
      "`equivalence` must be a named vector of finite, nonnegative numbers."
    )
  }
  normalized <- names(equivalence)
  normalized[toupper(normalized) %in% c("CS", "CC", "IS", "IC")] <-
    toupper(normalized[toupper(normalized) %in% c("CS", "CC", "IS", "IC")])
  normalized[tolower(normalized) == "axis_intercept"] <- "axis_intercept"
  normalized[tolower(normalized) == "axis_slope"] <- "axis_slope"
  names(equivalence) <- normalized
  if (anyDuplicated(names(equivalence))) {
    rlang::abort("`equivalence` margin names must be unique.")
  }
  missing <- setdiff(required, names(equivalence))
  extra <- setdiff(names(equivalence), required)
  if (length(missing) > 0L || length(extra) > 0L) {
    details <- c(
      if (length(missing) > 0L) {
        paste0("missing ", paste(missing, collapse = ", "))
      },
      if (length(extra) > 0L) {
        paste0("unknown ", paste(extra, collapse = ", "))
      }
    )
    rlang::abort(paste0(
      "`equivalence` must name exactly CS, CC, IS, IC, axis_intercept, ",
      "and axis_slope (",
      paste(details, collapse = "; "),
      ")."
    ))
  }
  equivalence[required]
}

yao_ma_parameters <- function(x, level) {
  if (is.data.frame(x)) {
    return(normalize_yao_ma_table(x))
  }

  if (inherits(x, "rsa_mplus")) {
    x <- x$model
  } else if (inherits(x, "rsa_mplus_workflow")) {
    if (!is.null(x$results)) {
      x <- x$results
    } else if (!is.null(x$mplus) && !is.null(x$mplus$results)) {
      x <- x$mplus
    } else {
      rlang::abort("The RSA-Mplus workflow does not contain fitted results.")
    }
  }

  if (!inherits(x, "mplus.model") && !inherits(x, "mplusObject")) {
    rlang::abort(
      paste(
        "`x` must be an RSA-Mplus result, fitted Mplus object, or parameter",
        "data frame."
      )
    )
  }
  if (!level %in% c(0.90, 0.95, 0.99)) {
    rlang::abort(
      "Mplus confidence intervals support `level` values 0.90, 0.95, or 0.99."
    )
  }

  estimates <- tryCatch(
    as.data.frame(stats::coef(x, type = "un", params = "new")),
    error = function(e) NULL
  )
  intervals <- tryCatch(
    as.data.frame(stats::confint(
      x,
      level = level,
      type = "un",
      params = "new"
    )),
    error = function(e) NULL
  )
  if (is.null(estimates) || nrow(estimates) == 0L) {
    rlang::abort("No Mplus `NEW` parameters were found.")
  }
  if (is.null(intervals) || nrow(intervals) == 0L) {
    rlang::abort(
      "No Mplus confidence intervals were found; request `CINTERVAL` in OUTPUT."
    )
  }

  estimates$parameter <- toupper(trimws(estimates$Label))
  intervals$parameter <- toupper(trimws(intervals$Label))
  matched <- match(estimates$parameter, intervals$parameter)
  data.frame(
    parameter = estimates$parameter,
    estimate = estimates$est,
    conf.low = intervals$LowerCI[matched],
    conf.high = intervals$UpperCI[matched],
    std.error = if ("se" %in% names(estimates)) estimates$se else NA_real_,
    p.value = if ("pval" %in% names(estimates)) estimates$pval else NA_real_,
    stringsAsFactors = FALSE
  ) |>
    normalize_yao_ma_table()
}

normalize_yao_ma_table <- function(x) {
  required <- c("parameter", "estimate", "conf.low", "conf.high")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    rlang::abort(sprintf(
      "Parameter data is missing column%s: %s.",
      if (length(missing) == 1L) "" else "s",
      paste(sprintf("`%s`", missing), collapse = ", ")
    ))
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE)
  out$parameter <- toupper(trimws(as.character(out$parameter)))
  numeric_columns <- intersect(
    c("estimate", "conf.low", "conf.high", "std.error", "p.value"),
    names(out)
  )
  if (any(!vapply(out[numeric_columns], is.numeric, logical(1)))) {
    rlang::abort("Estimate and interval columns must be numeric.")
  }
  if (anyNA(out$parameter) || any(!nzchar(out$parameter))) {
    rlang::abort("Parameter names must be non-missing and non-empty.")
  }
  if (anyDuplicated(out$parameter)) {
    duplicates <- unique(out$parameter[duplicated(out$parameter)])
    rlang::abort(sprintf(
      "Parameter names must be unique; duplicated: %s.",
      paste(duplicates, collapse = ", ")
    ))
  }
  out
}

validate_yao_ma_intervals <- function(parameters, required) {
  rows <- match(required, parameters$parameter)
  values <- parameters[rows, c("estimate", "conf.low", "conf.high")]
  if (
    any(!stats::complete.cases(values)) || any(!is.finite(as.matrix(values)))
  ) {
    rlang::abort("Required surface estimates and intervals must be finite.")
  }
  if (any(values$conf.low > values$conf.high)) {
    rlang::abort(
      "Each confidence-interval lower bound must not exceed its upper bound."
    )
  }
  invisible(parameters)
}

yao_ma_decision <- function(
  parameters,
  condition,
  parameter,
  null,
  margin,
  method
) {
  row <- parameters[parameters$parameter == parameter, , drop = FALSE]
  if (nrow(row) != 1L) {
    rlang::abort(sprintf(
      "Could not uniquely identify parameter `%s`.",
      parameter
    ))
  }
  low <- row$conf.low
  high <- row$conf.high
  state <- if (method == "significance") {
    if (low > null || high < null) "different" else "not_different"
  } else if (low >= null - margin && high <= null + margin) {
    "equivalent"
  } else if (low > null + margin || high < null - margin) {
    "different"
  } else {
    "indeterminate"
  }
  direction <- if (high < null) {
    "below"
  } else if (low > null) {
    "above"
  } else {
    "overlaps"
  }
  data.frame(
    condition = condition,
    parameter = parameter,
    estimate = row$estimate,
    conf.low = low,
    conf.high = high,
    null = null,
    margin = margin,
    state = state,
    direction = direction,
    stringsAsFactors = FALSE
  )
}

yao_ma_classification_summary <- function(
  parameters,
  decisions,
  valence,
  method,
  axis_available
) {
  decision <- function(condition) {
    decisions[decisions$condition == condition, , drop = FALSE]
  }
  is_different <- function(condition) decision(condition)$state == "different"

  ic <- decision("IC")
  desired_direction <- if (valence == "positive") "below" else "above"
  ic_desired <- is_different("IC") && ic$direction == desired_direction
  ic_opposite <- is_different("IC") && ic$direction != desired_direction

  if (method == "equivalence" && ic$state == "indeterminate") {
    return(yao_ma_summary_row(
      status = "indeterminate",
      label = "indeterminate congruence curvature",
      valence = valence,
      method = method
    ))
  }

  if (!ic_desired) {
    cs_different <- is_different("CS")
    is_state <- decision("IS")$state
    is_different_value <- identical(is_state, "different")
    label <- if (ic_opposite) {
      "outside typology: incongruence curvature has the opposite sign"
    } else if (is_different_value && !cs_different) {
      "outside typology: directional difference effect without congruence"
    } else if (cs_different && !is_different_value) {
      "outside typology: level effect without congruence"
    } else if (cs_different || is_different_value) {
      "outside typology: level and difference effects without congruence"
    } else {
      "no congruence effect"
    }
    return(yao_ma_summary_row(
      status = "outside_typology",
      label = label,
      valence = valence,
      method = method
    ))
  }

  fit_conditions <- c("CS", "CC")
  if (
    method == "equivalence" &&
      any(
        decisions$state[match(fit_conditions, decisions$condition)] ==
          "indeterminate"
      )
  ) {
    return(yao_ma_summary_row(
      status = "indeterminate",
      label = "indeterminate level-effect condition",
      valence = valence,
      method = method
    ))
  }

  if (!axis_available) {
    return(yao_ma_summary_row(
      status = "partial",
      label = "congruence detected; principal-axis subtype unavailable",
      valence = valence,
      method = method
    ))
  }

  axis_conditions <- c("axis_intercept", "axis_slope")
  if (
    method == "equivalence" &&
      any(
        decisions$state[match(axis_conditions, decisions$condition)] ==
          "indeterminate"
      )
  ) {
    return(yao_ma_summary_row(
      status = "indeterminate",
      label = "indeterminate principal-axis condition",
      valence = valence,
      method = method
    ))
  }

  intercept_different <- is_different("axis_intercept")
  slope_different <- is_different("axis_slope")
  cs_different <- is_different("CS")
  cc_different <- is_different("CC")

  base <- if (!slope_different && !intercept_different) {
    1L
  } else if (!slope_different && intercept_different) {
    5L
  } else if (slope_different && !intercept_different) {
    9L
  } else {
    13L
  }
  offset <- if (!cs_different && !cc_different) {
    0L
  } else if (cs_different && !cc_different) {
    1L
  } else if (!cs_different && cc_different) {
    2L
  } else {
    3L
  }
  type <- base + offset
  fit_line <- c("", " & LLE", " & SLE", " & CLE")[[offset + 1L]]
  archetype <- if (base %in% c(1L, 9L)) {
    "Exact correspondence"
  } else {
    "Commensurate compatibility"
  }
  contingency <- if (base >= 9L) " with contingency" else ""

  direction <- if (!intercept_different && !slope_different) {
    "none"
  } else if (slope_different) {
    "level_dependent"
  } else {
    intercept <- decision("axis_intercept")$estimate
    if (intercept > 0) "y_greater_than_x" else "x_greater_than_y"
  }

  yao_ma_summary_row(
    status = "classified",
    type = type,
    label = paste0(archetype, fit_line, contingency),
    in_typology = TRUE,
    direction = direction,
    valence = valence,
    method = method
  )
}

yao_ma_summary_row <- function(
  status,
  label,
  valence,
  method,
  type = NA_integer_,
  in_typology = FALSE,
  direction = NA_character_
) {
  data.frame(
    type = type,
    label = label,
    in_typology = in_typology,
    status = status,
    direction = direction,
    valence = valence,
    method = method,
    stringsAsFactors = FALSE
  )
}
