# Edwards-Parry congruence spline regression helpers.

# This file is intentionally written as a package-ready R script. Public
# functions have roxygen2 documentation blocks; internal helpers use ordinary
# comments and should remain unexported if this code is moved into a package.
#
# The nonlinear models follow Edwards and Parry (2018). The response is `z`,
# the wanted/ideal component is `x`, and the actual/current component is `y`.
#
# One seam:
#   z = b0 + b1*x + b2*y + b3*(y - c0 - c1*x) * I(y < c0 + c1*x) + e
#
# Two seams:
#   z = b0 + b1*x + b2*y
#       + b3*(y - c10 - c11*x) * I(y < c10 + c11*x)
#       + b4*(y - c20 - c21*x) * I(y < c20 + c21*x) + e

# Normalize a variable supplied as a string or symbol in internal calls.
as_variable_name <- function(x, expr, arg_name) {
  if (is.character(x) && length(x) == 1L) {
    return(x)
  }

  if (is.symbol(expr)) {
    return(as.character(expr))
  }

  stop(arg_name, " must be a single column name or unquoted column.", call. = FALSE)
}

# Normalize a variable supplied to a public function. The `tryCatch()` lets
# unquoted column names be captured without being evaluated in the caller.
resolve_public_variable <- function(value, expr, arg_name) {
  value <- tryCatch(value, error = function(e) NULL)
  if (is.character(value) && length(value) == 1L) {
    return(value)
  }
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  stop(arg_name, " must be a single column name or unquoted column.", call. = FALSE)
}

# Moore-Penrose inverse for Wald tests when constraint Jacobians are singular
# or nearly singular. This mirrors Stata's tendency to drop weak constraints.
pinv <- function(x, tol = sqrt(.Machine$double.eps)) {
  sv <- svd(x)
  keep <- sv$d > tol * max(sv$d)
  if (!any(keep)) {
    return(matrix(NA_real_, nrow = ncol(x), ncol = nrow(x)))
  }
  sv$v[, keep, drop = FALSE] %*%
    (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
}

# Normalize the `center`/`scale` arguments. Logical values are accepted for
# backward compatibility: `center = TRUE` historically meant variablewise
# centering, `center = FALSE` meant no centering. The same convention applies to
# `scale` (TRUE was never pooled in older code, so it maps to "none").
resolve_center_method <- function(center) {
  if (is.logical(center)) {
    return(if (isTRUE(center)) "variablewise" else "none")
  }
  match.arg(center, c("pooled", "variablewise", "none"))
}

resolve_scale_method <- function(scale) {
  if (is.logical(scale)) {
    return(if (isTRUE(scale)) "pooled" else "none")
  }
  match.arg(scale, c("pooled", "none"))
}

# Build the complete-case model frame used by all public fitting functions.
#
# Centering and scaling follow the `RSA` package conventions because congruence
# interpretation depends on `x` and `y` sharing one origin and one scale.
# Variablewise centering (the historical default) silently distorts the `X = Y`
# line, so the safe pooled options are now the defaults.
coerce_complete_model_frame <- function(data, wanted, actual, outcome,
                                        center = c("pooled", "variablewise", "none"),
                                        scale = c("pooled", "none")) {
  wanted <- as_variable_name(wanted, substitute(wanted), "wanted")
  actual <- as_variable_name(actual, substitute(actual), "actual")
  outcome <- as_variable_name(outcome, substitute(outcome), "outcome")
  center <- resolve_center_method(center)
  scale <- resolve_scale_method(scale)

  missing_cols <- setdiff(c(actual, wanted, outcome), names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  x_raw <- as.numeric(data[[wanted]])
  y_raw <- as.numeric(data[[actual]])
  z <- as.numeric(data[[outcome]])
  keep <- stats::complete.cases(x_raw, y_raw, z)

  x <- x_raw[keep]
  y <- y_raw[keep]

  x_center <- 0
  y_center <- 0
  if (center == "variablewise") {
    x_center <- mean(x)
    y_center <- mean(y)
  } else if (center == "pooled") {
    x_center <- y_center <- mean(c(x, y))
  }
  x <- x - x_center
  y <- y - y_center

  scale_value <- 1
  if (scale == "pooled") {
    scale_value <- stats::sd(c(x, y))
  }
  x <- x / scale_value
  y <- y / scale_value

  out <- data.frame(
    x = x,
    y = y,
    z = z[keep],
    .row = which(keep),
    check.names = FALSE
  )
  attr(out, "x_center") <- x_center
  attr(out, "y_center") <- y_center
  attr(out, "scale") <- scale_value
  attr(out, "center_method") <- center
  attr(out, "scale_method") <- scale
  out
}

#' Prepare variables for congruence regression models
#'
#' Create the centered component variables and auxiliary terms used by
#' Edwards-Parry absolute-difference, piecewise, and spline congruence models.
#'
#' @param data A data frame containing the component and outcome variables.
#' @param wanted,actual Column names, quoted or unquoted, for the two component
#'   variables. By field convention, `wanted` or ideal values are mapped to
#'   `x`; `actual` or current values are mapped to `y`.
#' @param outcome Column name, quoted or unquoted, for the response variable.
#' @param center Centering of the component variables, following the `RSA`
#'   package. One of `"pooled"` (subtract the shared mean of `x` and `y`;
#'   default), `"variablewise"` (subtract each component's own mean), or
#'   `"none"`. Logical values are accepted for backward compatibility:
#'   `TRUE` maps to `"variablewise"`, `FALSE` to `"none"`. Pooled centering
#'   preserves the `X = Y` congruence interpretation; use `"none"` for columns
#'   that are already centered, such as the downloaded `*HC` and `*WC` columns.
#' @param scale Scaling of the component variables. One of `"pooled"` (divide
#'   both `x` and `y` by their pooled standard deviation; default) or `"none"`.
#' @param hinge_offset Numeric half-width, in **working-scale** units, of the
#'   fixed two-seam hinges `hinge_upper`/`hinge_lower` (`Y = X ± hinge_offset`).
#'   Because the offset is applied after centering and scaling, under pooled
#'   `scale = "pooled"` an offset of `1` corresponds to one pooled standard
#'   deviation rather than one raw scale unit.
#'
#' @return A data frame with standardized columns `x`, `y`, `z`, `.row`, and
#'   helper columns for absolute difference, one-break piecewise regression,
#'   constrained piecewise regression, and fixed two-seam comparisons. The
#'   centering and scaling constants are recorded as attributes (`x_center`,
#'   `y_center`, `scale`, `center_method`, `scale_method`).
#'
#' @examples
#' \dontrun{
#' prepared <- prepare_congruence_data(my_data, wanted, actual, satisfaction)
#' }
#'
#' @export
prepare_congruence_data <- function(data, wanted, actual, outcome,
                                    center = c("pooled", "variablewise", "none"),
                                    scale = c("pooled", "none"),
                                    hinge_offset = 1) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  if (is.logical(center) && missing(scale)) scale <- "none"
  mf <- coerce_complete_model_frame(data, wanted, actual, outcome, center, scale)
  diff <- mf$y - mf$x
  w <- as.numeric(diff < 0)

  mf$abs_diff <- abs(diff)
  mf$w <- w
  mf$xw <- mf$x * w
  mf$yw <- mf$y * w
  mf$zd <- diff * w
  mf$hinge_upper <- spline_hinge(mf$x, mf$y, hinge_offset, 1)
  mf$hinge_lower <- spline_hinge(mf$x, mf$y, -hinge_offset, 1)
  attr(mf, "hinge_offset") <- hinge_offset
  mf
}

# Compute a continuous hinge term for a seam line y = c0 + c1 * x.
spline_hinge <- function(x, y, c0, c1) {
  (y - c0 - c1 * x) * as.numeric(y < c0 + c1 * x)
}

# Generate fitted values for one- or two-seam nonlinear spline surfaces.
spline_prediction <- function(par, x, y, n_seams = 1L) {
  if (n_seams == 1L) {
    return(par[["b0"]] + par[["b1"]] * x + par[["b2"]] * y +
      par[["b3"]] * spline_hinge(x, y, par[["c0"]], par[["c1"]]))
  }

  par[["b0"]] + par[["b1"]] * x + par[["b2"]] * y +
    par[["b3"]] * spline_hinge(x, y, par[["c10"]], par[["c11"]]) +
    par[["b4"]] * spline_hinge(x, y, par[["c20"]], par[["c21"]])
}

# Canonical full parameter names for a one- or two-seam spline.
spline_param_names <- function(n_seams) {
  if (n_seams == 1L) {
    c("b0", "b1", "b2", "b3", "c0", "c1")
  } else {
    c("b0", "b1", "b2", "b3", "b4", "c10", "c11", "c20", "c21")
  }
}

# Validate and normalize a `fix` argument (named numeric of held-constant
# parameters) against the canonical parameter set for `n_seams`.
normalize_fix <- function(fix, n_seams) {
  if (is.null(fix)) {
    return(NULL)
  }
  if (!is.numeric(fix) || is.null(names(fix))) {
    stop("`fix` must be a named numeric vector of parameters to hold constant.", call. = FALSE)
  }
  template <- spline_param_names(n_seams)
  bad <- setdiff(names(fix), template)
  if (length(bad) > 0L) {
    stop(
      "`fix` names not in the model: ", paste(bad, collapse = ", "),
      ". Valid names: ", paste(template, collapse = ", "), ".",
      call. = FALSE
    )
  }
  fix[intersect(template, names(fix))]
}

# Residual vector for nonlinear least-squares fitting.
residual_vector <- function(par, mf, n_seams) {
  mf$z - spline_prediction(par, mf$x, mf$y, n_seams)
}

# Residual sum of squares for nonlinear least-squares fitting.
rss_for <- function(par, mf, n_seams) {
  r <- residual_vector(par, mf, n_seams)
  sum(r^2)
}

# Approximate the parameter covariance matrix from the numerical Jacobian.
#
# When `fix` holds parameters constant, the Jacobian is taken only over the free
# parameters and the fixed rows/columns are returned as zeros, so downstream
# delta-method machinery (which indexes the full parameter vector) sees those
# parameters as having no sampling variability.
vcov_from_jacobian <- function(par, mf, n_seams, fix = NULL) {
  template <- names(par)
  free <- setdiff(template, names(fix))
  out <- matrix(0, length(template), length(template), dimnames = list(template, template))

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    out[] <- NA_real_
    return(out)
  }

  pred_fun <- function(p_free) {
    p <- par
    p[free] <- p_free
    spline_prediction(p, mf$x, mf$y, n_seams)
  }
  jac <- numDeriv::jacobian(pred_fun, unname(par[free]))
  rss <- rss_for(par, mf, n_seams)
  df_resid <- nrow(mf) - length(free)
  xtx_inv <- pinv(crossprod(jac))
  out[free, free] <- (rss / df_resid) * xtx_inv
  out
}

# Derive Stata-style starting values from constrained or unconstrained
# piecewise OLS models.
start_from_piecewise <- function(mf, n_seams = 1L, method = c("constrained", "unconstrained"),
                                 hinge_offset = 1) {
  method <- match.arg(method)

  if (n_seams == 1L && method == "constrained") {
    fit <- stats::lm(z ~ x + y + zd, data = prepare_hinges_for_lm(mf, hinge_offset))
    cf <- stats::coef(fit)
    return(c(
      b0 = unname(cf[["(Intercept)"]]),
      b1 = unname(cf[["x"]]),
      b2 = unname(cf[["y"]]),
      b3 = unname(cf[["zd"]]),
      c0 = 0,
      c1 = 1
    ))
  }

  if (n_seams == 1L && method == "unconstrained") {
    fit <- stats::lm(z ~ x + y + w + xw + yw, data = prepare_hinges_for_lm(mf, hinge_offset))
    cf <- stats::coef(fit)
    a3 <- unname(cf[["w"]])
    a4 <- unname(cf[["xw"]])
    a5 <- unname(cf[["yw"]])
    c0 <- if (is.finite(a5) && abs(a5) > .Machine$double.eps) -a3 / a5 else 0
    c1 <- if (is.finite(a5) && abs(a5) > .Machine$double.eps) -a4 / a5 else 1
    return(c(
      b0 = unname(cf[["(Intercept)"]]),
      b1 = unname(cf[["x"]]),
      b2 = unname(cf[["y"]]),
      b3 = a5,
      c0 = c0,
      c1 = c1
    ))
  }

  if (n_seams == 2L && method == "constrained") {
    fit <- stats::lm(z ~ x + y + h1_fixed + h2_fixed, data = prepare_hinges_for_lm(mf, hinge_offset))
    cf <- stats::coef(fit)
    return(c(
      b0 = unname(cf[["(Intercept)"]]),
      b1 = unname(cf[["x"]]),
      b2 = unname(cf[["y"]]),
      b3 = unname(cf[["h1_fixed"]]),
      b4 = unname(cf[["h2_fixed"]]),
      c10 = hinge_offset,
      c11 = 1,
      c20 = -hinge_offset,
      c21 = 1
    ))
  }

  one <- start_from_piecewise(mf, 1L, "unconstrained", hinge_offset)
  c(
    b0 = unname(one[["b0"]]),
    b1 = unname(one[["b1"]]),
    b2 = unname(one[["b2"]]),
    b3 = unname(one[["b3"]]) / 2,
    b4 = unname(one[["b3"]]) / 2,
    c10 = unname(one[["c0"]]),
    c11 = unname(one[["c1"]]),
    c20 = -hinge_offset,
    c21 = 1
  )
}

# Add fixed seam and one-break terms used to compute starting values.
prepare_hinges_for_lm <- function(mf, hinge_offset = 1) {
  out <- mf
  diff <- out$y - out$x
  out$w <- as.numeric(diff < 0)
  out$xw <- out$x * out$w
  out$yw <- out$y * out$w
  out$zd <- diff * out$w
  out$h1_fixed <- spline_hinge(out$x, out$y, hinge_offset, 1)
  out$h2_fixed <- spline_hinge(out$x, out$y, -hinge_offset, 1)
  out
}

# Build the candidate start matrix. Multiple starts are important because the
# two-seam surface can have local minima and label-switched seam solutions.
make_start_grid <- function(mf, n_seams, starts = NULL, multistart = TRUE, hinge_offset = 1,
                            fix = NULL) {
  base <- if (n_seams == 1L) {
    rbind(
      constrained = start_from_piecewise(mf, 1L, "constrained", hinge_offset),
      unconstrained = start_from_piecewise(mf, 1L, "unconstrained", hinge_offset)
    )
  } else {
    rbind(
      constrained = start_from_piecewise(mf, 2L, "constrained", hinge_offset),
      unconstrained = start_from_piecewise(mf, 2L, "unconstrained", hinge_offset)
    )
  }

  if (!is.null(starts)) {
    if (is.list(starts) && !is.numeric(starts)) {
      user <- do.call(rbind, lapply(starts, function(x) x[colnames(base)]))
    } else {
      user <- matrix(starts[colnames(base)], nrow = 1L)
      colnames(user) <- colnames(base)
    }
    base <- rbind(user = user, base)
  }

  if (isTRUE(multistart)) {
    scales <- pmax(abs(base[1L, ]), 0.25)
    jitters <- rbind(
      plus = base[1L, ] + 0.10 * scales,
      minus = base[1L, ] - 0.10 * scales
    )
    if (n_seams == 2L) {
      swapped <- base[1L, c("b0", "b1", "b2", "b4", "b3", "c20", "c21", "c10", "c11")]
      names(swapped) <- colnames(base)
      jitters <- rbind(jitters, swapped = swapped)
    }
    base <- rbind(base, jitters)
  }

  grid <- unique(as.data.frame(base))
  if (!is.null(fix)) {
    grid <- grid[, setdiff(colnames(grid), names(fix)), drop = FALSE]
  }
  grid
}

# Build the spline model formula for nlsLM, inlining any fixed parameters as
# numeric literals so the solver only estimates the free parameters.
spline_model_formula <- function(n_seams, fix = NULL) {
  rhs <- if (n_seams == 1L) {
    "b0 + b1 * x + b2 * y + b3 * (y - c0 - c1 * x) * (y < c0 + c1 * x)"
  } else {
    paste(
      "b0 + b1 * x + b2 * y +",
      "b3 * (y - c10 - c11 * x) * (y < c10 + c11 * x) +",
      "b4 * (y - c20 - c21 * x) * (y < c20 + c21 * x)"
    )
  }
  for (nm in names(fix)) {
    rhs <- gsub(
      paste0("\\b", nm, "\\b"),
      sprintf("(%.10g)", fix[[nm]]),
      rhs
    )
  }
  stats::as.formula(paste("z ~", rhs))
}

# Primary nonlinear solver using Levenberg-Marquardt via minpack.lm.
fit_with_nlslm <- function(start, mf, n_seams, control, fix = NULL) {
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("minpack.lm is not installed", call. = FALSE)
  }

  dat <- data.frame(x = mf$x, y = mf$y, z = mf$z)
  form <- spline_model_formula(n_seams, fix)

  if (!is.null(control$maxit) && is.null(control$maxiter)) {
    control$maxiter <- control$maxit
    control$maxit <- NULL
  }
  if (!is.null(control$maxiter) && control$maxiter > 1024) {
    control$maxiter <- 1024
  }
  ctl <- do.call(
    minpack.lm::nls.lm.control,
    utils::modifyList(list(maxiter = 1024), control)
  )
  fit <- minpack.lm::nlsLM(form, data = dat, start = as.list(start), control = ctl)
  par <- c(stats::coef(fit), fix)
  list(
    par = par,
    rss = rss_for(par, mf, n_seams),
    convergence = 0L,
    solver = "minpack.lm::nlsLM",
    raw = fit
  )
}

# Dependency-free fallback solver used when minpack.lm is unavailable or when
# callers explicitly set `prefer_nlslm = FALSE`.
fit_with_optim <- function(start, mf, n_seams, control, fix = NULL) {
  maxit <- control$maxit %||% 10000
  reltol <- control$reltol %||% 1e-12

  objective <- function(p) {
    names(p) <- names(start)
    rss_for(c(p, fix), mf, n_seams)
  }
  nm <- stats::optim(start, objective, method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = reltol)
  )
  bfgs <- stats::optim(nm$par, objective, method = "BFGS",
    control = list(maxit = maxit, reltol = reltol)
  )
  par <- bfgs$par
  names(par) <- names(start)
  par <- c(par, fix)
  list(
    par = par,
    rss = bfgs$value,
    convergence = bfgs$convergence,
    solver = "stats::optim(Nelder-Mead+BFGS)",
    raw = list(nelder_mead = nm, bfgs = bfgs)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Fit an Edwards-Parry congruence spline surface
#'
#' Estimate one- or two-seam spline regression surfaces for congruence
#' hypotheses using nonlinear least squares. The one-seam model matches
#' Equation 39 in Edwards and Parry (2018); the two-seam model extends it with
#' a second estimated seam term.
#'
#' @param data A data frame containing the component and outcome variables.
#' @param wanted,actual Column names, quoted or unquoted, for the two component
#'   variables. `wanted` is modeled as `x`; `actual` is modeled as `y`.
#' @param outcome Column name, quoted or unquoted, for the response variable.
#' @param n_seams Number of seams to estimate. Supported values are `1` and
#'   `2`.
#' @param center,scale Centering and scaling of the component variables, passed
#'   to the internal model-frame builder. Defaults to pooled centering and
#'   pooled scaling (the congruence-safe `RSA` convention). See
#'   [prepare_congruence_data()] for the accepted values and the logical
#'   backward-compatibility mapping.
#' @param hinge_offset Numeric half-width used for the fixed-seam starting
#'   values, in working-scale units. See [prepare_congruence_data()].
#' @param starts Optional named numeric vector or list of named numeric vectors
#'   with starting values. One-seam names are `b0`, `b1`, `b2`, `b3`, `c0`,
#'   and `c1`. Two-seam names are `b0`, `b1`, `b2`, `b3`, `b4`, `c10`, `c11`,
#'   `c20`, and `c21`.
#' @param fix Optional named numeric vector of parameters to hold constant
#'   rather than estimate, for example `c(c0 = 0, c1 = 1)` to fix the one-seam
#'   line of congruence at `Y = X`. Held-constant parameters are dropped from
#'   the optimization, excluded from the residual degrees of freedom, and given
#'   zero covariance, yielding a constrained fit suitable for a nested
#'   comparison against the unconstrained spline (see
#'   [compare_congruence_models()]).
#' @param multistart Logical. If `TRUE`, fit from Stata-style constrained and
#'   unconstrained starts plus small jittered variants and keep the lowest-RSS
#'   converged solution.
#' @param control List of solver control options. `maxit` is accepted for
#'   convenience and translated to `minpack.lm`'s `maxiter`.
#' @param prefer_nlslm Logical. If `TRUE` and `minpack.lm` is installed, use
#'   `minpack.lm::nlsLM()`. Otherwise use a base-R `optim()` fallback.
#'
#' @return A `congruence_spline` object, a list containing coefficients,
#'   covariance matrix, fitted values, residuals, RSS, R-squared, residual
#'   degrees of freedom, solver metadata, and all start attempts. Any
#'   held-constant parameters are recorded in the `fixed` element.
#'
#' @examples
#' \dontrun{
#' fit <- fit_spline_congruence(dat, wanted, actual, satisfaction, n_seams = 1)
#' coef(fit)
#' surface_features(fit)
#' }
#'
#' @export
fit_spline_congruence <- function(data, wanted, actual, outcome,
                                  n_seams = c(1L, 2L),
                                  center = c("pooled", "variablewise", "none"),
                                  scale = c("pooled", "none"),
                                  hinge_offset = 1,
                                  starts = NULL,
                                  fix = NULL,
                                  multistart = TRUE,
                                  control = list(maxit = 10000),
                                  prefer_nlslm = TRUE) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  n_seams <- as.integer(match.arg(as.character(n_seams), c("1", "2")))
  fix <- normalize_fix(fix, n_seams)
  if (is.logical(center) && missing(scale)) scale <- "none"
  mf <- coerce_complete_model_frame(data, wanted, actual, outcome, center, scale)
  start_grid <- make_start_grid(mf, n_seams, starts, multistart, hinge_offset, fix)

  attempts <- vector("list", nrow(start_grid))
  for (i in seq_len(nrow(start_grid))) {
    start <- unlist(start_grid[i, ], use.names = TRUE)
    attempts[[i]] <- tryCatch({
      if (isTRUE(prefer_nlslm) && requireNamespace("minpack.lm", quietly = TRUE)) {
        fit_with_nlslm(start, mf, n_seams, control, fix)
      } else {
        fit_with_optim(start, mf, n_seams, control, fix)
      }
    }, error = function(e) {
      list(error = conditionMessage(e), rss = Inf, start = start)
    })
    attempts[[i]]$start <- start
  }

  rss <- vapply(attempts, function(x) x$rss %||% Inf, numeric(1))
  if (!any(is.finite(rss))) {
    errors <- vapply(attempts, `[[`, "", "error")
    stop(
      "All nonlinear fits failed: ",
      paste(errors, collapse = "; "),
      call. = FALSE
    )
  }

  best <- attempts[[which.min(rss)]]
  par <- best$par[spline_param_names(n_seams)]
  n_free <- length(par) - length(fix)
  fitted <- spline_prediction(par, mf$x, mf$y, n_seams)
  resid <- mf$z - fitted
  tss <- sum((mf$z - mean(mf$z))^2)
  vcov <- vcov_from_jacobian(par, mf, n_seams, fix)

  out <- list(
    call = match.call(),
    n_seams = n_seams,
    coefficients = par,
    fixed = fix,
    vcov = vcov,
    data = mf,
    fitted.values = fitted,
    residuals = resid,
    rss = sum(resid^2),
    tss = tss,
    r.squared = 1 - sum(resid^2) / tss,
    df.residual = nrow(mf) - n_free,
    solver = best$solver,
    convergence = best$convergence,
    attempts = attempts,
    warnings = character()
  )
  class(out) <- "congruence_spline"
  out
}

#' Print a congruence spline model
#'
#' @param x A `congruence_spline` object.
#' @param ... Additional arguments passed to `print()`.
#'
#' @return `x`, invisibly.
#'
#' @export
print.congruence_spline <- function(x, ...) {
  cat("Congruence spline regression\n")
  cat("  seams:", x$n_seams, "\n")
  cat("  solver:", x$solver, "\n")
  cat("  RSS:", format(x$rss, digits = 7), "\n")
  cat("  R-squared:", format(x$r.squared, digits = 5), "\n\n")
  print(coef(x), ...)
  invisible(x)
}

#' Extract congruence spline coefficients
#'
#' @param object A `congruence_spline` object.
#' @param ... Unused.
#'
#' @return A named numeric vector of model coefficients.
#'
#' @export
coef.congruence_spline <- function(object, ...) {
  object$coefficients
}

#' Extract a congruence spline covariance matrix
#'
#' @param object A `congruence_spline` object.
#' @param ... Unused.
#'
#' @return A numeric covariance matrix. Values are `NA` if `numDeriv` is not
#'   installed.
#'
#' @export
vcov.congruence_spline <- function(object, ...) {
  object$vcov
}

# Seam (kink) parameter names for a one- or two-seam model.
seam_param_names <- function(n_seams) {
  if (n_seams == 1L) c("c0", "c1") else c("c10", "c11", "c20", "c21")
}

# Warn that delta-method SEs touching free seam parameters are approximate
# because the numeric Jacobian is taken at the non-differentiable kink.
warn_seam_delta <- function(fit) {
  free_seams <- setdiff(seam_param_names(fit$n_seams), names(fit$fixed))
  if (length(free_seams) > 0L) {
    warning(
      "Delta-method standard errors for seam parameter(s) ",
      paste(free_seams, collapse = ", "),
      " and quantities derived from them are approximate: the numeric Jacobian ",
      "is evaluated at the non-differentiable seam. Use bootstrap_spline() for ",
      "inference on seam features.",
      call. = FALSE
    )
  }
}

#' Log-likelihood of a fitted congruence spline
#'
#' Gaussian log-likelihood implied by the residual sum of squares, enabling
#' [AIC()] and [BIC()] for `congruence_spline` objects.
#'
#' @param object A `congruence_spline` object.
#' @param ... Unused.
#'
#' @return An object of class `logLik` with `df` and `nobs` attributes. The
#'   degrees of freedom count the free (estimated) parameters plus the residual
#'   variance.
#'
#' @export
logLik.congruence_spline <- function(object, ...) {
  n <- nobs(object)
  n_free <- length(object$coefficients) - length(object$fixed)
  ll <- -0.5 * n * (log(2 * pi) + log(object$rss / n) + 1)
  attr(ll, "df") <- n_free + 1L
  attr(ll, "nobs") <- n
  class(ll) <- "logLik"
  ll
}

#' Number of observations in a fitted congruence spline
#'
#' @param object A `congruence_spline` object.
#' @param ... Unused.
#'
#' @return The number of complete cases used to fit the model.
#'
#' @export
nobs.congruence_spline <- function(object, ...) {
  nrow(object$data)
}

#' Summarize a fitted congruence spline
#'
#' Assemble a coefficient table with delta-method standard errors, t
#' statistics, and p-values alongside fit statistics.
#'
#' @param object A `congruence_spline` object.
#' @param warn_seam Logical. If `TRUE` (default), warn that delta-method
#'   standard errors for seam parameters are approximate. See [spline_tests()].
#' @param ... Unused.
#'
#' @return An object of class `summary.congruence_spline`: a list with a
#'   `coefficients` data frame and the fit summary statistics.
#'
#' @export
summary.congruence_spline <- function(object, warn_seam = TRUE, ...) {
  if (isTRUE(warn_seam)) {
    warn_seam_delta(object)
  }
  est <- object$coefficients
  v <- object$vcov
  se <- sqrt(pmax(diag(v), 0))
  fixed_names <- names(object$fixed)
  se[names(est) %in% fixed_names] <- NA_real_
  tval <- est / se
  pval <- 2 * stats::pt(abs(tval), df = object$df.residual, lower.tail = FALSE)

  coefs <- data.frame(
    term = names(est),
    estimate = unname(est),
    std.error = unname(se),
    statistic = unname(tval),
    p.value = unname(pval),
    row.names = NULL
  )

  mf <- object$data
  out <- list(
    coefficients = coefs,
    n_seams = object$n_seams,
    fixed = object$fixed,
    r.squared = object$r.squared,
    rss = object$rss,
    df.residual = object$df.residual,
    nobs = nobs(object),
    solver = object$solver,
    aic = stats::AIC(object),
    center_method = attr(mf, "center_method"),
    scale_method = attr(mf, "scale_method")
  )
  class(out) <- "summary.congruence_spline"
  out
}

#' @export
print.summary.congruence_spline <- function(x, ...) {
  cat("Congruence spline regression\n")
  cat("  seams:", x$n_seams, "\n")
  if (length(x$fixed)) {
    cat("  fixed:", paste(sprintf("%s=%g", names(x$fixed), x$fixed), collapse = ", "), "\n")
  }
  cat("  centering:", x$center_method, " scaling:", x$scale_method, "\n")
  cat("  solver:", x$solver, "\n\n")
  coefs <- x$coefficients
  rownames(coefs) <- coefs$term
  coefs$term <- NULL
  stats::printCoefmat(as.matrix(coefs), has.Pvalue = TRUE, na.print = "")
  cat(sprintf(
    "\nResidual SE on %d df | R-squared: %.4f | RSS: %.4g | AIC: %.2f | n: %d\n",
    x$df.residual, x$r.squared, x$rss, x$aic, x$nobs
  ))
  invisible(x)
}

#' Compute derived surface features for a fitted spline model
#'
#' Derive interpretable features from the fitted surface, including side
#' slopes, along-seam slopes, symmetry expressions, and seam shifts along lines
#' parallel to `Y = -X`.
#'
#' @param fit A `congruence_spline` object.
#' @param lines Numeric vector giving the `k` values where lines of interest
#'   cross `Y = X` at `(k, k)`. For one-seam models, each line is
#'   `Y = 2k - X`.
#'
#' @return A named numeric vector of derived quantities. One-seam models return
#'   the full Edwards-Parry set currently implemented. Two-seam models return
#'   the seam shifts along `Y = -X`, the section slopes in `x` and `y` for each
#'   active-hinge region (base, seam1-only, seam2-only, both), and a crossing
#'   diagnostic: `seams_cross` (1 if the seams cross within the observed `x`
#'   range, else 0), `crossing_x`, and `n_sections` (3 if the seams cross inside
#'   the data, else 4).
#'
#' @examples
#' \dontrun{
#' surface_features(fit, lines = c(-1, 0, 1))
#' }
#'
#' @export
surface_features <- function(fit, lines = c(-1, 0, 1)) {
  stopifnot(inherits(fit, "congruence_spline"))
  p <- coef(fit)

  if (fit$n_seams == 1L) {
    c0 <- p[["c0"]]
    c1 <- p[["c1"]]
    shifts <- stats::setNames(
      vapply(lines, function(k) sqrt(2) * (k - c0 - c1 * k) / (c1 + 1), numeric(1)),
      paste0("shift_y=", 2 * lines, "-x")
    )
    return(c(
      right_intercept = p[["b0"]] - p[["b3"]] * c0,
      right_x_slope = p[["b1"]] - p[["b3"]] * c1,
      right_y_slope = p[["b2"]] + p[["b3"]],
      seam_intercept = p[["b0"]] + p[["b2"]] * c0,
      seam_slope = p[["b1"]] + p[["b2"]] * c1,
      equal_opposite_left = p[["b1"]] + p[["b2"]],
      equal_opposite_right = p[["b1"]] + p[["b2"]] + p[["b3"]] * (1 - c1),
      symmetry_x = 2 * p[["b1"]] - p[["b3"]] * c1,
      symmetry_y = 2 * p[["b2"]] + p[["b3"]],
      shifts
    ))
  }

  b1 <- p[["b1"]]
  b2 <- p[["b2"]]
  b3 <- p[["b3"]]
  b4 <- p[["b4"]]
  c11 <- p[["c11"]]
  c21 <- p[["c21"]]

  crossing_x <- if (abs(c11 - c21) > .Machine$double.eps) {
    (p[["c20"]] - p[["c10"]]) / (c11 - c21)
  } else {
    NA_real_
  }
  x_range <- range(fit$data$x)
  seams_cross <- is.finite(crossing_x) &&
    crossing_x >= x_range[1] && crossing_x <= x_range[2]

  c(
    seam1_shift_y_neg_x = sqrt(2) * p[["c10"]] / (c11 + 1),
    seam2_shift_y_neg_x = sqrt(2) * p[["c20"]] / (c21 + 1),
    base_x_slope = b1,
    base_y_slope = b2,
    seam1_x_slope = b1 - b3 * c11,
    seam1_y_slope = b2 + b3,
    seam2_x_slope = b1 - b4 * c21,
    seam2_y_slope = b2 + b4,
    both_x_slope = b1 - b3 * c11 - b4 * c21,
    both_y_slope = b2 + b3 + b4,
    crossing_x = crossing_x,
    seams_cross = as.numeric(seams_cross),
    n_sections = if (seams_cross) 3 else 4
  )
}

# Build a high-contrast diverging color ramp from a vector of colors or the
# name of an hcl.colors() palette.
spline_palette <- function(palette, n) {
  if (length(palette) == 1L && is.character(palette)) {
    return(grDevices::hcl.colors(n, palette))
  }
  grDevices::colorRampPalette(palette)(n)
}

# Map a fitted surface matrix to a (nrow-1) x (ncol-1) matrix of facet colors,
# shading each facet by the mean predicted outcome of its four corners.
surface_facet_colors <- function(z_mat, palette, n) {
  nx <- nrow(z_mat)
  ny <- ncol(z_mat)
  facet <- (z_mat[-nx, -ny] + z_mat[-1, -ny] + z_mat[-nx, -1] + z_mat[-1, -1]) / 4
  cols <- spline_palette(palette, n)
  if (diff(range(facet)) < .Machine$double.eps) {
    return(matrix(cols[1L], nx - 1L, ny - 1L))
  }
  idx <- as.integer(cut(facet, breaks = n, include.lowest = TRUE))
  matrix(cols[idx], nx - 1L, ny - 1L)
}

#' Plot a fitted congruence spline surface in 3D
#'
#' Draw a three-dimensional response surface for a fitted Edwards-Parry spline
#' model. The surface is evaluated on a regular grid over the observed `x` and
#' `y` ranges and plotted with base R's `persp()`. By default the facets are
#' shaded by predicted outcome on a high-contrast diverging palette (as the
#' `RSA` package does) so the surface shape is legible, and the estimated seam
#' line(s) are projected onto it.
#'
#' @param fit A `congruence_spline` object.
#' @param grid_size Integer number of grid points per axis. The default of `25`
#'   yields larger, higher-contrast facet tiles than a dense grid.
#' @param xlim,ylim Optional axis limits. If either is `NULL` and
#'   `equal_limits = TRUE`, both axes use one shared symmetric range.
#' @param equal_limits Logical. If `TRUE`, use equal limits for the `x` and
#'   `y` axes so the line of congruence is visually anchored by points such as
#'   `(-2, -2)` and `(2, 2)`.
#' @param theta,phi Viewing angles passed to `persp()`.
#' @param expand Expansion factor passed to `persp()`.
#' @param color_by Facet shading. `"outcome"` (default) colors each facet by its
#'   predicted outcome; `"none"` uses the flat `col` fill.
#' @param palette Either a vector of colors to interpolate or the name of a
#'   [grDevices::hcl.colors()] palette used when `color_by = "outcome"`.
#' @param n_color Number of color bins for outcome shading.
#' @param col Flat surface fill color used when `color_by = "none"`.
#' @param border Facet border color. Defaults to `NA` (no border) so colored
#'   tiles read as solid blocks.
#' @param ticktype Tick type passed to `persp()`.
#' @param xlab,ylab,zlab Axis labels.
#' @param main Plot title.
#' @param show_seams Logical. If `TRUE`, draw the estimated seam line(s) (the
#'   fitted ridge of the surface).
#' @param show_fit_line Logical. If `TRUE` (default), draw the purely
#'   theoretical line of congruence, `X = Y`, on the surface. Set to `FALSE` to
#'   drop it for models where congruence along `X = Y` is not the hypothesis.
#' @param show_congruence Deprecated alias for `show_fit_line`, retained for
#'   backward compatibility.
#' @param show_incongruence Logical. If `TRUE`, draw the line of incongruence,
#'   `X = -Y`, on the fitted surface.
#' @param seam_col,seam_lwd Seam line color and width.
#' @param congruence_col,congruence_lty,congruence_lwd Line style for `X = Y`.
#' @param incongruence_col,incongruence_lty,incongruence_lwd Line style for
#'   `X = -Y`.
#' @param ... Additional arguments passed to `persp()`.
#'
#' @return Invisibly returns a list containing the grid vectors, fitted surface
#'   matrix, `persp()` transformation matrix, and plotted seam coordinates.
#'
#' @examples
#' \dontrun{
#' plot_spline_surface(fit)
#' plot_spline_surface(fit, show_fit_line = FALSE)
#' }
#'
#' @export
plot_spline_surface <- function(fit,
                                grid_size = 25,
                                xlim = NULL,
                                ylim = NULL,
                                equal_limits = TRUE,
                                theta = -35,
                                phi = 25,
                                expand = 0.65,
                                color_by = c("outcome", "none"),
                                palette = c("#a50026", "#d73027", "#f46d43",
                                            "#fdae61", "#fee08b", "#ffffbf",
                                            "#d9ef8b", "#a6d96a", "#66bd63",
                                            "#1a9850", "#006837"),
                                n_color = 16,
                                col = "lightblue",
                                border = NA,
                                ticktype = "detailed",
                                xlab = "Wanted (centered)",
                                ylab = "Actual (centered)",
                                zlab = "Outcome",
                                main = "Congruence spline surface",
                                show_seams = TRUE,
                                show_fit_line = TRUE,
                                show_congruence = NULL,
                                show_incongruence = FALSE,
                                seam_col = "#2166ac",
                                seam_lwd = 2.5,
                                congruence_col = "grey25",
                                congruence_lty = 2,
                                congruence_lwd = 1.2,
                                incongruence_col = "grey45",
                                incongruence_lty = 3,
                                incongruence_lwd = 1.1,
                                ...) {
  stopifnot(inherits(fit, "congruence_spline"))
  color_by <- match.arg(color_by)
  if (!is.null(show_congruence)) {
    show_fit_line <- show_congruence
  }

  grid_size <- as.integer(grid_size)
  if (!is.finite(grid_size) || grid_size < 5L) {
    stop("grid_size must be an integer of at least 5.", call. = FALSE)
  }

  if (is.null(xlim)) {
    xlim <- range(fit$data$x)
  }
  if (is.null(ylim)) {
    ylim <- range(fit$data$y)
  }
  if (isTRUE(equal_limits)) {
    lim <- range(xlim, ylim)
    lim <- max(abs(lim))
    xlim <- c(-lim, lim)
    ylim <- c(-lim, lim)
  }

  x_seq <- seq(xlim[1], xlim[2], length.out = grid_size)
  y_seq <- seq(ylim[1], ylim[2], length.out = grid_size)
  grid <- expand.grid(x = x_seq, y = y_seq)
  z_hat <- spline_prediction(coef(fit), grid$x, grid$y, fit$n_seams)
  z_mat <- matrix(z_hat, nrow = length(x_seq), ncol = length(y_seq))

  facet_col <- col
  if (color_by == "outcome") {
    facet_col <- surface_facet_colors(z_mat, palette, n_color)
  }

  trans <- graphics::persp(
    x_seq,
    y_seq,
    z_mat,
    theta = theta,
    phi = phi,
    expand = expand,
    col = facet_col,
    border = border,
    ticktype = ticktype,
    xlab = xlab,
    ylab = ylab,
    zlab = zlab,
    main = main,
    ...
  )

  draw_surface_line <- function(xs, ys, col, lwd, lty) {
    keep <- xs >= min(x_seq) & xs <= max(x_seq) &
      ys >= min(y_seq) & ys <= max(y_seq)
    if (sum(keep) < 2L) {
      return(data.frame(x = numeric(), y = numeric(), z = numeric()))
    }
    xs <- xs[keep]
    ys <- ys[keep]
    zs <- spline_prediction(coef(fit), xs, ys, fit$n_seams)
    projected <- grDevices::trans3d(xs, ys, zs, pmat = trans)
    graphics::lines(projected, col = col, lwd = lwd, lty = lty)
    data.frame(x = xs, y = ys, z = zs)
  }

  congruence <- NULL
  incongruence <- NULL
  if (isTRUE(show_fit_line)) {
    congruence <- draw_surface_line(
      x_seq,
      x_seq,
      col = congruence_col,
      lwd = congruence_lwd,
      lty = congruence_lty
    )
  }
  if (isTRUE(show_incongruence)) {
    incongruence <- draw_surface_line(
      x_seq,
      -x_seq,
      col = incongruence_col,
      lwd = incongruence_lwd,
      lty = incongruence_lty
    )
  }

  seams <- list()
  if (isTRUE(show_seams)) {
    seam_params <- if (fit$n_seams == 1L) {
      list(c(c0 = coef(fit)[["c0"]], c1 = coef(fit)[["c1"]]))
    } else {
      list(
        c(c0 = coef(fit)[["c10"]], c1 = coef(fit)[["c11"]]),
        c(c0 = coef(fit)[["c20"]], c1 = coef(fit)[["c21"]])
      )
    }

    for (i in seq_along(seam_params)) {
      c0 <- seam_params[[i]][["c0"]]
      c1 <- seam_params[[i]][["c1"]]
      xs <- x_seq
      ys <- c0 + c1 * xs
      seams[[i]] <- draw_surface_line(
        xs,
        ys,
        col = seam_col,
        lwd = seam_lwd,
        lty = 1
      )
    }
  }

  invisible(list(
    x = x_seq,
    y = y_seq,
    z = z_mat,
    transform = trans,
    congruence = congruence,
    incongruence = incongruence,
    seams = seams
  ))
}

#' Plot a fitted congruence spline surface as a 2D contour map
#'
#' Draw the fitted spline surface as a filled `ggplot2` heatmap with contour
#' lines, the estimated seam line(s), and (optionally) the line of congruence.
#' A 2D view is often clearer than [plot_spline_surface()] for reading where the
#' seams fall. Requires the `ggplot2` package.
#'
#' @param fit A `congruence_spline` object.
#' @param grid_size Integer number of grid points per axis.
#' @param xlim,ylim Optional axis limits. If either is `NULL` and
#'   `equal_limits = TRUE`, both axes use one shared symmetric range.
#' @param equal_limits Logical. If `TRUE`, use one shared symmetric range for
#'   both axes.
#' @param palette Either a vector of colors to interpolate or the name of a
#'   [grDevices::hcl.colors()] palette for the fill gradient.
#' @param bins Number of contour bins.
#' @param show_seams Logical. If `TRUE`, overlay the estimated seam line(s).
#' @param show_fit_line Logical. If `TRUE` (default), overlay the theoretical
#'   line of congruence `X = Y`.
#' @param seam_col,seam_lwd Seam line color and width.
#' @param fit_line_col,fit_line_lty Line style for `X = Y`.
#' @param xlab,ylab,fill_lab,main Axis, legend, and title labels.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' plot_spline_contour(fit)
#' }
#'
#' @export
plot_spline_contour <- function(fit,
                                grid_size = 80,
                                xlim = NULL,
                                ylim = NULL,
                                equal_limits = TRUE,
                                palette = c("#a50026", "#d73027", "#f46d43",
                                            "#fdae61", "#fee08b", "#ffffbf",
                                            "#d9ef8b", "#a6d96a", "#66bd63",
                                            "#1a9850", "#006837"),
                                bins = 12,
                                show_seams = TRUE,
                                show_fit_line = TRUE,
                                seam_col = "#2166ac",
                                seam_lwd = 1.1,
                                fit_line_col = "grey25",
                                fit_line_lty = "dashed",
                                xlab = "Wanted (centered)",
                                ylab = "Actual (centered)",
                                fill_lab = "Outcome",
                                main = "Congruence spline surface") {
  stopifnot(inherits(fit, "congruence_spline"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_spline_contour().", call. = FALSE)
  }

  if (is.null(xlim)) xlim <- range(fit$data$x)
  if (is.null(ylim)) ylim <- range(fit$data$y)
  if (isTRUE(equal_limits)) {
    lim <- max(abs(range(xlim, ylim)))
    xlim <- c(-lim, lim)
    ylim <- c(-lim, lim)
  }

  x_seq <- seq(xlim[1], xlim[2], length.out = grid_size)
  y_seq <- seq(ylim[1], ylim[2], length.out = grid_size)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$z <- spline_prediction(coef(fit), grid$x, grid$y, fit$n_seams)
  cols <- spline_palette(palette, 256)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, z = z)) +
    ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = TRUE) +
    ggplot2::geom_contour(colour = "grey20", bins = bins, linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(colours = cols, name = fill_lab) +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = xlab, y = ylab, title = main)

  if (isTRUE(show_fit_line)) {
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0,
      colour = fit_line_col, linetype = fit_line_lty
    )
  }
  if (isTRUE(show_seams)) {
    seam_params <- if (fit$n_seams == 1L) {
      list(c(c0 = coef(fit)[["c0"]], c1 = coef(fit)[["c1"]]))
    } else {
      list(
        c(c0 = coef(fit)[["c10"]], c1 = coef(fit)[["c11"]]),
        c(c0 = coef(fit)[["c20"]], c1 = coef(fit)[["c21"]])
      )
    }
    for (sp in seam_params) {
      p <- p + ggplot2::geom_abline(
        slope = sp[["c1"]], intercept = sp[["c0"]],
        colour = seam_col, linewidth = seam_lwd
      )
    }
  }
  p
}

# Compute a scalar delta-method test from a function of model parameters.
delta_stat <- function(fit, expr_fun, null = 0, name = "test") {
  p <- coef(fit)
  v <- vcov(fit)
  estimate <- expr_fun(p)
  if (anyNA(v)) {
    return(data.frame(
      term = name,
      estimate = estimate,
      se = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    ))
  }
  grad <- numDeriv::grad(function(theta) {
    names(theta) <- names(p)
    expr_fun(theta)
  }, unname(p))
  se <- sqrt(drop(t(grad) %*% v %*% grad))
  statistic <- (estimate - null) / se
  data.frame(
    term = name,
    estimate = estimate,
    se = se,
    statistic = statistic,
    p.value = 2 * stats::pt(abs(statistic), df = fit$df.residual, lower.tail = FALSE),
    row.names = NULL
  )
}

# Compute a joint Wald F test for parameter functions.
wald_joint <- function(fit, funcs, nulls, name) {
  p <- coef(fit)
  v <- vcov(fit)
  est <- vapply(funcs, function(f) f(p), numeric(1))
  if (anyNA(v)) {
    return(data.frame(term = name, df = length(est), statistic = NA_real_, p.value = NA_real_))
  }
  jac <- numDeriv::jacobian(function(theta) {
    names(theta) <- names(p)
    vapply(funcs, function(f) f(theta), numeric(1))
  }, unname(p))
  middle <- jac %*% v %*% t(jac)
  diff <- est - nulls
  stat <- drop(t(diff) %*% pinv(middle) %*% diff) / length(diff)
  data.frame(
    term = name,
    df = length(diff),
    statistic = stat,
    p.value = stats::pf(stat, df1 = length(diff), df2 = fit$df.residual, lower.tail = FALSE),
    row.names = NULL
  )
}

#' Run Edwards-Parry-style Wald tests for a spline model
#'
#' Compute normal-theory delta and Wald tests for quantities used to interpret
#' congruence spline surfaces. One-seam models receive the most complete set of
#' tests, following the Stata code distributed with Edwards and Parry (2018).
#' Two-seam models receive the omnibus and seam-contribution tests present in
#' the Stata scripts.
#'
#' @param fit A `congruence_spline` object.
#' @param warn_seam Logical. If `TRUE` (default), warn that the delta-method
#'   standard errors underlying these tests are approximate for quantities
#'   involving the seam parameters, because the numeric Jacobian is evaluated at
#'   the non-differentiable seam. Edwards and Parry (2018) recommend
#'   [bootstrap_spline()] for inference on seam features.
#'
#' @return A list with two data frames: `scalar` for one-degree-of-freedom
#'   delta-method tests and `joint` for multi-constraint Wald F tests.
#'
#' @examples
#' \dontrun{
#' tests <- spline_tests(fit)
#' tests$joint
#' }
#'
#' @export
spline_tests <- function(fit, warn_seam = TRUE) {
  stopifnot(inherits(fit, "congruence_spline"))
  if (isTRUE(warn_seam)) {
    warn_seam_delta(fit)
  }

  if (fit$n_seams == 1L) {
    scalar <- rbind(
      delta_stat(fit, function(p) p[["b0"]] - p[["b3"]] * p[["c0"]], 0, "right_intercept"),
      delta_stat(fit, function(p) p[["b1"]] - p[["b3"]] * p[["c1"]], 0, "right_x_slope"),
      delta_stat(fit, function(p) p[["b2"]] + p[["b3"]], 0, "right_y_slope"),
      delta_stat(fit, function(p) p[["b0"]] + p[["b2"]] * p[["c0"]], 0, "seam_intercept"),
      delta_stat(fit, function(p) p[["b1"]] + p[["b2"]] * p[["c1"]], 0, "seam_slope"),
      delta_stat(fit, function(p) sqrt(2) * p[["c0"]] / (p[["c1"]] + 1), 0, "shift_y=-x"),
      delta_stat(fit, function(p) {
        -sqrt(2) * (1 - p[["c0"]] - p[["c1"]]) / (p[["c1"]] + 1)
      }, 0, "shift_y=2-x"),
      delta_stat(fit, function(p) {
        sqrt(2) * (1 + p[["c0"]] - p[["c1"]]) / (p[["c1"]] + 1)
      }, 0, "shift_y=-2-x"),
      delta_stat(fit, function(p) p[["b1"]] + p[["b2"]], 0, "equal_opposite_left"),
      delta_stat(fit, function(p) {
        p[["b1"]] + p[["b2"]] + p[["b3"]] * (1 - p[["c1"]])
      }, 0, "equal_opposite_right"),
      delta_stat(fit, function(p) 2 * p[["b1"]] - p[["b3"]] * p[["c1"]], 0, "symmetry_x"),
      delta_stat(fit, function(p) 2 * p[["b2"]] + p[["b3"]], 0, "symmetry_y")
    )
    joint <- rbind(
      wald_joint(fit, list(
        function(p) p[["b3"]], function(p) p[["c0"]], function(p) p[["c1"]]
      ), c(0, 0, 0), "deviation_from_no_seam"),
      wald_joint(fit, list(
        function(p) p[["b1"]] + p[["b2"]],
        function(p) 2 * p[["b1"]] - p[["b3"]],
        function(p) p[["c0"]],
        function(p) p[["c1"]]
      ), c(0, 0, 0, 1), "absolute_difference_constraints"),
      wald_joint(fit, list(
        function(p) p[["c0"]], function(p) p[["c1"]]
      ), c(0, 1), "seam_equals_y_equals_x")
    )
    return(list(scalar = scalar, joint = joint))
  }

  list(
    scalar = data.frame(),
    joint = rbind(
      wald_joint(fit, lapply(names(coef(fit))[-1], function(nm) {
        force(nm)
        function(p) p[[nm]]
      }), rep(0, length(coef(fit)) - 1), "omnibus_non_intercept"),
      wald_joint(fit, list(
        function(p) p[["b3"]], function(p) p[["b4"]],
        function(p) p[["c10"]], function(p) p[["c11"]],
        function(p) p[["c20"]], function(p) p[["c21"]]
      ), rep(0, 6), "all_seam_terms"),
      wald_joint(fit, list(
        function(p) p[["b3"]], function(p) p[["c10"]], function(p) p[["c11"]]
      ), rep(0, 3), "seam1_terms"),
      wald_joint(fit, list(
        function(p) p[["b4"]], function(p) p[["c20"]], function(p) p[["c21"]]
      ), rep(0, 3), "seam2_terms")
    )
  )
}

#' Fit OLS comparison models for congruence analyses
#'
#' Fit the absolute-difference, linear, one-break piecewise, constrained
#' one-seam piecewise, and optionally fixed two-seam OLS models used as
#' comparisons and starting-value sources for spline regression.
#'
#' @param data A data frame containing the component and outcome variables.
#' @param wanted,actual Column names, quoted or unquoted, for the two component
#'   variables. `wanted` is modeled as `x`; `actual` is modeled as `y`.
#' @param outcome Column name, quoted or unquoted, for the response variable.
#' @param center,scale Centering and scaling of the component variables, passed
#'   to [prepare_congruence_data()]. Defaults to pooled centering and pooled
#'   scaling.
#' @param hinge_offset Numeric half-width, in working-scale units, for the
#'   fixed two-seam OLS seam lines `Y = X ± hinge_offset`. See
#'   [prepare_congruence_data()].
#' @param n_seams Number of seams to include in comparison models. If `2`, add
#'   the fixed two-seam OLS model using seam lines `Y = X + hinge_offset` and
#'   `Y = X - hinge_offset`.
#' @param constrained Logical. If `FALSE`, omit the constrained one-seam OLS
#'   model from the returned list.
#'
#' @return A named list of `lm` objects with class `congruence_piecewise`.
#'
#' @examples
#' \dontrun{
#' fits <- fit_piecewise_congruence(dat, wanted, actual, satisfaction)
#' tidy_lm_summary(fits)
#' }
#'
#' @export
fit_piecewise_congruence <- function(data, wanted, actual, outcome,
                                     center = c("pooled", "variablewise", "none"),
                                     scale = c("pooled", "none"),
                                     hinge_offset = 1,
                                     n_seams = c(1L, 2L),
                                     constrained = TRUE) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  n_seams <- as.integer(match.arg(as.character(n_seams), c("1", "2")))
  if (is.logical(center) && missing(scale)) scale <- "none"
  mf <- prepare_congruence_data(data, wanted, actual, outcome, center, scale, hinge_offset)

  fits <- list(
    absolute_difference = stats::lm(z ~ abs_diff, data = mf),
    linear = stats::lm(z ~ x + y, data = mf),
    one_break = stats::lm(z ~ x + y + w + xw + yw, data = mf),
    constrained_one_seam = stats::lm(z ~ x + y + zd, data = mf)
  )
  if (n_seams == 2L) {
    fits$fixed_two_seam <- stats::lm(z ~ x + y + hinge_upper + hinge_lower, data = mf)
  }
  if (!isTRUE(constrained)) {
    fits$constrained_one_seam <- NULL
  }

  class(fits) <- "congruence_piecewise"
  fits
}

#' Summarize OLS comparison models
#'
#' Convert the list returned by `fit_piecewise_congruence()` to a compact data
#' frame with estimates and model R-squared values.
#'
#' @param fits A `congruence_piecewise` object or named list of `lm` objects.
#'
#' @return A data frame with columns `model`, `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, and `r.squared`.
#'
#' @export
tidy_lm_summary <- function(fits) {
  do.call(rbind, lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    smry <- summary(fit)
    coef_tab <- stats::coef(smry)
    data.frame(
      model = nm,
      term = rownames(coef_tab),
      estimate = unname(coef_tab[, "Estimate"]),
      std.error = unname(coef_tab[, "Std. Error"]),
      statistic = unname(coef_tab[, 3L]),
      p.value = unname(coef_tab[, 4L]),
      r.squared = unname(smry$r.squared),
      row.names = NULL
    )
  }))
}

# Residual sum of squares for an `lm` or `congruence_spline` model.
model_rss <- function(m) {
  if (inherits(m, "congruence_spline")) m$rss else sum(stats::residuals(m)^2)
}

# Model R-squared for an `lm` or `congruence_spline` model.
model_r2 <- function(m) {
  if (inherits(m, "congruence_spline")) m$r.squared else summary(m)$r.squared
}

#' Compare nested congruence models
#'
#' Build a sequential nested-model comparison table for any mix of `lm`
#' comparison models (from [fit_piecewise_congruence()]) and spline fits (from
#' [fit_spline_congruence()]). Each adjacent pair contributes a change in
#' R-squared and an F test on the change in residual sum of squares, alongside
#' per-model AIC.
#'
#' Pass the models in increasing complexity, for example the Edwards-Parry chain
#' absolute-difference \eqn{\subset} linear \eqn{\subset} one-break piecewise
#' \eqn{\subset} spline, or a one-seam spline followed by a two-seam spline. If a
#' larger model has a higher RSS than the simpler model nested within it
#' (a sign of a nonlinear local minimum), its F and p-value are returned as `NA`
#' and a warning recommends keeping the simpler model.
#'
#' @param ... Two or more fitted models in nested order. Names, if supplied, are
#'   used as row labels.
#'
#' @return A data frame with one row per model: `model`, `npar`, `df.residual`,
#'   `rss`, `r.squared`, `AIC`, and the comparison columns `df`, `deltaR2`, `F`,
#'   and `p.value` (the first row's comparison columns are `NA`).
#'
#' @examples
#' \dontrun{
#' ols <- fit_piecewise_congruence(dat, wanted, actual, satisfaction)
#' spline <- fit_spline_congruence(dat, wanted, actual, satisfaction)
#' compare_congruence_models(
#'   absdiff = ols$absolute_difference,
#'   linear = ols$linear,
#'   piecewise = ols$one_break,
#'   spline = spline
#' )
#' }
#'
#' @export
compare_congruence_models <- function(...) {
  models <- list(...)
  if (length(models) < 2L) {
    stop("Provide at least two nested models to compare.", call. = FALSE)
  }
  ok <- vapply(models, function(m) inherits(m, c("lm", "congruence_spline")), logical(1))
  if (!all(ok)) {
    stop("All models must be `lm` or `congruence_spline` objects.", call. = FALSE)
  }

  labels <- names(models)
  if (is.null(labels)) labels <- rep("", length(models))
  blank <- !nzchar(labels)
  labels[blank] <- paste0("model", seq_along(models))[blank]

  n <- vapply(models, function(m) as.integer(stats::nobs(m)), integer(1))
  if (length(unique(n)) > 1L) {
    warning("Models were fit on different numbers of observations; comparisons may be invalid.", call. = FALSE)
  }
  dfres <- vapply(models, function(m) as.numeric(stats::df.residual(m)), numeric(1))
  rss <- vapply(models, model_rss, numeric(1))
  r2 <- vapply(models, model_r2, numeric(1))
  aic <- vapply(models, function(m) as.numeric(stats::AIC(m)), numeric(1))
  npar <- n - dfres

  K <- length(models)
  df_diff <- F_stat <- p_val <- dR2 <- rep(NA_real_, K)
  r2_drop <- FALSE
  for (i in 2:K) {
    rss0 <- rss[i - 1L]
    rss1 <- rss[i]
    df0 <- dfres[i - 1L]
    df1 <- dfres[i]
    dR2[i] <- r2[i] - r2[i - 1L]
    df_diff[i] <- df0 - df1
    if (df_diff[i] <= 0 || df1 <= 0) {
      next
    }
    if (rss1 > rss0 + 1e-9) {
      r2_drop <- TRUE
      next
    }
    F_stat[i] <- ((rss0 - rss1) / df_diff[i]) / (rss1 / df1)
    p_val[i] <- stats::pf(F_stat[i], df_diff[i], df1, lower.tail = FALSE)
  }
  if (r2_drop) {
    warning(
      "A larger model fit worse (higher RSS) than the simpler model nested ",
      "within it, likely a nonlinear local minimum; its F/p are NA. Prefer the ",
      "simpler model.",
      call. = FALSE
    )
  }

  data.frame(
    model = labels,
    npar = npar,
    df.residual = dfres,
    rss = rss,
    r.squared = r2,
    AIC = aic,
    df = df_diff,
    deltaR2 = dR2,
    F = F_stat,
    p.value = p_val,
    row.names = NULL
  )
}

#' @rdname compare_congruence_models
#' @param object A `congruence_spline` object (the first model).
#' @export
anova.congruence_spline <- function(object, ...) {
  models <- c(list(object), list(...))
  do.call(compare_congruence_models, models)
}

#' Bootstrap coefficients and surface features for a spline model
#'
#' Refit the spline model to bootstrap resamples and compute confidence
#' intervals for coefficients and derived surface features via
#' [boot::boot.ci()]. Edwards and Parry (2018) recommend the bootstrap over
#' delta-method inference because the seam parameters sit at a non-differentiable
#' kink (see [spline_tests()]).
#'
#' @param fit A `congruence_spline` object.
#' @param R Number of bootstrap resamples. Edwards and Parry used 10,000 in the
#'   example analyses; smaller values are useful for smoke tests.
#' @param conf Confidence level for intervals.
#' @param type Interval type passed to [boot::boot.ci()]: `"perc"` (percentile,
#'   the default), `"basic"`, `"norm"` (normal approximation), or `"bca"`
#'   (bias-corrected and accelerated). If `boot.ci()` cannot compute the
#'   requested interval for a degenerate statistic, the percentile interval is
#'   used as a fallback.
#'
#' @return A data frame with columns `term`, `estimate`, `lower`, and `upper`.
#'   The proportion of bootstrap resamples whose refit failed (and were dropped)
#'   is attached as the `fail_rate` attribute; a non-zero failure rate also
#'   triggers a warning.
#'
#' @examples
#' \dontrun{
#' bootstrap_spline(fit, R = 100, type = "bca")
#' }
#'
#' @export
bootstrap_spline <- function(fit, R = 10000, conf = 0.95,
                             type = c("perc", "basic", "norm", "bca")) {
  stopifnot(inherits(fit, "congruence_spline"))
  type <- match.arg(type)
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for bootstrap_spline().", call. = FALSE)
  }
  comp_name <- c(perc = "percent", basic = "basic", norm = "normal", bca = "bca")[[type]]

  mf <- fit$data
  stat <- function(data, indices) {
    boot_mf <- data[indices, , drop = FALSE]
    starts <- as.list(coef(fit))
    refit <- tryCatch(
      fit_spline_congruence(
        boot_mf, "x", "y", "z",
        n_seams = fit$n_seams,
        center = FALSE,
        scale = "none",
        starts = starts,
        fix = fit$fixed,
        multistart = FALSE,
        prefer_nlslm = requireNamespace("minpack.lm", quietly = TRUE)
      ),
      error = function(e) NULL
    )
    if (is.null(refit)) {
      return(rep(NA_real_, length(coef(fit)) + length(surface_features(fit))))
    }
    c(coef(refit), surface_features(refit))
  }

  boot_obj <- boot::boot(mf, stat, R = R)
  original <- c(coef(fit), surface_features(fit))
  alpha <- (1 - conf) / 2

  fail_rate <- mean(!stats::complete.cases(boot_obj$t))
  if (fail_rate > 0) {
    warning(
      sprintf(
        "%.1f%% of bootstrap refits failed and were dropped from the intervals.",
        100 * fail_rate
      ),
      call. = FALSE
    )
  }

  ci_for <- function(i) {
    vals <- boot_obj$t[, i]
    ok <- is.finite(vals)
    percentile <- function() {
      if (sum(ok) < 2L) {
        return(c(NA_real_, NA_real_))
      }
      unname(stats::quantile(vals[ok], probs = c(alpha, 1 - alpha), na.rm = TRUE))
    }
    if (sum(ok) < 2L || stats::sd(vals[ok]) == 0) {
      return(percentile())
    }
    ci <- tryCatch(
      suppressWarnings(boot::boot.ci(boot_obj, conf = conf, type = type, index = i)),
      error = function(e) NULL
    )
    comp <- ci[[comp_name]]
    if (is.null(comp)) {
      return(percentile())
    }
    nc <- ncol(comp)
    c(comp[1L, nc - 1L], comp[1L, nc])
  }

  cis <- t(vapply(seq_along(original), ci_for, numeric(2)))
  out <- data.frame(
    term = names(original),
    estimate = unname(original),
    lower = cis[, 1L],
    upper = cis[, 2L],
    row.names = NULL
  )
  attr(out, "fail_rate") <- fail_rate
  out
}
