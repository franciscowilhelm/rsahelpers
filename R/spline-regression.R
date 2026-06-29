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

# Build the complete-case model frame used by all public fitting functions.
coerce_complete_model_frame <- function(data, wanted, actual, outcome, center = TRUE) {
  wanted <- as_variable_name(wanted, substitute(wanted), "wanted")
  actual <- as_variable_name(actual, substitute(actual), "actual")
  outcome <- as_variable_name(outcome, substitute(outcome), "outcome")

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
  if (isTRUE(center)) {
    x <- x - mean(x)
    y <- y - mean(y)
  }

  data.frame(
    x = x,
    y = y,
    z = z[keep],
    .row = which(keep),
    check.names = FALSE
  )
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
#' @param center Logical. If `TRUE`, subtract each component's mean before
#'   constructing model terms. Use `FALSE` when passing variables that are
#'   already centered, such as the downloaded `*HC` and `*WC` columns.
#'
#' @return A data frame with standardized columns `x`, `y`, `z`, `.row`, and
#'   helper columns for absolute difference, one-break piecewise regression,
#'   constrained piecewise regression, and fixed two-seam comparisons.
#'
#' @examples
#' \dontrun{
#' prepared <- prepare_congruence_data(my_data, wanted, actual, satisfaction)
#' }
#'
#' @export
prepare_congruence_data <- function(data, wanted, actual, outcome, center = TRUE) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  mf <- coerce_complete_model_frame(data, wanted, actual, outcome, center)
  diff <- mf$y - mf$x
  w <- as.numeric(diff < 0)

  mf$abs_diff <- abs(diff)
  mf$w <- w
  mf$xw <- mf$x * w
  mf$yw <- mf$y * w
  mf$zd <- diff * w
  mf$hinge_upper <- (mf$y - 1 - mf$x) * as.numeric(mf$y < 1 + mf$x)
  mf$hinge_lower <- (mf$y + 1 - mf$x) * as.numeric(mf$y < -1 + mf$x)
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
vcov_from_jacobian <- function(par, mf, n_seams) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    return(matrix(NA_real_, length(par), length(par), dimnames = list(names(par), names(par))))
  }

  pred_fun <- function(p) {
    names(p) <- names(par)
    spline_prediction(p, mf$x, mf$y, n_seams)
  }
  jac <- numDeriv::jacobian(pred_fun, unname(par))
  rss <- rss_for(par, mf, n_seams)
  df_resid <- nrow(mf) - length(par)
  xtx_inv <- pinv(crossprod(jac))
  out <- (rss / df_resid) * xtx_inv
  dimnames(out) <- list(names(par), names(par))
  out
}

# Derive Stata-style starting values from constrained or unconstrained
# piecewise OLS models.
start_from_piecewise <- function(mf, n_seams = 1L, method = c("constrained", "unconstrained")) {
  method <- match.arg(method)

  if (n_seams == 1L && method == "constrained") {
    fit <- stats::lm(z ~ x + y + zd, data = prepare_hinges_for_lm(mf))
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
    fit <- stats::lm(z ~ x + y + w + xw + yw, data = prepare_hinges_for_lm(mf))
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
    fit <- stats::lm(z ~ x + y + h1_fixed + h2_fixed, data = prepare_hinges_for_lm(mf))
    cf <- stats::coef(fit)
    return(c(
      b0 = unname(cf[["(Intercept)"]]),
      b1 = unname(cf[["x"]]),
      b2 = unname(cf[["y"]]),
      b3 = unname(cf[["h1_fixed"]]),
      b4 = unname(cf[["h2_fixed"]]),
      c10 = 1,
      c11 = 1,
      c20 = -1,
      c21 = 1
    ))
  }

  one <- start_from_piecewise(mf, 1L, "unconstrained")
  c(
    b0 = unname(one[["b0"]]),
    b1 = unname(one[["b1"]]),
    b2 = unname(one[["b2"]]),
    b3 = unname(one[["b3"]]) / 2,
    b4 = unname(one[["b3"]]) / 2,
    c10 = unname(one[["c0"]]),
    c11 = unname(one[["c1"]]),
    c20 = -1,
    c21 = 1
  )
}

# Add fixed seam and one-break terms used to compute starting values.
prepare_hinges_for_lm <- function(mf) {
  out <- mf
  diff <- out$y - out$x
  out$w <- as.numeric(diff < 0)
  out$xw <- out$x * out$w
  out$yw <- out$y * out$w
  out$zd <- diff * out$w
  out$h1_fixed <- spline_hinge(out$x, out$y, 1, 1)
  out$h2_fixed <- spline_hinge(out$x, out$y, -1, 1)
  out
}

# Build the candidate start matrix. Multiple starts are important because the
# two-seam surface can have local minima and label-switched seam solutions.
make_start_grid <- function(mf, n_seams, starts = NULL, multistart = TRUE) {
  base <- if (n_seams == 1L) {
    rbind(
      constrained = start_from_piecewise(mf, 1L, "constrained"),
      unconstrained = start_from_piecewise(mf, 1L, "unconstrained")
    )
  } else {
    rbind(
      constrained = start_from_piecewise(mf, 2L, "constrained"),
      unconstrained = start_from_piecewise(mf, 2L, "unconstrained")
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

  unique(as.data.frame(base))
}

# Primary nonlinear solver using Levenberg-Marquardt via minpack.lm.
fit_with_nlslm <- function(start, mf, n_seams, control) {
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("minpack.lm is not installed", call. = FALSE)
  }

  dat <- data.frame(x = mf$x, y = mf$y, z = mf$z)
  if (n_seams == 1L) {
    form <- z ~ b0 + b1 * x + b2 * y + b3 * (y - c0 - c1 * x) * (y < c0 + c1 * x)
  } else {
    form <- z ~ b0 + b1 * x + b2 * y +
      b3 * (y - c10 - c11 * x) * (y < c10 + c11 * x) +
      b4 * (y - c20 - c21 * x) * (y < c20 + c21 * x)
  }

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
  par <- stats::coef(fit)
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
fit_with_optim <- function(start, mf, n_seams, control) {
  maxit <- control$maxit %||% 10000
  reltol <- control$reltol %||% 1e-12

  objective <- function(p) {
    names(p) <- names(start)
    rss_for(p, mf, n_seams)
  }
  nm <- stats::optim(start, objective, method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = reltol)
  )
  bfgs <- stats::optim(nm$par, objective, method = "BFGS",
    control = list(maxit = maxit, reltol = reltol)
  )
  par <- bfgs$par
  names(par) <- names(start)
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
#' @param center Logical. If `TRUE`, subtract each component's mean before
#'   fitting. Use `FALSE` for pre-centered columns.
#' @param starts Optional named numeric vector or list of named numeric vectors
#'   with starting values. One-seam names are `b0`, `b1`, `b2`, `b3`, `c0`,
#'   and `c1`. Two-seam names are `b0`, `b1`, `b2`, `b3`, `b4`, `c10`, `c11`,
#'   `c20`, and `c21`.
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
#'   degrees of freedom, solver metadata, and all start attempts.
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
                                  center = TRUE,
                                  starts = NULL,
                                  multistart = TRUE,
                                  control = list(maxit = 10000),
                                  prefer_nlslm = TRUE) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  n_seams <- as.integer(match.arg(as.character(n_seams), c("1", "2")))
  mf <- coerce_complete_model_frame(data, wanted, actual, outcome, center)
  start_grid <- make_start_grid(mf, n_seams, starts, multistart)

  attempts <- vector("list", nrow(start_grid))
  for (i in seq_len(nrow(start_grid))) {
    start <- unlist(start_grid[i, ], use.names = TRUE)
    attempts[[i]] <- tryCatch({
      if (isTRUE(prefer_nlslm) && requireNamespace("minpack.lm", quietly = TRUE)) {
        fit_with_nlslm(start, mf, n_seams, control)
      } else {
        fit_with_optim(start, mf, n_seams, control)
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
  par <- best$par
  fitted <- spline_prediction(par, mf$x, mf$y, n_seams)
  resid <- mf$z - fitted
  tss <- sum((mf$z - mean(mf$z))^2)
  vcov <- vcov_from_jacobian(par, mf, n_seams)

  out <- list(
    call = match.call(),
    n_seams = n_seams,
    coefficients = par,
    vcov = vcov,
    data = mf,
    fitted.values = fitted,
    residuals = resid,
    rss = sum(resid^2),
    tss = tss,
    r.squared = 1 - sum(resid^2) / tss,
    df.residual = nrow(mf) - length(par),
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
#'   seam shifts along `Y = -X` for each seam.
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

  c(
    seam1_shift_y_neg_x = sqrt(2) * p[["c10"]] / (p[["c11"]] + 1),
    seam2_shift_y_neg_x = sqrt(2) * p[["c20"]] / (p[["c21"]] + 1)
  )
}

#' Plot a fitted congruence spline surface in 3D
#'
#' Draw a three-dimensional response surface for a fitted Edwards-Parry spline
#' model. The surface is evaluated on a regular grid over the observed `x` and
#' `y` ranges and plotted with base R's `persp()`. Estimated seam line(s) are
#' projected onto the fitted surface and overlaid in red.
#'
#' @param fit A `congruence_spline` object.
#' @param grid_size Integer number of grid points per axis.
#' @param xlim,ylim Optional axis limits. If either is `NULL` and
#'   `equal_limits = TRUE`, both axes use one shared symmetric range.
#' @param equal_limits Logical. If `TRUE`, use equal limits for the `x` and
#'   `y` axes so the line of congruence is visually anchored by points such as
#'   `(-2, -2)` and `(2, 2)`.
#' @param theta,phi Viewing angles passed to `persp()`.
#' @param expand Expansion factor passed to `persp()`.
#' @param col Surface fill color.
#' @param border Surface border color.
#' @param ticktype Tick type passed to `persp()`.
#' @param xlab,ylab,zlab Axis labels.
#' @param main Plot title.
#' @param show_seams Logical. If `TRUE`, draw estimated seam line(s).
#' @param show_congruence Logical. If `TRUE`, draw the line of congruence,
#'   `X = Y`, on the fitted surface.
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
#' }
#'
#' @export
plot_spline_surface <- function(fit,
                                grid_size = 60,
                                xlim = NULL,
                                ylim = NULL,
                                equal_limits = TRUE,
                                theta = -35,
                                phi = 25,
                                expand = 0.65,
                                col = "lightblue",
                                border = "grey70",
                                ticktype = "detailed",
                                xlab = "Wanted (centered)",
                                ylab = "Actual (centered)",
                                zlab = "Outcome",
                                main = "Congruence spline surface",
                                show_seams = TRUE,
                                show_congruence = TRUE,
                                show_incongruence = FALSE,
                                seam_col = "red3",
                                seam_lwd = 2,
                                congruence_col = "black",
                                congruence_lty = 1,
                                congruence_lwd = 2,
                                incongruence_col = "grey30",
                                incongruence_lty = 2,
                                incongruence_lwd = 1.5,
                                ...) {
  stopifnot(inherits(fit, "congruence_spline"))

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

  trans <- graphics::persp(
    x_seq,
    y_seq,
    z_mat,
    theta = theta,
    phi = phi,
    expand = expand,
    col = col,
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
  if (isTRUE(show_congruence)) {
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
spline_tests <- function(fit) {
  stopifnot(inherits(fit, "congruence_spline"))

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
#' @param center Logical. If `TRUE`, subtract each component's mean before
#'   constructing model terms.
#' @param n_seams Number of seams to include in comparison models. If `2`, add
#'   the fixed two-seam OLS model using seam lines `Y = X + 1` and `Y = X - 1`.
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
                                     center = TRUE,
                                     n_seams = c(1L, 2L),
                                     constrained = TRUE) {
  wanted <- resolve_public_variable(wanted, substitute(wanted), "wanted")
  actual <- resolve_public_variable(actual, substitute(actual), "actual")
  outcome <- resolve_public_variable(outcome, substitute(outcome), "outcome")
  n_seams <- as.integer(match.arg(as.character(n_seams), c("1", "2")))
  mf <- prepare_congruence_data(data, wanted, actual, outcome, center)

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
#' @return A data frame with columns `model`, `term`, `estimate`, and
#'   `r.squared`.
#'
#' @export
tidy_lm_summary <- function(fits) {
  do.call(rbind, lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    data.frame(
      model = nm,
      term = names(stats::coef(fit)),
      estimate = unname(stats::coef(fit)),
      r.squared = unname(summary(fit)$r.squared),
      row.names = NULL
    )
  }))
}

#' Bootstrap coefficients and surface features for a spline model
#'
#' Refit the spline model to bootstrap resamples and compute percentile
#' confidence intervals for coefficients and derived surface features.
#'
#' @param fit A `congruence_spline` object.
#' @param R Number of bootstrap resamples. Edwards and Parry used 10,000 in the
#'   example analyses; smaller values are useful for smoke tests.
#' @param conf Confidence level for intervals.
#' @param type Interval type. Currently accepted for API compatibility; this
#'   implementation returns percentile intervals for all supported values.
#'
#' @return A data frame with columns `term`, `estimate`, `lower`, and `upper`.
#'
#' @examples
#' \dontrun{
#' bootstrap_spline(fit, R = 100)
#' }
#'
#' @export
bootstrap_spline <- function(fit, R = 10000, conf = 0.95, type = c("perc", "basic", "norm")) {
  stopifnot(inherits(fit, "congruence_spline"))
  type <- match.arg(type)
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for bootstrap_spline().", call. = FALSE)
  }

  mf <- fit$data
  stat <- function(data, indices) {
    boot_mf <- data[indices, , drop = FALSE]
    starts <- as.list(coef(fit))
    refit <- tryCatch(
      fit_spline_congruence(
        boot_mf, "x", "y", "z",
        n_seams = fit$n_seams,
        center = FALSE,
        starts = starts,
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
  cis <- t(vapply(seq_along(original), function(i) {
    vals <- boot_obj$t[, i]
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      return(c(lower = NA_real_, upper = NA_real_))
    }
    stats::quantile(vals, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  }, numeric(2)))
  colnames(cis) <- c("lower", "upper")
  data.frame(
    term = names(original),
    estimate = unname(original),
    lower = cis[, "lower"],
    upper = cis[, "upper"],
    row.names = NULL
  )
}
