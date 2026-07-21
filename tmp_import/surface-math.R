# surface-math.R -- pure algebra for quadratic PE fit response surfaces.
#
# Exported from the rsasim project (lib/surface_math.R). Only the pieces the
# classification and coverage exports depend on are carried over; the surface
# predictors and the local-quadratic-approximation helpers were left behind.
#
# Conventions:
#   Quadratic model:  O = b0 + b1*P + b2*E + b3*P^2 + b4*P*E + b5*E^2
#   L/D basis:        L = P + E,  D = P - E
#                     O = b0 + r*L + s*L^2 + m*D - q*D^2 + c*L*D
#   D > 0 = deficiency (P > E), D < 0 = excess (E > P).
# A "b-tibble" throughout is a one-row tibble/data.frame with columns
# b0, b1, b2, b3, b4, b5.

# ---- L/D basis <-> b-coefficients -------------------------------------------

#' Convert L/D-basis coefficients to polynomial b-coefficients
#'
#' @param b0,r,s,m,q,c L/D-basis coefficients (see file header).
#' @return A one-row [tibble::tibble] with columns `b0`--`b5`.
#' @export
ld_to_b <- function(b0 = 0, r = 0, s = 0, m = 0, q = 0, c = 0) {
  tibble::tibble(
    b0 = b0,
    b1 = r + m,
    b2 = r - m,
    b3 = s - q + c,
    b4 = 2 * s + 2 * q,
    b5 = s - q - c
  )
}

#' Convert polynomial b-coefficients to the L/D basis
#'
#' @param b A b-tibble (columns `b0`--`b5`).
#' @return A one-row [tibble::tibble] with columns `b0, r, s, m, q, c`.
#' @export
b_to_ld <- function(b) {
  tibble::tibble(
    b0 = b[["b0"]],
    r = (b[["b1"]] + b[["b2"]]) / 2,
    s = (b[["b3"]] + b[["b4"]] + b[["b5"]]) / 4,
    m = (b[["b1"]] - b[["b2"]]) / 2,
    q = (b[["b4"]] - b[["b3"]] - b[["b5"]]) / 4,
    c = (b[["b3"]] - b[["b5"]]) / 2
  )
}

# ---- RSA parameters ---------------------------------------------------------

#' Response-surface a-parameters
#'
#' @param b A b-tibble (columns `b0`--`b5`).
#' @return A one-row [tibble::tibble] with `a1`--`a5`.
#' @export
rsa_a_parameters <- function(b) {
  tibble::tibble(
    a1 = b[["b1"]] + b[["b2"]],
    a2 = b[["b3"]] + b[["b4"]] + b[["b5"]],
    a3 = b[["b1"]] - b[["b2"]],
    a4 = b[["b3"]] - b[["b4"]] + b[["b5"]],
    a5 = b[["b3"]] - b[["b5"]]
  )
}

#' Lateral shift quantity
#'
#' Location `t*` of the extremum along the line of incongruence (LOIC),
#' parameterized as `P = t, E = -t`. In the L/D basis, `t* = m / (4q)`.
#'
#' @param b A b-tibble (columns `b0`--`b5`).
#' @param tol Curvature tolerance below which the shift is undefined (`NA`).
#' @return A length-1 numeric, or `NA_real_` when `a4` is ~0.
#' @export
lateral_shift_quantity <- function(b, tol = 1e-10) {
  a <- rsa_a_parameters(b)
  if (abs(a[["a4"]]) < tol) {
    return(NA_real_)
  }
  -a[["a3"]] / (2 * a[["a4"]])
}

#' Stationary point, principal axes, and eigenvalues of a quadratic surface
#'
#' Formulas follow the RSA package (`getPar`) / Edwards (2002). Degenerate
#' cases the standard formulas cannot handle are flagged in the `degenerate`
#' column: `"plane"` (no second-order terms), `"isotropic"` (circular
#' paraboloid; stationary point defined, axes not), and `"parabolic"`
#' (rank-1 curvature: a ridge/valley line instead of a stationary point --
#' the case for every pure squared-difference surface).
#'
#' @param b A b-tibble (columns `b0`--`b5`).
#' @param tol Numerical tolerance for the degeneracy tests.
#' @return A one-row [tibble::tibble] with `X0, Y0, p10, p11, p20, p21,
#'   PA1_curv, PA2_curv, l1, l2, degenerate`.
#' @export
rsa_p_parameters <- function(b, tol = 1e-10) {
  b1 <- b[["b1"]]; b2 <- b[["b2"]]
  b3 <- b[["b3"]]; b4 <- b[["b4"]]; b5 <- b[["b5"]]

  out <- tibble::tibble(
    X0 = NA_real_, Y0 = NA_real_,
    p10 = NA_real_, p11 = NA_real_,
    p20 = NA_real_, p21 = NA_real_,
    PA1_curv = NA_real_, PA2_curv = NA_real_,
    l1 = NA_real_, l2 = NA_real_,
    degenerate = NA_character_
  )

  if (max(abs(c(b3, b4, b5))) < tol) {
    out$degenerate <- "plane"
    return(out)
  }

  disc <- sqrt((b3 - b5)^2 + b4^2)
  out$l1 <- ((b3 + b5) + disc) / 2
  out$l2 <- ((b3 + b5) - disc) / 2

  det <- 4 * b3 * b5 - b4^2

  if (disc < tol) {
    # b3 == b5, b4 == 0: axes directionless, but the stationary point exists.
    out$degenerate <- "isotropic"
    out$X0 <- (b2 * b4 - 2 * b1 * b5) / det
    out$Y0 <- (b1 * b4 - 2 * b2 * b3) / det
    return(out)
  }

  if (abs(b4) >= tol) {
    out$p11 <- (b5 - b3 + disc) / b4
    out$p21 <- (b5 - b3 - disc) / b4
  } else {
    # Axis-aligned surface: FPA follows the direction of the larger eigenvalue.
    # Inf encodes a vertical axis (line X = const).
    out$p11 <- if (b3 > b5) 0 else Inf
    out$p21 <- if (b3 > b5) Inf else 0
  }
  out$PA1_curv <- if (is.finite(out$p11)) b3 + b4 * out$p11 + b5 * out$p11^2 else b5
  out$PA2_curv <- if (is.finite(out$p21)) b3 + b4 * out$p21 + b5 * out$p21^2 else b5

  if (abs(det) >= tol) {
    out$X0 <- (b2 * b4 - 2 * b1 * b5) / det
    out$Y0 <- (b1 * b4 - 2 * b2 * b3) / det
    if (is.finite(out$p11)) out$p10 <- out$Y0 - out$p11 * out$X0
    if (is.finite(out$p21)) out$p20 <- out$Y0 - out$p21 * out$X0
    return(out)
  }

  # Parabolic surface: ridge/valley line instead of a stationary point.
  out$degenerate <- "parabolic"
  l_nz <- if (abs(out$l1) >= abs(out$l2)) out$l1 else out$l2
  ev <- c(b4 / 2, l_nz - b3)
  if (max(abs(ev)) < tol) ev <- c(l_nz - b5, b4 / 2)
  u <- ev[1]; w <- ev[2]
  denom <- u * b4 + 2 * w * b5
  if (abs(denom) >= tol && is.finite(out$p11)) {
    out$p10 <- -(u * b1 + w * b2) / denom
  }
  out
}
