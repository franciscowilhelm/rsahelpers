# classify-unconstrained.R -- Yao & Ma (2023) congruence-type classification
# from the *unconstrained* second-order polynomial.
#
# Exported from the rsasim project (lib/estimate.R). Two tiers, both reading
# the a-/p-parameters off a single full quadratic fit:
#   * classify_yao_ma()     -- tolerance-based, from point estimates (fast,
#                              for the fast lm() path and Monte Carlo loops);
#   * classify_yao_ma_rsa() -- significance-based, from a fitted RSA::RSA()
#                              object (robust lavaan CIs; optional bootstrap
#                              CIs for the ratio-based p10/p11).
#
# Requires surface-math.R (rsa_a_parameters, rsa_p_parameters,
# lateral_shift_quantity).

# ---- Fast lm() path ---------------------------------------------------------

#' Fit the full second-order polynomial by OLS
#'
#' @param data A data frame with centered predictors `P`, `E` and response
#'   `outcome`.
#' @return An [stats::lm] object.
#' @export
fit_quadratic_lm <- function(data) {
  stats::lm(outcome ~ P + E + I(P^2) + I(P * E) + I(E^2), data = data)
}

#' Extract polynomial b-coefficients from a fitted quadratic lm
#'
#' @param model An [stats::lm] fit from [fit_quadratic_lm].
#' @return A one-row b-tibble (`b0`--`b5`).
#' @export
extract_b_coefficients <- function(model) {
  cf <- stats::coef(model)
  tibble::tibble(
    b0 = unname(cf[["(Intercept)"]]),
    b1 = unname(cf[["P"]]),
    b2 = unname(cf[["E"]]),
    b3 = unname(cf[["I(P^2)"]]),
    b4 = unname(cf[["I(P * E)"]]),
    b5 = unname(cf[["I(E^2)"]])
  )
}

#' Assemble b-coefficients, a-parameters, p-parameters, and the lateral shift
#'
#' @param b A b-tibble (`b0`--`b5`).
#' @return A one-row tibble binding `b`, [rsa_a_parameters], [rsa_p_parameters],
#'   and `lsq` = [lateral_shift_quantity].
#' @export
estimates_from_b <- function(b) {
  dplyr::bind_cols(
    b,
    rsa_a_parameters(b),
    rsa_p_parameters(b),
    tibble::tibble(lsq = lateral_shift_quantity(b))
  )
}

#' Estimates row directly from a fitted quadratic lm
#'
#' @param model An [stats::lm] fit from [fit_quadratic_lm].
#' @return See [estimates_from_b].
#' @export
estimates_from_lm <- function(model) {
  estimates_from_b(extract_b_coefficients(model))
}

#' Negate all polynomial coefficients (valence flip)
#'
#' Turns a negative-valence (valley) surface into the positive-valence (ridge)
#' surface the typology conditions are written for. Negation swaps the
#' principal axes (`p11(-b) = p21(b)`, `p10(-b) = p20(b)`), so this is not a
#' mere sign reversal.
#'
#' @param b A b-tibble (`b0`--`b5`).
#' @return A b-tibble with every coefficient negated.
#' @export
flip_valence_b <- function(b) {
  tibble::as_tibble(lapply(b, function(x) -x))
}

# ---- Yao & Ma typology from decision flags ----------------------------------

#' Map Yao & Ma decision flags to a Table-4 type
#'
#' Table 4 core: `type = base(misfit-line, contingency) + offset(fit-line)`.
#' `base` from SDE (`p10 = 0`) vs CDE (`p10 != 0`) and contingency
#' (`p11 != 1`); `offset` from the fit line -- Flat / LLE (`a1 != 0`) /
#' SLE (`a2 != 0`) / CLE (both).
#'
#' @param a4_negative Logical: inverted-U along the LOIC (`a4 < 0`). If not
#'   `TRUE`, the surface is outside the typology and a descriptive label is
#'   returned.
#' @param a1_nonzero,a2_nonzero Logical fit-line flags.
#' @param p10_nonzero,p11_not_one Logical misfit-line flags.
#' @param a3_nonzero Logical: directional difference effect (used only to label
#'   out-of-typology surfaces).
#' @return A one-row [tibble::tibble] with `type` (integer, `NA` outside the
#'   typology), `label`, and `in_typology` (logical).
#' @export
yao_ma_type_from_flags <- function(a4_negative, a1_nonzero, a2_nonzero,
                                   p10_nonzero, p11_not_one,
                                   a3_nonzero = FALSE) {
  if (is.na(a4_negative) || !a4_negative) {
    label <- if (isTRUE(a3_nonzero) && !a1_nonzero) {
      "outside typology: directional difference effect (a3 without a4)"
    } else if (a1_nonzero && !isTRUE(a3_nonzero)) {
      "outside typology: level effect without congruence (a1 without a4)"
    } else if (a1_nonzero || isTRUE(a3_nonzero)) {
      "outside typology: level and difference effects without congruence"
    } else {
      "no congruence effect"
    }
    return(tibble::tibble(type = NA_integer_, label = label, in_typology = FALSE))
  }

  if (is.na(p10_nonzero) || is.na(p11_not_one) || is.na(a1_nonzero) || is.na(a2_nonzero)) {
    return(tibble::tibble(
      type = NA_integer_,
      label = "unclassifiable: degenerate p-parameters",
      in_typology = FALSE
    ))
  }

  base <- if (!p11_not_one && !p10_nonzero) 1L
  else if (!p11_not_one && p10_nonzero) 5L
  else if (p11_not_one && !p10_nonzero) 9L
  else 13L

  offset <- if (!a1_nonzero && !a2_nonzero) 0L
  else if (a1_nonzero && !a2_nonzero) 1L
  else if (!a1_nonzero && a2_nonzero) 2L
  else 3L

  type <- base + offset
  fit_line <- c("", " & LLE", " & SLE", " & CLE")[offset + 1]
  archetype <- if (base %in% c(1L, 9L)) "Exact correspondence" else "Commensurate compatibility"
  contingency <- if (base >= 9L) " with contingency" else ""
  tibble::tibble(
    type = type,
    label = paste0(archetype, fit_line, contingency),
    in_typology = TRUE
  )
}

# ---- Tier 1: tolerance-based classification of point estimates ---------------

#' Classify a point-estimate surface into a Yao & Ma type (tolerance tier)
#'
#' For negatively valenced outcomes the congruence signature is `a4 > 0`; the
#' typology applies to the sign-flipped surface, so the negated coefficients
#' are classified. Isotropic surfaces (undefined FPA direction) adopt the
#' LOC-parallel axis: `p11 := 1`, `p10 := Y0 - X0`.
#'
#' Tolerance defaults suit the centered 1-5 scale: true effects are
#' `|a| >= 0.4` and `|p-shifts| >= 0.24`, while sampling SDs at `n >= 500` are
#' well below half these cutoffs, so `eps` sits between noise and signal.
#'
#' @param est Output of [estimates_from_b] / [estimates_from_lm]
#'   (b-coefficients must be present).
#' @param valence `"positive"` or `"negative"`.
#' @param eps_a,eps_p10,eps_p11 Tolerances for the a-parameters, `p10`, and
#'   `p11 - 1`.
#' @return The [yao_ma_type_from_flags] tibble with a `tier` column.
#' @export
classify_yao_ma <- function(est, valence = "positive",
                            eps_a = 0.10, eps_p10 = 0.15, eps_p11 = 0.15) {
  if (valence == "negative") {
    est <- estimates_from_b(flip_valence_b(est[c("b0", "b1", "b2", "b3", "b4", "b5")]))
  }
  p10 <- est[["p10"]]
  p11 <- est[["p11"]]
  if (!is.null(est[["degenerate"]]) && isTRUE(est[["degenerate"]] == "isotropic")) {
    p11 <- 1
    p10 <- est[["Y0"]] - est[["X0"]]
  }

  a4_negative <- est[["a4"]] < -eps_a
  res <- yao_ma_type_from_flags(
    a4_negative = a4_negative,
    a1_nonzero = abs(est[["a1"]]) > eps_a,
    a2_nonzero = abs(est[["a2"]]) > eps_a,
    p10_nonzero = if (a4_negative) abs(p10) > eps_p10 else FALSE,
    p11_not_one = if (a4_negative) abs(p11 - 1) > eps_p11 else FALSE,
    a3_nonzero = abs(est[["a3"]]) > eps_a
  )
  res$tier <- "tolerance"
  res
}

# ---- RSA package path (significance tier) ------------------------------------

#' Fit RSA::RSA() with project defaults
#'
#' Data must contain centered `P`, `E` and `outcome`. Predictors are assumed
#' pre-centered (`center = "none"`); simulated data is clean (`out.rm = FALSE`,
#' also much faster).
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param models Passed to [RSA::RSA] (default `"full"`).
#' @param ... Further arguments to [RSA::RSA] (e.g. `add =`).
#' @return An `RSA` object.
#' @export
fit_rsa <- function(data, models = "full", ...) {
  RSA::RSA(outcome ~ P * E, data = as.data.frame(data),
    models = models, center = "none", scale = "none",
    out.rm = FALSE, verbose = FALSE, ...
  )
}

#' Tidy one model's coefficient table from a fitted RSA object
#'
#' Labels are normalized to `PA1_curv`-style naming.
#'
#' @param fit An `RSA` object from [fit_rsa].
#' @param model Model name within `fit` (default `"full"`).
#' @return A [tibble::tibble] with `param, est, se, pvalue, ci_lower, ci_upper`.
#' @export
tidy_rsa_params <- function(fit, model = "full") {
  cf <- RSA::getPar(fit, "coef", model = model)
  tibble::tibble(
    param = sub("\\.", "_", cf$label),
    est = cf$est,
    se = cf$se,
    pvalue = cf$pvalue,
    ci_lower = cf$ci.lower,
    ci_upper = cf$ci.upper
  )
}

#' Classify a fitted RSA surface into a Yao & Ma type (significance tier)
#'
#' "Nonzero" = the CI excludes 0; contingency = the `p11` CI excludes 1;
#' congruence = the `a4` CI lies entirely below 0. With `boot_ci = TRUE`,
#' percentile-bootstrap CIs replace the delta-method CIs for the ratio-based
#' `p10`/`p11` (the Yao & Ma recommendation for showcase fits).
#'
#' @param fit An `RSA` object (must include the `"full"` model, or the model
#'   named in `model`).
#' @param model Model name to read parameters from (default `"full"`).
#' @param boot_ci Logical; use bootstrap CIs for `p10`/`p11`.
#' @param R Bootstrap replicates when `boot_ci = TRUE`.
#' @return The [yao_ma_type_from_flags] tibble with a `tier` column.
#' @export
classify_yao_ma_rsa <- function(fit, model = "full", boot_ci = FALSE, R = 1000) {
  params <- tidy_rsa_params(fit, model = model)
  get <- function(name) params[params$param == name, ]

  if (boot_ci) {
    bci <- as.data.frame(stats::confint(fit, model = model, method = "boot", R = R))
    bci$param <- sub("\\.", "_", rownames(bci))
    for (nm in c("p10", "p11")) {
      row <- bci[bci$param == nm, ]
      if (nrow(row) == 1) {
        params[params$param == nm, c("ci_lower", "ci_upper")] <- row[1, 1:2]
      }
    }
    get <- function(name) params[params$param == name, ]
  }

  excludes <- function(name, value) {
    row <- get(name)
    nrow(row) == 1 && !is.na(row$ci_lower) &&
      (row$ci_lower > value || row$ci_upper < value)
  }

  res <- yao_ma_type_from_flags(
    a4_negative = {
      row <- get("a4")
      nrow(row) == 1 && !is.na(row$ci_upper) && row$ci_upper < 0
    },
    a1_nonzero = excludes("a1", 0),
    a2_nonzero = excludes("a2", 0),
    p10_nonzero = excludes("p10", 0),
    p11_not_one = excludes("p11", 1),
    a3_nonzero = excludes("a3", 0)
  )
  res$tier <- if (boot_ci) "significance (bootstrap p10/p11)" else "significance"
  res
}
