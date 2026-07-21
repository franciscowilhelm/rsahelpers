# classify-constrained.R -- Yao & Ma (2023) classification via the a priori
# *constrained* polynomial models (Table-4 equality conditions as estimable
# constraint sets), plus the constrained-model comparison machinery.
#
# Exported from the rsasim project (lib/model_comparison.R). classify_yao_ma()
# (classify-unconstrained.R) reads the type off one full fit; this route
# instead fits the family of constrained models -- each encoding one type's
# equality conditions -- and selects among them by AICc / constraint
# tenability. classify_yao_ma_constrained() is a NEW convenience wrapper (not
# in the original lib/) that returns the retained model's declared type.
#
# Two fitting paths:
#   * ls     -- constrained least squares (fast; no extra dependencies);
#   * lavaan -- RSA::RSA() named models / custom `add =` constraints (robust
#               MLR fits with scaled chi-square difference tests).
#
# Requires classify-unconstrained.R (fit_quadratic_lm, extract_b_coefficients,
# fit_rsa).

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Constraint lexicon ------------------------------------------------------

# Each linear constraint is w'(b1..b5) == 0.
.b_constraint_lexicon <- list(
  b1_zero = list(w = c(1, 0, 0, 0, 0), lavaan = "b1 == 0"),
  b2_zero = list(w = c(0, 1, 0, 0, 0), lavaan = "b2 == 0"),
  b3_zero = list(w = c(0, 0, 1, 0, 0), lavaan = "b3 == 0"),
  b4_zero = list(w = c(0, 0, 0, 1, 0), lavaan = "b4 == 0"),
  b5_zero = list(w = c(0, 0, 0, 0, 1), lavaan = "b5 == 0"),
  a1_zero = list(w = c(1, 1, 0, 0, 0), lavaan = "b1 + b2 == 0"),
  a3_zero = list(w = c(1, -1, 0, 0, 0), lavaan = "b1 - b2 == 0"),
  a2_zero = list(w = c(0, 0, 1, 1, 1), lavaan = "b3 + b4 + b5 == 0"),
  a5_zero = list(w = c(0, 0, 1, 0, -1), lavaan = "b3 - b5 == 0")
)

# FPA-through-origin (p10 == 0 with p11 free): the stationary point lies on the
# first-principal-axis line through the origin. Cleared of the determinant and
# b4 so no division appears.
.fpa_origin_lavaan <- paste0(
  "b4*(b1*b4 - 2*b2*b3) - ",
  "(b5 - b3 + sqrt((b3 - b5)^2 + b4^2))*(b2*b4 - 2*b1*b5) == 0"
)

#' Residual of the FPA-through-origin constraint
#'
#' @param b A b-tibble (`b0`--`b5`).
#' @return A length-1 numeric; 0 when the FPA passes through the origin.
#' @export
fpa_origin_residual <- function(b) {
  disc <- sqrt((b[["b3"]] - b[["b5"]])^2 + b[["b4"]]^2)
  b[["b4"]] * (b[["b1"]] * b[["b4"]] - 2 * b[["b2"]] * b[["b3"]]) -
    (b[["b5"]] - b[["b3"]] + disc) *
      (b[["b2"]] * b[["b4"]] - 2 * b[["b1"]] * b[["b5"]])
}

# ---- Model specs -------------------------------------------------------------

#' Construct one constrained-model spec
#'
#' @param name,label Identifier and human-readable label.
#' @param linear Character vector of constraint names from the internal lexicon
#'   (`b1_zero`, `a1_zero`, `a2_zero`, `a3_zero`, `a5_zero`, ...).
#' @param fpa_origin Logical; add the nonlinear FPA-through-origin constraint.
#' @param yao_ma_type Integer Yao & Ma type this spec encodes (`NA` for rivals).
#' @param rsa_named Name of the exactly-coinciding RSA-package model, or `NA`.
#' @return A named list describing the spec.
#' @export
surface_model_spec <- function(name, label, linear = character(0),
                               fpa_origin = FALSE,
                               yao_ma_type = NA_integer_,
                               rsa_named = NA_character_) {
  stopifnot(all(linear %in% names(.b_constraint_lexicon)))
  list(
    name = name, label = label, linear = linear, fpa_origin = fpa_origin,
    yao_ma_type = as.integer(yao_ma_type), rsa_named = rsa_named
  )
}

#' Registry of a priori constrained models
#'
#' The 16 Yao & Ma types (type 16 = `full`) plus first-order rivals
#' (`null`, `additive`, `diff`, `mean`, `IA`). `rsa_named` marks models that
#' coincide exactly with an RSA-package named model; the named rotated models
#' (SRSQD/SRRR) are parabolic and strictly narrower than the Table-4
#' contingency conditions, so types 9-15 have no exact named counterpart.
#'
#' @return A named list of specs, keyed by spec name.
#' @export
surface_model_specs <- function() {
  specs <- list(
    surface_model_spec("type_01", "Type 1: squared difference",
      linear = c("b1_zero", "b2_zero", "a5_zero", "a2_zero"),
      yao_ma_type = 1, rsa_named = "SQD"
    ),
    surface_model_spec("type_02", "Type 2: rising ridge",
      linear = c("a3_zero", "a5_zero", "a2_zero"),
      yao_ma_type = 2, rsa_named = "RR"
    ),
    surface_model_spec("type_03", "Type 3: squared difference + SLE",
      linear = c("b1_zero", "b2_zero", "a5_zero"), yao_ma_type = 3
    ),
    surface_model_spec("type_04", "Type 4: rising ridge + CLE",
      linear = c("a3_zero", "a5_zero"), yao_ma_type = 4
    ),
    surface_model_spec("type_05", "Type 5: shifted squared difference",
      linear = c("a1_zero", "a5_zero", "a2_zero"),
      yao_ma_type = 5, rsa_named = "SSQD"
    ),
    surface_model_spec("type_06", "Type 6: shifted rising ridge",
      linear = c("a5_zero", "a2_zero"), yao_ma_type = 6, rsa_named = "SRR"
    ),
    surface_model_spec("type_07", "Type 7: shifted squared difference + SLE",
      linear = c("a1_zero", "a5_zero"), yao_ma_type = 7
    ),
    surface_model_spec("type_08", "Type 8: shifted ridge + CLE",
      linear = "a5_zero", yao_ma_type = 8
    ),
    surface_model_spec("type_09", "Type 9: rotated, ridge through origin",
      linear = c("b1_zero", "b2_zero", "a2_zero"), yao_ma_type = 9
    ),
    surface_model_spec("type_10", "Type 10: rotated rising, FPA through origin",
      linear = "a2_zero", fpa_origin = TRUE, yao_ma_type = 10
    ),
    surface_model_spec("type_11", "Type 11: rotated + SLE, origin optimum",
      linear = c("b1_zero", "b2_zero"), yao_ma_type = 11
    ),
    surface_model_spec("type_12", "Type 12: rotated + CLE, FPA through origin",
      fpa_origin = TRUE, yao_ma_type = 12
    ),
    surface_model_spec("type_13", "Type 13: rotated + shifted, flat fit line",
      linear = c("a1_zero", "a2_zero"), yao_ma_type = 13
    ),
    surface_model_spec("type_14", "Type 14: rotated + shifted + LLE",
      linear = "a2_zero", yao_ma_type = 14
    ),
    surface_model_spec("type_15", "Type 15: rotated + shifted + SLE",
      linear = "a1_zero", yao_ma_type = 15
    ),
    surface_model_spec("full", "Full second-order polynomial (type 16)",
      yao_ma_type = 16, rsa_named = "full"
    ),
    surface_model_spec("null", "Intercept only",
      linear = c("b1_zero", "b2_zero", "b3_zero", "b4_zero", "b5_zero"),
      rsa_named = "null"
    ),
    surface_model_spec("additive", "Additive main effects",
      linear = c("b3_zero", "b4_zero", "b5_zero"), rsa_named = "additive"
    ),
    surface_model_spec("diff", "Difference plane",
      linear = c("a1_zero", "b3_zero", "b4_zero", "b5_zero"),
      rsa_named = "diff"
    ),
    surface_model_spec("mean", "Mean plane",
      linear = c("a3_zero", "b3_zero", "b4_zero", "b5_zero"),
      rsa_named = "mean"
    ),
    surface_model_spec("IA", "Main effects + interaction",
      linear = c("b3_zero", "b5_zero"), rsa_named = "IA"
    )
  )
  stats::setNames(specs, vapply(specs, `[[`, "", "name"))
}

#' Linear-constraint matrix of a spec (rows = constraints on `b1`--`b5`)
#' @param spec A spec from [surface_model_specs].
#' @return A numeric matrix with 5 columns.
#' @export
spec_constraint_matrix <- function(spec) {
  if (length(spec$linear) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 5))
  }
  do.call(rbind, lapply(.b_constraint_lexicon[spec$linear], `[[`, "w"))
}

#' lavaan `add =` syntax for a spec (NA for the unconstrained model)
#' @param spec A spec from [surface_model_specs].
#' @return A length-1 character, or `NA_character_`.
#' @export
spec_lavaan_add <- function(spec) {
  lines <- vapply(.b_constraint_lexicon[spec$linear], `[[`, "", "lavaan")
  if (spec$fpa_origin) lines <- c(lines, .fpa_origin_lavaan)
  if (length(lines) == 0) {
    return(NA_character_)
  }
  paste(lines, collapse = "\n")
}

#' Number of free mean-structure parameters implied by a spec
#' @param spec A spec from [surface_model_specs].
#' @return An integer (`b0` plus free `b` directions).
#' @export
spec_n_free <- function(spec) {
  r <- if (length(spec$linear) == 0) 0L else qr(spec_constraint_matrix(spec))$rank
  5L - r - as.integer(spec$fpa_origin) + 1L
}

# Does a spec imply the FPA passes through the origin (directly, via a pinned
# stationary point, or because the FPA is the LOC itself)?
.spec_implies_fpa_origin <- function(spec) {
  if (spec$fpa_origin) {
    return(TRUE)
  }
  C <- spec_constraint_matrix(spec)
  if (nrow(C) == 0) {
    return(FALSE)
  }
  spans <- function(rows) qr(C)$rank == qr(rbind(C, rows))$rank
  spans(rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0))) ||
    spans(rbind(c(1, -1, 0, 0, 0), c(0, 0, 1, 0, -1)))
}

#' Is model `a` nested in model `b`?
#'
#' `a` is the more constrained model: every constraint of `b` is implied by
#' `a`'s constraints.
#'
#' @param a,b Specs from [surface_model_specs].
#' @return Logical.
#' @export
spec_is_nested <- function(a, b) {
  Ca <- spec_constraint_matrix(a)
  Cb <- spec_constraint_matrix(b)
  lin_ok <- nrow(Cb) == 0 ||
    (nrow(Ca) > 0 && qr(Ca)$rank == qr(rbind(Ca, Cb))$rank)
  fpa_ok <- !b$fpa_origin || .spec_implies_fpa_origin(a)
  lin_ok && fpa_ok
}

# ---- Fast constrained least squares -----------------------------------------

.quad_design <- function(data) {
  cbind(P = data$P, E = data$E, P2 = data$P^2, PE = data$P * data$E, E2 = data$E^2)
}

# Orthonormal basis of the null space {v: Cv = 0}; p columns in C.
.null_space <- function(C, p, tol = 1e-10) {
  if (nrow(C) == 0) {
    return(diag(p))
  }
  s <- svd(C, nu = 0, nv = p)
  r <- sum(s$d > tol * max(s$d))
  if (r >= p) {
    return(matrix(numeric(0), nrow = p, ncol = 0))
  }
  s$v[, seq(r + 1, p), drop = FALSE]
}

#' Constrained least-squares fit of one spec
#'
#' Linear constraint sets reduce to `lm` on the null-space-transformed design.
#' The FPA-through-origin constraint is handled by profiling: for a fixed
#' quadratic part the constraint makes `b1` a multiple of `b2`, so the linear
#' part solves in closed form and only the free quadratic directions are
#' optimized numerically.
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param spec A spec from [surface_model_specs].
#' @return `list(name, b, ssr, n, k_mean, converged)`; `k_mean` counts free
#'   mean-structure parameters.
#' @export
fit_constrained_ls <- function(data, spec) {
  y <- data$outcome
  X <- .quad_design(data)
  n <- length(y)

  if (!spec$fpa_origin) {
    N <- .null_space(spec_constraint_matrix(spec), p = 5)
    Xc <- if (ncol(N) > 0) cbind(1, X %*% N) else matrix(1, n, 1)
    fit <- stats::lm.fit(Xc, y)
    bvec <- if (ncol(N) > 0) as.numeric(N %*% fit$coefficients[-1]) else rep(0, 5)
    return(list(
      name = spec$name,
      b = tibble::tibble(
        b0 = fit$coefficients[[1]], b1 = bvec[1], b2 = bvec[2],
        b3 = bvec[3], b4 = bvec[4], b5 = bvec[5]
      ),
      ssr = sum(fit$residuals^2), n = n,
      k_mean = 1L + ncol(N), converged = TRUE
    ))
  }

  # FPA-through-origin: any additional linear constraints must involve only the
  # quadratic coefficients (true for type 10's a2 = 0).
  C <- spec_constraint_matrix(spec)
  stopifnot(nrow(C) == 0 || all(abs(C[, 1:2]) < 1e-12))
  Nq <- .null_space(C[, 3:5, drop = FALSE], p = 3)

  profile_ssr <- function(theta) {
    q <- as.numeric(Nq %*% theta)
    b3 <- q[1]; b4 <- q[2]; b5 <- q[3]
    if (abs(b4) < 1e-8) {
      return(list(ssr = Inf))
    }
    p11 <- (b5 - b3 + sqrt((b3 - b5)^2 + b4^2)) / b4
    denom <- b4 + 2 * p11 * b5
    if (abs(denom) < 1e-8) {
      return(list(ssr = Inf))
    }
    kappa <- (p11 * b4 + 2 * b3) / denom
    Z <- cbind(1, kappa * data$P + data$E)
    f <- stats::lm.fit(Z, y - as.numeric(X[, 3:5, drop = FALSE] %*% q))
    list(
      ssr = sum(f$residuals^2),
      b0 = f$coefficients[[1]], b2 = f$coefficients[[2]], kappa = kappa
    )
  }

  b_full <- extract_b_coefficients(fit_quadratic_lm(data))
  q_full <- c(b_full$b3, b_full$b4, b_full$b5)
  starts <- list(as.numeric(crossprod(Nq, q_full)))
  if (abs(q_full[2]) < 0.05) {
    starts <- c(starts, lapply(c(0.1, -0.1), function(nudge) {
      as.numeric(crossprod(Nq, c(q_full[1], nudge, q_full[3])))
    }))
  }
  big <- sum((y - mean(y))^2) * 10
  objective <- function(th) {
    v <- profile_ssr(th)$ssr
    if (is.finite(v)) v else big
  }

  best <- NULL
  for (theta0 in starts) {
    if (!is.finite(profile_ssr(theta0)$ssr)) next
    opt <- tryCatch(
      stats::optim(theta0, objective,
        method = "Nelder-Mead",
        control = list(maxit = 2000, reltol = 1e-12)
      ),
      error = function(e) NULL
    )
    if (is.null(opt)) next
    sol <- profile_ssr(opt$par)
    if (!is.finite(sol$ssr)) next
    if (is.null(best) || sol$ssr < best$sol$ssr) {
      best <- list(opt = opt, sol = sol)
    }
  }

  k_mean <- 2L + ncol(Nq) # b0, b2, free quadratic directions
  if (is.null(best)) {
    return(list(
      name = spec$name,
      b = tibble::tibble(
        b0 = NA_real_, b1 = NA_real_, b2 = NA_real_,
        b3 = NA_real_, b4 = NA_real_, b5 = NA_real_
      ),
      ssr = NA_real_, n = n, k_mean = k_mean, converged = FALSE
    ))
  }
  sol <- best$sol
  q <- as.numeric(Nq %*% best$opt$par)
  list(
    name = spec$name,
    b = tibble::tibble(
      b0 = sol$b0, b1 = sol$kappa * sol$b2, b2 = sol$b2,
      b3 = q[1], b4 = q[2], b5 = q[3]
    ),
    ssr = sol$ssr, n = n, k_mean = k_mean,
    converged = best$opt$convergence == 0
  )
}

# ---- Information criteria ----------------------------------------------------

#' Gaussian AICc from a residual sum of squares
#' @param ssr Residual sum of squares.
#' @param n Sample size.
#' @param k_mean Free mean-structure parameters (residual variance added).
#' @return A numeric AICc.
#' @export
aicc_gaussian <- function(ssr, n, k_mean) {
  k <- k_mean + 1
  n * (log(2 * pi * ssr / n) + 1) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
}

#' Akaike weights from a vector of AICc values
#' @param aicc Numeric vector of AICc values.
#' @return Numeric vector of weights summing to 1.
#' @export
aicc_weights <- function(aicc) {
  d <- aicc - min(aicc, na.rm = TRUE)
  w <- exp(-d / 2)
  w / sum(w, na.rm = TRUE)
}

# ---- Candidate-set comparison (fast path) ------------------------------------

#' Compare a candidate set of constrained models by least squares
#'
#' Returns the comparison table: AICc (+ weights), R^2, and the F test of each
#' model's constraints against the full polynomial. `full` and `null` are
#' always included.
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param candidates Spec names to compare (default: all specs).
#' @param specs Spec registry (default [surface_model_specs]).
#' @return A tibble ordered by AICc; the fits are attached as `attr(., "fits")`.
#' @export
compare_surface_models_ls <- function(data, candidates = NULL,
                                      specs = surface_model_specs()) {
  candidates <- unique(c(candidates %||% names(specs), "full", "null"))
  stopifnot(all(candidates %in% names(specs)))
  fits <- lapply(specs[candidates], function(sp) fit_constrained_ls(data, sp))

  full <- fits[["full"]]
  tss <- sum((data$outcome - mean(data$outcome))^2)
  tab <- tibble::tibble(
    model = candidates,
    k_mean = vapply(fits, `[[`, 0L, "k_mean"),
    converged = vapply(fits, `[[`, TRUE, "converged"),
    ssr = vapply(fits, `[[`, 0, "ssr"),
    r2 = 1 - ssr / tss,
    aicc = aicc_gaussian(ssr, full$n, k_mean),
    df_vs_full = full$k_mean - k_mean,
    p_vs_full = purrr::map2_dbl(ssr, df_vs_full, function(s, df1) {
      if (df1 <= 0) {
        return(NA_real_)
      }
      f <- ((s - full$ssr) / df1) / (full$ssr / (full$n - full$k_mean))
      stats::pf(f, df1, full$n - full$k_mean, lower.tail = FALSE)
    })
  )
  tab$aicc[!tab$converged] <- NA_real_
  tab$delta_aicc <- tab$aicc - min(tab$aicc, na.rm = TRUE)
  tab$weight <- aicc_weights(tab$aicc)
  tab <- tab[order(tab$aicc), ]
  attr(tab, "fits") <- fits
  tab
}

# ---- Candidate-set comparison (lavaan / RSA showcase path) -------------------

#' Fit a candidate set through RSA/lavaan (MLR)
#'
#' Named models share one RSA() call; each custom constraint set gets its own
#' call (RSA appends `add` to every model, so custom constraints cannot ride
#' along with `full`).
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param candidates Spec names to fit.
#' @param specs Spec registry (default [surface_model_specs]).
#' @param ... Passed to [fit_rsa].
#' @return `list(models = named lavaan fits, rsa = named RSA objects)`.
#' @export
fit_candidates_lavaan <- function(data, candidates,
                                  specs = surface_model_specs(), ...) {
  candidates <- unique(c(candidates, "full", "null"))
  stopifnot(all(candidates %in% names(specs)))
  cand_specs <- specs[candidates]
  named <- vapply(cand_specs, function(sp) !is.na(sp$rsa_named), TRUE)

  base_models <- unique(c(
    vapply(cand_specs[named], `[[`, "", "rsa_named"), "full", "null"
  ))
  base_fit <- fit_rsa(data, models = base_models, ...)

  models <- list()
  rsa <- list()
  for (nm in candidates) {
    sp <- cand_specs[[nm]]
    if (!is.na(sp$rsa_named)) {
      models[[nm]] <- base_fit$models[[sp$rsa_named]]
      rsa[[nm]] <- base_fit
    } else {
      cf <- fit_rsa(data, models = "full", add = spec_lavaan_add(sp), ...)
      models[[nm]] <- cf$models$full
      rsa[[nm]] <- cf
    }
  }
  list(models = models, rsa = rsa)
}

.lavaan_aicc <- function(m) {
  fm <- lavaan::fitMeasures(m, c("aic", "npar", "ntotal"))
  k <- fm[["npar"]]
  n <- fm[["ntotal"]]
  fm[["aic"]] + 2 * k * (k + 1) / (n - k - 1)
}

.lavaan_r2 <- function(m, dv = "outcome") {
  r2 <- tryCatch(lavaan::lavInspect(m, "rsquare"), error = function(e) NULL)
  if (is.null(r2) || !dv %in% names(r2)) {
    return(NA_real_)
  }
  as.numeric(r2[[dv]])
}

# Scaled (Satorra-Bentler) chi-square difference test against the full model.
.lavaan_p_vs_full <- function(m, m_full) {
  out <- tryCatch(
    lavaan::lavTestLRT(m, m_full),
    error = function(e) NULL, warning = function(w) {
      suppressWarnings(lavaan::lavTestLRT(m, m_full))
    }
  )
  if (is.null(out) || nrow(out) < 2) {
    return(list(chisq = NA_real_, df = NA_real_, p = NA_real_))
  }
  list(
    chisq = out[2, "Chisq diff"], df = out[2, "Df diff"],
    p = out[2, "Pr(>Chisq)"]
  )
}

#' Compare a candidate set through RSA/lavaan
#'
#' @inheritParams fit_candidates_lavaan
#' @return A tibble ordered by AICc; fits attached as `attr(., "fits")`.
#' @export
compare_surface_models_lavaan <- function(data, candidates = NULL,
                                          specs = surface_model_specs(), ...) {
  candidates <- unique(c(candidates %||% names(specs), "full", "null"))
  fits <- fit_candidates_lavaan(data, candidates, specs = specs, ...)
  m_full <- fits$models[["full"]]

  rows <- purrr::map(candidates, function(nm) {
    m <- fits$models[[nm]]
    lrt <- if (nm == "full") {
      list(chisq = NA_real_, df = NA_real_, p = NA_real_)
    } else {
      .lavaan_p_vs_full(m, m_full)
    }
    tibble::tibble(
      model = nm,
      npar = as.numeric(lavaan::fitMeasures(m, "npar")),
      converged = lavaan::lavInspect(m, "converged"),
      r2 = .lavaan_r2(m),
      aicc = .lavaan_aicc(m),
      chisq_vs_full = lrt$chisq, df_vs_full = lrt$df, p_vs_full = lrt$p
    )
  })
  tab <- dplyr::bind_rows(rows)
  tab$aicc[!tab$converged] <- NA_real_
  tab$delta_aicc <- tab$aicc - min(tab$aicc, na.rm = TRUE)
  tab$weight <- aicc_weights(tab$aicc)
  tab <- tab[order(tab$aicc), ]
  attr(tab, "fits") <- fits
  tab
}

# ---- Decision workflows ------------------------------------------------------

#' Default candidate set for a theory model
#'
#' The theory itself, its nested simplifications from the type family, `full`,
#' and `null`.
#'
#' @param theory A spec name.
#' @param specs Spec registry (default [surface_model_specs]).
#' @return A character vector of spec names.
#' @export
default_candidate_set <- function(theory, specs = surface_model_specs()) {
  th <- specs[[theory]]
  nested <- names(specs)[vapply(specs, function(sp) {
    grepl("^type_", sp$name) && sp$name != theory && spec_is_nested(sp, th)
  }, TRUE)]
  unique(c(theory, nested, "full", "null"))
}

#' Three-move constrained-model test for one theory-implied surface
#'
#' 1. constraint tenability (chi-square/F difference vs. full);
#' 2. candidate selection (AICc weights across the candidate set);
#' 3. predictive value (R^2 of the retained model vs. null).
#' The retained model is the AICc-best among candidates whose constraints are
#' tenable (`p_vs_full > alpha`; the full model is tenable by definition).
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param theory A spec name, or a one-row data frame with `theory_spec` and
#'   `valence` columns (negative-valence outcomes are reverse-coded before
#'   fitting).
#' @param rivals Extra spec names to add to the candidate set.
#' @param method `"lavaan"` or `"ls"`.
#' @param alpha Tenability threshold.
#' @param specs Spec registry (default [surface_model_specs]).
#' @param ... Passed to the comparison function.
#' @return A list describing the decision (retained model, tenability, R^2
#'   gain, the full decision table, and the fits).
#' @export
test_surface_theory <- function(data, theory, rivals = NULL,
                                method = c("lavaan", "ls"),
                                alpha = 0.05,
                                specs = surface_model_specs(), ...) {
  method <- match.arg(method)

  valence_flipped <- FALSE
  if (is.data.frame(theory)) {
    stopifnot(nrow(theory) == 1, "theory_spec" %in% names(theory))
    if (identical(theory$valence, "negative")) {
      data$outcome <- -data$outcome
      valence_flipped <- TRUE
    }
    theory <- theory$theory_spec
  }
  stopifnot(is.character(theory), theory %in% names(specs))
  if (!is.null(rivals)) stopifnot(all(rivals %in% names(specs)))

  candidates <- unique(c(default_candidate_set(theory, specs), rivals))
  tab <- switch(method,
    lavaan = compare_surface_models_lavaan(data, candidates, specs = specs, ...),
    ls = compare_surface_models_ls(data, candidates, specs = specs)
  )

  tab$tenable <- is.na(tab$p_vs_full) | tab$p_vs_full > alpha
  eligible <- tab[tab$tenable & !is.na(tab$aicc), ]
  retained <- if (nrow(eligible) > 0) eligible$model[[1]] else "full"

  r2_null <- tab$r2[tab$model == "null"]
  list(
    theory = theory,
    candidates = candidates,
    decision = tab,
    retained = retained,
    theory_tenable = tab$tenable[tab$model == theory],
    theory_retained = identical(retained, theory),
    r2_gain_vs_null = tab$r2[tab$model == retained] - r2_null,
    alpha = alpha,
    method = method,
    valence_flipped = valence_flipped,
    fits = attr(tab, "fits")
  )
}

#' Classify a dataset into a Yao & Ma type via constrained model selection
#'
#' The constrained-route counterpart to [classify_yao_ma]: instead of reading
#' the type off one full fit, this fits all 16 Table-4 type models (types 1-15
#' plus `full` = type 16), keeps those whose constraints are tenable against
#' the full polynomial, and returns the declared type of the AICc-best tenable
#' model. This is a convenience wrapper not present in the original rsasim
#' `lib/`; it composes [compare_surface_models_ls] /
#' [compare_surface_models_lavaan] with the spec `yao_ma_type` map.
#'
#' @param data A data frame with `P`, `E`, `outcome`.
#' @param valence `"positive"` or `"negative"` (negative outcomes are
#'   reverse-coded before fitting, matching [test_surface_theory]).
#' @param method `"ls"` (fast, default) or `"lavaan"`.
#' @param alpha Tenability threshold for the constraint tests.
#' @param specs Spec registry (default [surface_model_specs]).
#' @param ... Passed to the comparison function.
#' @return A list with `type` (integer Yao & Ma type of the retained model),
#'   `retained` (spec name), `label`, and `decision` (the ranked table).
#' @export
classify_yao_ma_constrained <- function(data, valence = "positive",
                                        method = c("ls", "lavaan"),
                                        alpha = 0.05,
                                        specs = surface_model_specs(), ...) {
  method <- match.arg(method)
  if (identical(valence, "negative")) data$outcome <- -data$outcome

  candidates <- unique(c(
    names(specs)[grepl("^type_", names(specs))], "full", "null"
  ))
  tab <- switch(method,
    ls = compare_surface_models_ls(data, candidates, specs = specs),
    lavaan = compare_surface_models_lavaan(data, candidates, specs = specs, ...)
  )

  tab$tenable <- is.na(tab$p_vs_full) | tab$p_vs_full > alpha
  eligible <- tab[tab$tenable & !is.na(tab$aicc) & tab$model != "null", ]
  retained <- if (nrow(eligible) > 0) eligible$model[[1]] else "full"

  list(
    type = specs[[retained]]$yao_ma_type,
    retained = retained,
    label = specs[[retained]]$label,
    decision = tab
  )
}
