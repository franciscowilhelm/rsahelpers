#' Monte Carlo intervals for response surface parameters
#'
#' Draws the polynomial weights `b1` to `b5` of a fitted Mplus response
#' surface from their asymptotic sampling distribution (a multivariate normal
#' with the estimates as means and the `TECH3` covariance matrix) and computes
#' the surface parameters on each draw. This gives intervals for parameters
#' whose delta-method standard errors are poor, such as ratios, without
#' refitting the model, and draws on which any further function of the surface
#' can be evaluated.
#'
#' @inheritParams RSA_mplus
#' @param n Number of draws.
#' @param probs Lower and upper quantile of the reported intervals.
#' @param seed Optional seed, set locally for the draws.
#'
#' @details
#' The model must have been estimated by maximum likelihood with `TECH1` and
#' `TECH3` in the `OUTPUT` command, for example
#' `blocks = list(OUTPUT = "STDYX CINTERVAL TECH1 TECH3;")` in
#' [create_rsa_mplus_model()]. Bayesian models have no `TECH3`; use their
#' posterior draws (`BPARAMETERS`) instead.
#'
#' Terms that share a parameter in the model, for example through a shared
#' label, share its draws. The restrictions of `rsa_model` are imposed on the
#' draws and on the estimates: terms the model sets to zero are 0, and terms
#' it defines from others (such as `b4 = -2*b5` in `"SRR"`) are computed from
#' them. This also removes the rounding of the printed estimates, which do not
#' satisfy the restrictions exactly.
#'
#' Besides `b1` to `b5`, the draws hold the surface parameters `a1` to `a5`
#' (`cs`, `cc`, `is`, `ic`, and `a5` in the generated `MODEL CONSTRAINT`),
#' and the parameters specific to `rsa_model` in [create_rsa_mplus_model()]:
#' the ridge parameters `c`, `s`, `bm`, and `a4r`, or the stationary point and
#' principal axes. Ratios can be undefined or extreme in single draws, so
#' their intervals are quantiles and the standard deviation in `se` may be
#' large.
#'
#' Mplus prints estimates and covariances with limited precision. The draws
#' are therefore centred on estimates rounded to three decimals.
#'
#' @return An object of class `rsa_mplus_draws`, a list with
#'   * `draws`: a data frame with one row per draw.
#'   * `summary`: a data frame with the `parameter`, its estimate `est`
#'     computed from the point estimates, the standard deviation `se` of the
#'     draws, and the `lower` and `upper` quantiles.
#'   * `rsa_model`, `n`, and `probs`.
#' @seealso [RSA_mplus()], [rsa_mplus_supported_range()]
#' @export
#'
#' @examples
#' \dontrun{
#' fitted <- rsa_mplus_workflow(
#'   data = dat,
#'   measurement = measurement,
#'   rsa_model = "SRR",
#'   blocks = list(OUTPUT = "TECH1 TECH3;"),
#'   modelout = "rsa_SRR.inp",
#'   run = TRUE
#' )
#' intervals <- rsa_mplus_draws(fitted, seed = 1)
#' intervals$summary
#' }
rsa_mplus_draws <- function(
  model,
  outcome = NULL,
  pred_x = NULL,
  pred_y = NULL,
  pred_x2 = NULL,
  pred_xy = NULL,
  pred_y2 = NULL,
  rsa_model = NULL,
  n = 10000L,
  probs = c(0.025, 0.975),
  seed = NULL
) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    rlang::abort("Package `MASS` must be installed to draw surface parameters.")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    rlang::abort("`n` must be a single positive number.")
  }
  if (
    !is.numeric(probs) ||
      length(probs) != 2L ||
      anyNA(probs) ||
      any(probs < 0 | probs > 1) ||
      probs[[1]] >= probs[[2]]
  ) {
    rlang::abort("`probs` must give a lower and an upper quantile.")
  }

  surface <- RSA_mplus(
    model,
    outcome = outcome,
    pred_x = pred_x,
    pred_y = pred_y,
    pred_x2 = pred_x2,
    pred_xy = pred_xy,
    pred_y2 = pred_y2,
    rsa_model = rsa_model,
    b0 = 0,
    include_new = FALSE,
    plot = FALSE
  )
  rsa_model <- surface$rsa_model
  results <- rsa_mplus_results(surface$model)

  specification <- results$tech1$parameterSpecification
  if (!is.null(specification) && is.null(specification$beta)) {
    specification <- specification[[1]]
  }
  covariance <- results$tech3$paramCov
  if (is.null(specification) || is.null(covariance)) {
    rlang::abort(c(
      "The Mplus output has no TECH1 or TECH3 section.",
      i = "Request them with `OUTPUT: TECH1 TECH3;`, e.g. `blocks = list(OUTPUT = \"TECH1 TECH3;\")`.",
      i = "Bayesian models have no TECH3; use their posterior draws instead."
    ))
  }
  covariance <- as.matrix(covariance)
  upper <- upper.tri(covariance)
  covariance[upper] <- t(covariance)[upper]

  terms <- names(rsa_mplus_coefficient_labels)
  labels <- unname(rsa_mplus_coefficient_labels)
  estimates <- stats::setNames(unname(surface$coefficients[terms]), labels)

  # Parameter number of each term in TECH1; 0 for fixed or absent terms.
  matched <- surface$regression_parameters
  predictors <- stats::setNames(
    sub("^.*<-", "", matched$Label),
    rsa_mplus_coefficient_labels[matched$term]
  )
  numbers <- stats::setNames(integer(length(labels)), labels)
  for (label in names(predictors)) {
    numbers[[label]] <- rsa_mplus_parameter_number(
      specification,
      outcome = surface$outcome,
      predictor = predictors[[label]]
    )
  }

  spec <- if (rsa_model %in% names(rsa_mplus_model_restrictions)) {
    rsa_mplus_model_spec(rsa_model)
  } else {
    list(zero = character(), linear = list(), nonlinear = character())
  }
  dependent <- c(spec$zero, names(spec$linear), names(spec$nonlinear))
  drawn <- setdiff(labels[numbers > 0L], dependent)
  ids <- unique(numbers[drawn])
  if (any(ids > nrow(covariance)) || anyNA(covariance[ids, ids])) {
    rlang::abort(
      "TECH3 does not hold the covariances of the response surface terms."
    )
  }

  n <- as.integer(n)
  draws <- matrix(
    rep(estimates, each = n),
    nrow = n,
    dimnames = list(NULL, labels)
  )
  if (length(ids) > 0L) {
    first <- drawn[match(ids, numbers[drawn])]
    sample_draws <- function() {
      MASS::mvrnorm(
        n,
        mu = estimates[first],
        Sigma = covariance[ids, ids, drop = FALSE]
      )
    }
    sampled <- if (is.null(seed)) {
      sample_draws()
    } else {
      rsa_mplus_with_seed(seed, sample_draws())
    }
    sampled <- matrix(sampled, nrow = n)
    for (label in drawn) {
      draws[, label] <- sampled[, match(numbers[[label]], ids)]
    }
  }

  draws <- rsa_mplus_derive_surface(as.data.frame(draws), spec)
  point <- rsa_mplus_derive_surface(as.data.frame(as.list(estimates)), spec)

  finite <- function(x) x[is.finite(x)]
  summary <- data.frame(
    parameter = names(draws),
    est = unlist(point[1L, ], use.names = FALSE),
    se = vapply(draws, function(x) stats::sd(finite(x)), numeric(1)),
    lower = vapply(
      draws,
      function(x) stats::quantile(finite(x), probs[[1]], names = FALSE),
      numeric(1)
    ),
    upper = vapply(
      draws,
      function(x) stats::quantile(finite(x), probs[[2]], names = FALSE),
      numeric(1)
    )
  )
  rownames(summary) <- NULL

  structure(
    list(
      draws = draws,
      summary = summary,
      rsa_model = rsa_model,
      n = n,
      probs = probs
    ),
    class = "rsa_mplus_draws"
  )
}

#' @export
print.rsa_mplus_draws <- function(x, digits = 3, ...) {
  cat("<rsa_mplus_draws>\n")
  cat(sprintf(
    "RSA model: %s; %d draws; quantiles %s and %s\n",
    x$rsa_model,
    x$n,
    format(x$probs[[1]]),
    format(x$probs[[2]])
  ))
  shown <- x$summary
  numeric_columns <- vapply(shown, is.numeric, logical(1))
  shown[numeric_columns] <- lapply(shown[numeric_columns], round, digits)
  print(shown, row.names = FALSE)
  invisible(x)
}

#' Supported region of a latent response surface
#'
#' Gives the range of the line of congruence and the line of incongruence
#' that the data support, from the model-implied distribution of the latent
#' X and Y factors. Factor scores are shrunk towards the mean, so their range
#' is narrower than the latent distribution to which the surface refers; this
#' function uses the estimated latent means, variances, and covariance
#' instead and assumes bivariate normality.
#'
#' @param model A fitted model as in [RSA_mplus()] (a workflow, an Mplus
#'   output path, or a parsed Mplus model), or a list of them. The moments of
#'   several models are averaged, so that surfaces of the same X and Y factors
#'   for different outcomes share one region.
#' @param pred_x,pred_y Names of the X and Y factors. Taken from a workflow
#'   when `NULL`.
#' @param probs Lower and upper quantile that bound the region.
#'
#' @details
#' A point `(X, Y)` is at position `(X + Y)/2` along the line of congruence
#' and `(X - Y)/2` along the line of incongruence, the coordinates in which
#' `a1` to `a4` describe the surface. Their means and variances follow from
#' the latent moments. A factor mean that the model does not estimate is 0,
#' the Mplus default.
#'
#' @return A list with the quantiles `loc` and `loic`, and `moments`, the
#'   (averaged) latent means `mx` and `my`, variances `vx` and `vy`, and
#'   covariance `cxy`.
#' @seealso [rsa_mplus_draws()]
#' @export
#'
#' @examples
#' \dontrun{
#' rsa_mplus_supported_range(fitted)
#' rsa_mplus_supported_range("model.out", pred_x = "X", pred_y = "Y")
#' }
rsa_mplus_supported_range <- function(
  model,
  pred_x = NULL,
  pred_y = NULL,
  probs = c(0.05, 0.95)
) {
  if (
    !is.numeric(probs) ||
      length(probs) != 2L ||
      anyNA(probs) ||
      any(probs <= 0 | probs >= 1)
  ) {
    rlang::abort("`probs` must give a lower and an upper quantile.")
  }
  single <- is.character(model) ||
    inherits(model, c("rsa_mplus_workflow", "mplus.model", "mplusObject"))
  models <- if (single) list(model) else model
  if (!is.list(models) || length(models) == 0L) {
    rlang::abort("`model` must be a fitted model or a list of fitted models.")
  }

  moments <- lapply(models, function(one) {
    x <- pred_x
    y <- pred_y
    if (inherits(one, "rsa_mplus_workflow")) {
      if (is.null(x)) x <- one$spec$metadata$pred_x
      if (is.null(y)) y <- one$spec$metadata$pred_y
    }
    if (is.null(x) || is.null(y)) {
      rlang::abort("`pred_x` and `pred_y` are needed to find the latent factors.")
    }
    rsa_mplus_latent_moments(rsa_mplus_results(one), x, y)
  })
  moments <- colMeans(do.call(rbind, moments))

  loc_mean <- (moments[["mx"]] + moments[["my"]]) / 2
  loic_mean <- (moments[["mx"]] - moments[["my"]]) / 2
  total <- moments[["vx"]] + moments[["vy"]]
  loc_sd <- sqrt((total + 2 * moments[["cxy"]]) / 4)
  loic_sd <- sqrt((total - 2 * moments[["cxy"]]) / 4)
  list(
    loc = stats::qnorm(probs, loc_mean, loc_sd),
    loic = stats::qnorm(probs, loic_mean, loic_sd),
    moments = moments
  )
}

# Parsed Mplus results of a workflow, a parsed model, or an output path.
rsa_mplus_results <- function(model) {
  if (inherits(model, "rsa_mplus_workflow")) {
    results <- model$results
    if (is.null(results)) results <- model$mplus$results
    if (is.null(results)) {
      rlang::abort("The RSA-Mplus workflow does not contain fitted results.")
    }
    return(results)
  }
  if (inherits(model, "mplusObject")) {
    return(model$results)
  }
  if (inherits(model, "mplus.model")) {
    return(model)
  }
  if (is.character(model) && length(model) == 1L && !is.na(model)) {
    require_rsa_mplus_namespace("MplusAutomation")
    return(MplusAutomation::readModels(target = model, quiet = TRUE))
  }
  rlang::abort(
    "`model` must be a file path, Mplus model object, or RSA-Mplus workflow."
  )
}

# Number of the regression of `outcome` on `predictor` in the TECH1 parameter
# specification: BETA for latent predictors, GAMMA for observed ones.
rsa_mplus_parameter_number <- function(specification, outcome, predictor) {
  for (name in c("beta", "gamma")) {
    matrix <- specification[[name]]
    if (is.null(matrix)) {
      next
    }
    row <- match(toupper(outcome), toupper(rownames(matrix)))
    column <- match(toupper(predictor), toupper(colnames(matrix)))
    if (!is.na(row) && !is.na(column)) {
      number <- suppressWarnings(as.integer(matrix[row, column]))
      if (!is.na(number) && number > 0L) {
        return(number)
      }
    }
  }
  0L
}

# Imposes the restrictions of a model on b1-b5 and adds the surface parameters
# a1-a5 and the model's ridge or principal-axis parameters, evaluating the
# same expressions that are written to MODEL CONSTRAINT.
rsa_mplus_derive_surface <- function(b, spec) {
  evaluate <- function(expression, values) {
    expression <- gsub("SQRT", "sqrt", expression, fixed = TRUE)
    suppressWarnings(eval(str2lang(expression), values))
  }
  b[spec$zero] <- 0
  for (term in names(spec$linear)) {
    weights <- spec$linear[[term]]
    b[[term]] <- as.vector(as.matrix(b[names(weights)]) %*% weights)
  }
  for (term in names(spec$nonlinear)) {
    b[[term]] <- evaluate(spec$nonlinear[[term]], b)
  }

  surface_names <- c(cs = "a1", cc = "a2", is = "a3", ic = "a4", a5 = "a5")
  for (name in names(rsa_mplus_surface_weights)) {
    weights <- rsa_mplus_surface_weights[[name]]
    b[[surface_names[[name]]]] <- as.vector(
      as.matrix(b[names(weights)]) %*% weights
    )
  }
  # p10 and p20 refer to x0, y0, and the slopes, so the order matters.
  extra <- c(spec$principal, spec$derived)
  for (name in names(extra)) {
    b[[name]] <- evaluate(extra[[name]], b)
  }
  b
}

rsa_mplus_latent_moments <- function(results, x, y) {
  parameters <- results$parameters$unstandardized
  if (is.null(parameters)) {
    rlang::abort("The Mplus output has no unstandardized parameter estimates.")
  }
  get <- function(header, param, default = NA_real_) {
    value <- parameters$est[
      toupper(parameters$paramHeader) == toupper(header) &
        toupper(parameters$param) == toupper(param)
    ]
    if (length(value) == 1L) as.numeric(value) else default
  }
  covariance <- get(paste0(y, ".WITH"), x)
  if (is.na(covariance)) {
    covariance <- get(paste0(x, ".WITH"), y)
  }
  moments <- c(
    mx = get("Means", x, 0),
    my = get("Means", y, 0),
    vx = get("Variances", x),
    vy = get("Variances", y),
    cxy = covariance
  )
  if (anyNA(moments)) {
    rlang::abort(sprintf(
      "Could not find the latent variances and covariance of `%s` and `%s`.",
      x,
      y
    ))
  }
  moments
}

# Evaluates `expr` with a seed and restores the state of the generator.
rsa_mplus_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  )
  set.seed(seed)
  expr
}
