# Parsed results in the shape of MplusAutomation::readModels().
mock_rsa_results <- function(
  LL,
  k,
  scaling = 1,
  n = 400,
  r2 = NA,
  parameters = NULL,
  warnings = list(),
  errors = list()
) {
  summaries <- data.frame(Observations = n, Parameters = k)
  if (!is.na(LL)) {
    summaries$LL <- LL
    summaries$LLCorrectionFactor <- scaling
  }
  list(
    summaries = summaries,
    parameters = list(
      unstandardized = parameters,
      r2 = if (!is.na(r2)) data.frame(param = c("Z1", "Z"), est = c(0.8, r2))
    ),
    warnings = warnings,
    errors = errors
  )
}

mock_rsa_workflow <- function(rsa_model, results) {
  new_rsa_mplus_workflow(
    spec = list(
      rsa_model = rsa_model,
      roles = c(x = "X", y = "Y", outcome = "Z"),
      metadata = list(pred_y2 = "YS")
    ),
    results = results,
    status = "run"
  )
}

# Fits of the rising-ridge data of the opt-in Mplus test (Mplus 8.10).
rising_ridge_fits <- function() {
  data.frame(
    model = c("full", "SRRR", "SRR", "SRSQD", "IA", "null"),
    n = 400,
    k = c(33, 32, 31, 31, 31, 28),
    LL = c(-3533.771, -3534.260, -3542.259, -3550.395, -3639.773, -3670.066),
    scaling = c(0.9875, 0.9914, 0.9972, 0.9959, 1.0993, 1.0230),
    R2 = c(0.600, 0.596, 0.585, 0.551, 0.156, 0.000),
    b5 = c(NA, -0.273, -0.362, -0.288, NA, NA),
    ridge_distance = c(NA, 0.166, 0.156, 0.140, NA, NA),
    converged = TRUE,
    status = "ok"
  )
}

no_messages <- function() {
  data.frame(model = character(), type = character(), message = character())
}

congruence_sim_path <- function() {
  path <- system.file("extdata", "congruence_sim.out", package = "rsahelpers")
  if (!nzchar(path)) {
    path <- normalizePath(
      file.path("inst", "extdata", "congruence_sim.out"),
      mustWork = TRUE
    )
  }
  path
}

test_that("the nesting relations follow the generated restrictions", {
  set.seed(11)
  labels <- unname(rsa_mplus_coefficient_labels)
  evaluate <- function(expressions, b) {
    vapply(expressions, function(x) eval(str2lang(x), as.list(b)), numeric(1))
  }
  draw <- function(spec) {
    b <- stats::setNames(numeric(5), labels)
    b[spec$free] <- stats::runif(length(spec$free), 0.2, 1) *
      sample(c(-1, 1), length(spec$free), replace = TRUE)
    for (term in names(spec$linear)) {
      weights <- spec$linear[[term]]
      b[[term]] <- sum(weights * b[names(weights)])
    }
    b[names(spec$nonlinear)] <- evaluate(spec$nonlinear, b)
    b
  }
  satisfies <- function(b, spec) {
    linear <- vapply(
      spec$linear,
      function(weights) sum(weights * b[names(weights)]),
      numeric(1)
    )
    restricted <- c(
      stats::setNames(numeric(length(spec$zero)), spec$zero),
      linear,
      evaluate(spec$nonlinear, b)
    )
    isTRUE(all.equal(unname(b[names(restricted)]), unname(restricted)))
  }

  models <- names(rsa_mplus_model_restrictions)
  expect_setequal(names(rsa_mplus_model_parents), models)
  for (child in models) {
    child_spec <- rsa_mplus_model_spec(child)
    for (parent in rsa_mplus_model_parents[[child]]) {
      parent_spec <- rsa_mplus_model_spec(parent)
      expect_lt(length(child_spec$free), length(parent_spec$free))
      for (i in 1:5) {
        expect_true(
          satisfies(draw(child_spec), parent_spec),
          label = paste(child, "in", parent)
        )
      }
    }
  }

  expect_length(rsa_mplus_model_ancestors("full"), 0L)
  expect_setequal(rsa_mplus_model_ancestors("null"), setdiff(models, "null"))
  expect_setequal(
    rsa_mplus_model_ancestors("SQD"),
    c("RR", "SRR", "SSQD", "SRSQD", "SRRR", "a5zero", "full")
  )
  expect_true(all(
    c("addx2", "addy2", "noIA") %in% rsa_mplus_model_ancestors("additive")
  ))
  expect_false("a5zero" %in% rsa_mplus_model_ancestors("noIA"))
  expect_false("SRRR" %in% rsa_mplus_model_ancestors("onlyx2"))
  expect_true("SRSQD" %in% rsa_mplus_model_ancestors("onlyy2"))
})

test_that("models are ranked by AICc and tested against the full model", {
  fits <- rising_ridge_fits()
  comparison <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)

  expect_s3_class(comparison, "rsa_mplus_comparison")
  expect_identical(
    comparison$model,
    c("SRRR", "full", "SRR", "SRSQD", "IA", "null")
  )
  aicc <- -2 * fits$LL + 2 * fits$k + 2 * fits$k * (fits$k + 1) /
    (fits$n - fits$k - 1)
  expect_equal(comparison$AICc, aicc[match(comparison$model, fits$model)])
  expect_equal(comparison$AICc[[2]], 7139.673, tolerance = 1e-6)
  expect_equal(comparison$delta_AICc, comparison$AICc - min(aicc))
  expect_equal(sum(comparison$weight), 1)
  expect_equal(comparison$cum_weight, cumsum(comparison$weight))
  expect_equal(
    comparison$weight[[1]] / comparison$weight[-1],
    comparison$evidence_ratio[-1]
  )
  expect_true(is.na(comparison$evidence_ratio[[1]]))

  # SRRR against full by hand: cd = (32*0.9914 - 33*0.9875)/(32 - 33).
  srrr <- comparison[comparison$model == "SRRR", ]
  expect_equal(srrr$lrt_df, 1)
  expect_equal(srrr$lrt_cd, 0.8627)
  expect_equal(srrr$lrt_chisq, 2 * 0.489 / 0.8627)
  expect_equal(srrr$lrt_p, stats::pchisq(2 * 0.489 / 0.8627, 1, lower.tail = FALSE))
  full <- comparison[comparison$model == "full", ]
  expect_true(all(is.na(unlist(full[c("lrt_df", "lrt_cd", "lrt_chisq", "lrt_p")]))))

  # IA has a larger scaling factor than the full model, so cd is negative.
  ia <- comparison[comparison$model == "IA", ]
  expect_equal(ia$lrt_cd, (31 * 1.0993 - 33 * 0.9875) / -2)
  expect_true(is.na(ia$lrt_chisq) && is.na(ia$lrt_p))
  expect_identical(ia$flags, "cd <= 0")
  expect_identical(comparison$flags[comparison$model != "IA"], character(5))

  expect_equal(comparison$R2, fits$R2[match(comparison$model, fits$model)])
  expect_true(all(comparison$converged))
  expect_snapshot(print(comparison))

  # Without a scaling correction factor (ESTIMATOR = ML), c = 1.
  fits$scaling <- NA_real_
  ml <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)
  expect_equal(ml$lrt_cd[ml$model == "SRRR"], 1)
  expect_equal(ml$lrt_chisq[ml$model == "IA"], 2 * (3639.773 - 3533.771))
})

test_that("likelihood-ratio tests can use another reference model", {
  fits <- rising_ridge_fits()
  comparison <- rsa_mplus_comparison_table(
    fits,
    no_messages(),
    0.01,
    3,
    reference = "SRRR"
  )

  # SRR against SRRR by hand: cd = (31*0.9972 - 32*0.9914)/(31 - 32).
  srr <- comparison[comparison$model == "SRR", ]
  cd <- 32 * 0.9914 - 31 * 0.9972
  expect_equal(srr$lrt_df, 1)
  expect_equal(srr$lrt_cd, cd)
  expect_equal(srr$lrt_chisq, 2 * (3542.259 - 3534.260) / cd)

  # Only the models nested in the reference are tested.
  tested <- comparison$model[!is.na(comparison$lrt_df)]
  expect_setequal(tested, c("SRR", "SRSQD", "null"))
  # The ranking does not depend on the reference.
  expect_identical(
    comparison$model,
    rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)$model
  )

  workflows <- list(
    full = mock_rsa_workflow("full", mock_rsa_results(-100, 10)),
    IA = mock_rsa_workflow("IA", mock_rsa_results(-101, 8))
  )
  expect_error(
    compare_rsa_mplus_models(workflows, reference = "a5zero"),
    "reference model `a5zero` is not among the fitted models",
    fixed = TRUE
  )
  expect_error(
    compare_rsa_mplus_models(workflows, reference = "cubic"),
    "`reference` must name one model",
    fixed = TRUE
  )
})

test_that("nesting violations and lines outside the data are flagged", {
  fits <- rising_ridge_fits()
  fits$LL[fits$model == "SRSQD"] <- -3534.21
  expect_snapshot(
    comparison <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)
  )
  expect_identical(
    comparison$flags[comparison$model == "SRSQD"],
    "LL above SRRR by 0.050"
  )
  expect_no_warning(rsa_mplus_comparison_table(fits, no_messages(), 0.1, 3))

  # A rotated ridge model that drifted towards b5 = 0 on a plane.
  fits <- rising_ridge_fits()
  fits$b5[fits$model == "SRSQD"] <- 0.0002
  fits$ridge_distance[fits$model == "SRSQD"] <- 9.88
  fits$b5[fits$model == "SRR"] <- 0.0001
  fits$ridge_distance[fits$model == "SRRR"] <- 2.9
  comparison <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)
  flags <- stats::setNames(comparison$flags, comparison$model)
  expect_identical(flags[["SRSQD"]], "b5 = 0; line 9.9 SD from the data")
  expect_identical(flags[["SRR"]], "")
  expect_identical(flags[["SRRR"]], "")
  comparison <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 2)
  expect_identical(
    comparison$flags[comparison$model == "SRRR"],
    "line 2.9 SD from the data"
  )
})

test_that("models without a log-likelihood are listed last and not ranked", {
  fits <- rising_ridge_fits()
  fits$LL[fits$model == "SRR"] <- NA_real_
  fits$converged[fits$model == "SRR"] <- FALSE
  fits$status[fits$model == "SRR"] <- "not converged"
  messages <- data.frame(
    model = c("SRR", "IA", "IA"),
    type = c("error", "warning", "warning"),
    message = c("THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY.", "a", "b")
  )
  comparison <- rsa_mplus_comparison_table(fits, messages, 0.01, 3)

  expect_identical(comparison$model[[6]], "SRR")
  last <- comparison[6, ]
  expect_true(all(is.na(unlist(last[c("AICc", "weight", "cum_weight", "lrt_p")]))))
  expect_false(last$converged)
  expect_identical(last$flags, "not converged; 1 Mplus error")
  expect_identical(
    comparison$flags[comparison$model == "IA"],
    "2 Mplus warnings; cd <= 0"
  )
  expect_equal(sum(comparison$weight, na.rm = TRUE), 1)
  expect_identical(attr(comparison, "messages"), messages)

  fits <- rising_ridge_fits()
  fits$n[fits$model == "null"] <- 399
  expect_warning(
    rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3),
    "different numbers of observations"
  )
})

test_that("fit summaries read the parsed Mplus results", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(12)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  full <- read_rsa_mplus_model(
    create_rsa_mplus_model(
      dat,
      list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
    ),
    output = congruence_sim_path()
  )
  comparison <- compare_rsa_mplus_models(list(full))
  summaries <- full$results$summaries
  expect_identical(comparison$model, "full")
  expect_equal(comparison$k, summaries$Parameters)
  expect_equal(comparison$LL, summaries$LL)
  expect_equal(comparison$AICc, summaries$AICC, tolerance = 1e-6)
  expect_equal(comparison$R2, 0.32)
  expect_equal(comparison$weight, 1)

  # Ridge parameters, means, and variances of a drifted SRSQD fit.
  parameters <- data.frame(
    paramHeader = c(
      "Z.ON", "Z.ON", "Means", "Means", "Variances", "Variances",
      "New.Additional.Parameters", "New.Additional.Parameters"
    ),
    param = c("XY", "YS", "X", "Y", "X", "Y", "C", "S"),
    est = c(-0.006, 0, 0.2, -0.2, 0.9, 1.1, 72.787, -6.995)
  )
  srsqd <- mock_rsa_workflow(
    "SRSQD",
    mock_rsa_results(-3479.888, 31, 0.9587, r2 = 0.313, parameters = parameters)
  )
  summary <- rsa_mplus_fit_summary(srsqd, "SRSQD")$fit
  expect_equal(summary$b5, 0)
  expect_equal(summary$R2, 0.313)
  expect_equal(
    summary$ridge_distance,
    abs(-6.995 * 0.2 + 72.787 + 0.2) / sqrt(1 + 6.995^2) / sqrt(1)
  )
  expect_true(summary$converged)

  stopped <- mock_rsa_workflow(
    "SRSQD",
    mock_rsa_results(
      NA,
      31,
      warnings = list(c("WARNING:  THE BEST LOGLIKELIHOOD", "WAS NOT REPLICATED.")),
      errors = list("THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY.")
    )
  )
  stopped <- rsa_mplus_fit_summary(stopped, "SRSQD")
  expect_false(stopped$fit$converged)
  expect_identical(stopped$fit$status, "not converged")
  expect_identical(stopped$messages$type, c("error", "warning"))
  expect_identical(
    stopped$messages$message[[2]],
    "WARNING:  THE BEST LOGLIKELIHOOD WAS NOT REPLICATED."
  )

  unrun <- rsa_mplus_fit_summary(NULL, "SRR")$fit
  expect_identical(unrun$status, "not run")
  expect_false(unrun$converged)
})

test_that("comparison input is validated", {
  srrr <- mock_rsa_workflow("SRRR", mock_rsa_results(-3534.26, 32))
  expect_error(compare_rsa_mplus_models(srrr), "named list of fitted workflows")
  expect_error(compare_rsa_mplus_models(list()), "named list of fitted workflows")
  expect_error(
    compare_rsa_mplus_models(list(SRRR = srrr, other = srrr)),
    "`SRRR` is repeated",
    fixed = TRUE
  )
  unnamed <- new_rsa_mplus_workflow(results = mock_rsa_results(-3534.26, 32))
  expect_error(
    compare_rsa_mplus_models(list(unnamed)),
    "Name each workflow by its model"
  )
  expect_identical(
    compare_rsa_mplus_models(list(SRRR = unnamed))$model,
    "SRRR"
  )
  expect_error(compare_rsa_mplus_models(list(srrr), boundary = 0), "`boundary`")
})

test_that("the family runner continues after failures and starts from the full model", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(13)
  dat <- as.data.frame(matrix(stats::rnorm(700), ncol = 7))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2", "age")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  output_dir <- withr::local_tempdir()
  modelout <- file.path(output_dir, "rsa.inp")

  # The full model reads a real output, the additive model fails in Mplus,
  # and the label `c` on a covariate clashes with the ridge parameter C of
  # SRR.
  local_mocked_bindings(
    run_rsa_mplus_model = function(x, ...) {
      rsa_model <- x$spec$rsa_model
      if (rsa_model == "full") {
        x <- read_rsa_mplus_model(x, output = congruence_sim_path())
      } else if (rsa_model == "additive") {
        rlang::abort("Mplus reported an error.")
      } else {
        x$results <- mock_rsa_results(-5290, 31, n = 500)
      }
      x$status <- "run"
      x
    }
  )
  expect_snapshot(
    family <- fit_rsa_mplus_family(
      dat,
      measurement,
      structural = "Z ~ X + Y + c*age",
      models = c("IA", "additive", "full", "SRR"),
      modelout = modelout,
      Mplus_command = "mplus"
    )
  )

  expect_s3_class(family, "rsa_mplus_family")
  expect_identical(names(family$models), c("full", "IA", "additive", "SRR"))
  expect_null(family$models$SRR)
  expect_identical(family$models$IA$status, "run")
  expect_identical(family$models$additive$status, "failed")
  expect_identical(family$failures$model, c("additive", "SRR"))
  expect_identical(family$failures$stage, c("run", "create"))

  # The restricted models start from the estimates of the full model.
  expect_identical(family$models$full$spec$start, stats::setNames(numeric(), character()))
  expect_identical(family$models$IA$spec$start, c(b1 = 0.326, b2 = 0.166, b4 = 0.069))
  expect_match(family$models$additive$mplus$MODEL, "Z ON Y*0.166 (b2);", fixed = TRUE)

  expect_identical(
    basename(family$files$input),
    paste0("rsa_", c("full", "IA", "additive", "SRR"), ".inp")
  )
  expect_identical(basename(family$files$data), "rsa.dat")
  expect_setequal(
    list.files(output_dir),
    c("rsa.dat", "rsa_full.inp", "rsa_IA.inp", "rsa_additive.inp")
  )
  expect_match(
    paste(readLines(family$files$input[["IA"]]), collapse = "\n"),
    "rsa.dat\";",
    fixed = TRUE
  )

  comparison <- compare_rsa_mplus_models(family)
  expect_identical(comparison$model, c("full", "IA", "additive", "SRR"))
  expect_identical(comparison$flags[3:4], c("failed", "failed"))
  messages <- attr(comparison, "messages")
  expect_identical(messages$type[1:2], c("failure", "failure"))
  expect_match(messages$message[[2]], "^create: Parameter label `c` is reserved")

  # Existing files are not overwritten, and nothing is written.
  expect_error(
    fit_rsa_mplus_family(dat, measurement, models = "full", modelout = modelout),
    "Refusing to overwrite existing file"
  )
})

test_that("the family runner validates its arguments before writing", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(x1 = 1:5, x2 = 2:6, y1 = 3:7, y2 = 4:8, z1 = 5:9, z2 = 6:10)
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  output_dir <- withr::local_tempdir()
  modelout <- file.path(output_dir, "rsa.inp")
  family <- function(...) {
    fit_rsa_mplus_family(
      dat,
      measurement,
      modelout = modelout,
      Mplus_command = "mplus",
      ...
    )
  }

  expect_snapshot(error = TRUE, {
    family(models = c("full", "weak", "cubic"))
    family(models = c("full", "full"))
    family(models = "SRRR")
    family(start = "SRRR")
    fit_rsa_mplus_family(dat, measurement, modelout = "rsa")
  })
  expect_error(
    fit_rsa_mplus_family(
      dat,
      measurement,
      modelout = file.path(output_dir, "missing", "rsa.inp")
    ),
    "Output directory does not exist"
  )
  expect_error(
    fit_rsa_mplus_family(
      dat,
      list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = "w"),
      modelout = modelout,
      Mplus_command = "mplus"
    ),
    "No model of the family could be generated"
  )
  expect_length(list.files(output_dir), 0L)
})
