find_rsa_mplus_example <- function(file = "congruence_sim.out") {
  installed_path <- system.file("extdata", file, package = "rsahelpers")
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  normalizePath(file.path("inst", "extdata", file), mustWork = TRUE)
}

find_rsa_mplus_bayes_example <- function() {
  find_rsa_mplus_example("congruence_bayes.out")
}

test_that("RSA_mplus extracts unstandardized coefficients from Mplus output", {
  skip_if_not_installed("MplusAutomation")

  model_path <- find_rsa_mplus_example()

  expect_snapshot(
    result <- RSA_mplus(
      model = model_path,
      outcome = "Z",
      pred_x = "X",
      pred_y = "Y",
      pred_x2 = "XS",
      pred_xy = "XY",
      pred_y2 = "YS",
      new_labels = c("CS", "CC", "IS", "IC", "A5"),
      plot = FALSE
    )
  )

  expect_s3_class(result, "rsa_mplus")
  expect_equal(
    unname(result$coefficients[c("x", "y", "x2", "xy", "y2", "b0")]),
    c(0.326, 0.166, -0.058, 0.069, -0.085, 0),
    tolerance = 1e-8
  )
  expect_equal(result$regression_parameters$term, c("x", "y", "x2", "xy", "y2"))
  expect_equal(result$new_parameters$Label, c("CS", "CC", "IS", "IC", "A5"))
  expect_null(result$plot)

  # Frequentist output carries no credibility interval, so the annotation
  # falls back to p-values (see the Bayesian test for the other branch).
  expect_disjoint(
    names(result$new_parameters),
    c("lower_2.5ci", "upper_2.5ci")
  )
})

test_that("RSA_mplus accepts an mplus.model object and standardized coefficients", {
  skip_if_not_installed("MplusAutomation")

  model_path <- find_rsa_mplus_example()
  mplus_model <- MplusAutomation::readModels(model_path, quiet = TRUE)

  expect_snapshot(
    result <- RSA_mplus(
      model = mplus_model,
      outcome = "Z",
      pred_x = "X",
      pred_y = "Y",
      pred_x2 = "XS",
      pred_xy = "XY",
      pred_y2 = "YS",
      coef_type = "stdyx",
      include_new = FALSE,
      plot = FALSE
    )
  )

  expect_equal(
    unname(result$coefficients[c("x", "y", "x2", "xy", "y2")]),
    c(0.432, 0.200, -0.080, 0.086, -0.097),
    tolerance = 1e-8
  )
  expect_null(result$new_parameters)
})

test_that("RSA_mplus annotates Mplus surface parameters with APA stars", {
  skip_if_not_installed("MplusAutomation")
  skip_if(!suppressWarnings(requireNamespace("RSA", quietly = TRUE)))

  result <- suppressWarnings(RSA_mplus(
    model = find_rsa_mplus_example(),
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    new_labels = c("CS", "CC", "IS", "IC", "A5"),
    plot = TRUE
  ))

  expect_s3_class(result$plot, "trellis")
  expect_identical(
    result$plot$panel.args.common$SPs,
    "a1: 0.49***    a2: -0.07*    a3: 0.16*    a4: -0.21*    a5: 0.03"
  )

  no_parameters <- suppressWarnings(RSA_mplus(
    model = find_rsa_mplus_example(),
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    include_new = FALSE,
    plot = TRUE,
    param = FALSE
  ))
  expect_identical(
    no_parameters$plot$panel.args.common$SPs,
    "a1: 0.49    a2: -0.07    a3: 0.16    a4: -0.21    a5: 0.03"
  )
})

test_that("surface annotation applies all thresholds and requires p-values", {
  parameters <- data.frame(
    Label = c("CS", "CC", "IS", "IC", "A5"),
    est = c(1, 2, 3, 4, 5),
    pval = c(0.05, 0.01, 0.001, 0.051, 0.50)
  )

  expect_identical(
    rsa_mplus_parameter_annotation(parameters),
    "a1: 1.00*    a2: 2.00**    a3: 3.00***    a4: 4.00    a5: 5.00"
  )

  parameters$pval[[1]] <- NA_real_
  expect_null(rsa_mplus_parameter_annotation(parameters))
  expect_null(rsa_mplus_parameter_annotation(parameters[-1, ]))
})

test_that("Bayesian output is starred from credibility intervals, not p-values", {
  skip_if_not_installed("MplusAutomation")

  result <- RSA_mplus(
    model = find_rsa_mplus_bayes_example(),
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    new_labels = c("CS", "CC", "IS", "IC", "A5"),
    plot = FALSE
  )

  expect_equal(
    result$new_parameters$lower_2.5ci,
    c(0.441, -0.124, -0.141, -0.598, -0.164),
    tolerance = 1e-8
  )
  expect_equal(
    result$new_parameters$upper_2.5ci,
    c(0.630, 0.040, 0.153, -0.221, 0.012),
    tolerance = 1e-8
  )

  # A5 has a one-tailed posterior p-value of .047 but a credibility interval
  # covering zero, so the p-value rule and the interval rule disagree.
  expect_equal(result$new_parameters$pval[[5]], 0.047, tolerance = 1e-8)
  expect_identical(
    rsa_mplus_parameter_annotation(result$new_parameters),
    "a1: 0.54*    a2: -0.04    a3: 0.01    a4: -0.41*    a5: -0.08"
  )
})

test_that("surface annotation prefers credibility intervals over p-values", {
  parameters <- data.frame(
    Label = c("CS", "CC", "IS", "IC", "A5"),
    est = c(1, 2, 3, 4, 5),
    pval = c(0, 0, 0, 0, 0),
    lower_2.5ci = c(0.5, -1, 2, -4, -0.1),
    upper_2.5ci = c(1.5, 1, 4, -3, 0.1)
  )

  expect_identical(
    rsa_mplus_parameter_annotation(parameters),
    "a1: 1.00*    a2: 2.00    a3: 3.00*    a4: 4.00*    a5: 5.00"
  )

  # Incomplete intervals fall back to the p-value thresholds.
  parameters$lower_2.5ci[[1]] <- NA_real_
  expect_identical(
    rsa_mplus_parameter_annotation(parameters),
    "a1: 1.00***    a2: 2.00***    a3: 3.00***    a4: 4.00***    a5: 5.00***"
  )
})
