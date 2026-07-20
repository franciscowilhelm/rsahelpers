find_rsa_mplus_example <- function() {
  installed_path <- system.file(
    "extdata",
    "congruence_sim.out",
    package = "rsahelpers"
  )
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  normalizePath(
    file.path("inst", "extdata", "congruence_sim.out"),
    mustWork = TRUE
  )
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
