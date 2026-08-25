yao_ma_test_parameters <- function(
  cs = FALSE,
  cc = FALSE,
  intercept = FALSE,
  slope = FALSE,
  valence = "positive"
) {
  axis <- if (valence == "positive") c("P10", "P11") else c("P20", "P21")
  estimates <- c(
    if (cs) 0.4 else 0,
    if (cc) 0.4 else 0,
    0,
    if (valence == "positive") -0.5 else 0.5,
    if (intercept) 0.4 else 0,
    if (slope) 1.4 else 1
  )
  nulls <- c(0, 0, 0, 0, 0, 1)
  different <- c(cs, cc, FALSE, TRUE, intercept, slope)
  half_width <- ifelse(different, 0.1, 0.2)
  data.frame(
    parameter = c("CS", "CC", "IS", "IC", axis),
    estimate = estimates,
    conf.low = ifelse(different, estimates - half_width, nulls - half_width),
    conf.high = ifelse(different, estimates + half_width, nulls + half_width)
  )
}

test_that("all 16 Yao-Ma types map from the expected conditions", {
  combinations <- expand.grid(
    cs = c(FALSE, TRUE),
    cc = c(FALSE, TRUE),
    intercept = c(FALSE, TRUE),
    slope = c(FALSE, TRUE)
  )

  observed <- vapply(
    seq_len(nrow(combinations)),
    function(i) {
      args <- combinations[i, ]
      result <- classify_yao_ma(
        yao_ma_test_parameters(
          cs = args$cs,
          cc = args$cc,
          intercept = args$intercept,
          slope = args$slope
        ),
        valence = "positive"
      )
      result$summary$type
    },
    integer(1)
  )

  offset <- combinations$cs + 2L * combinations$cc
  base <- 1L + 4L * combinations$intercept + 8L * combinations$slope
  expect_setequal(observed, 1:16)
  expect_equal(observed, base + offset)
})

test_that("negative valence selects the second principal axis", {
  parameters <- yao_ma_test_parameters(
    intercept = TRUE,
    valence = "negative"
  )
  result <- classify_yao_ma(parameters, valence = "negative")

  expect_equal(result$summary$type, 5L)
  expect_identical(result$summary$direction, "y_greater_than_x")
  expect_equal(
    result$decisions$parameter[result$decisions$condition == "axis_slope"],
    "P21"
  )
})

test_that("equivalence margins allow indeterminate decisions", {
  parameters <- yao_ma_test_parameters(cs = FALSE)
  parameters[parameters$parameter == "CS", c("conf.low", "conf.high")] <-
    c(-0.2, 0.2)
  margins <- c(
    CS = 0.1,
    CC = 0.1,
    IS = 0.1,
    IC = 0.1,
    axis_intercept = 0.1,
    axis_slope = 0.1
  )
  result <- classify_yao_ma(
    parameters,
    valence = "positive",
    equivalence = margins
  )

  expect_identical(result$summary$status, "indeterminate")
  expect_identical(
    result$decisions$state[result$decisions$condition == "CS"],
    "indeterminate"
  )
})

test_that("equivalence margins can establish a complete type", {
  parameters <- yao_ma_test_parameters()
  nulls <- c(CS = 0, CC = 0, IS = 0, P10 = 0, P11 = 1)
  for (parameter in names(nulls)) {
    row <- parameters$parameter == parameter
    parameters[row, c("estimate", "conf.low", "conf.high")] <-
      c(
        nulls[[parameter]],
        nulls[[parameter]] - 0.05,
        nulls[[parameter]] + 0.05
      )
  }
  margins <- c(
    CS = 0.1,
    CC = 0.1,
    IS = 0.1,
    IC = 0.1,
    axis_intercept = 0.1,
    axis_slope = 0.1
  )
  result <- classify_yao_ma(
    parameters,
    valence = "positive",
    equivalence = margins
  )

  expect_equal(result$summary$type, 1L)
  expect_identical(result$summary$status, "classified")
  expect_setequal(
    result$decisions$state,
    c("different", "equivalent")
  )
})

test_that("a missing principal axis returns a partial result", {
  parameters <- yao_ma_test_parameters()
  parameters <- parameters[!grepl("^P", parameters$parameter), ]
  result <- classify_yao_ma(parameters, valence = "positive")

  expect_identical(result$summary$status, "partial")
  expect_identical(is.na(result$summary$type), TRUE)
  expect_match(result$summary$label, "principal-axis subtype unavailable")
})

test_that("opposite curvature lies outside the typology", {
  parameters <- yao_ma_test_parameters()
  parameters[
    parameters$parameter == "IC",
    c("estimate", "conf.low", "conf.high")
  ] <-
    c(0.5, 0.3, 0.7)
  result <- classify_yao_ma(parameters, valence = "positive")

  expect_identical(result$summary$status, "outside_typology")
  expect_match(result$summary$label, "opposite sign")
})

test_that("classification validates parameter tables and margins", {
  expect_snapshot(
    error = TRUE,
    classify_yao_ma(data.frame(parameter = "CS"), valence = "positive")
  )
  expect_snapshot(
    error = TRUE,
    classify_yao_ma(
      yao_ma_test_parameters(),
      valence = "positive",
      equivalence = c(CS = 0.1)
    )
  )
})

test_that("cached Mplus output is classified without refitting", {
  skip_if_not_installed("MplusAutomation")
  model_path <- system.file(
    "extdata",
    "mplus-response-surfaces.out",
    package = "rsahelpers"
  )
  if (!nzchar(model_path)) {
    model_path <- normalizePath(
      file.path("inst", "extdata", "mplus-response-surfaces.out"),
      mustWork = TRUE
    )
  }
  model <- MplusAutomation::readModels(model_path, quiet = TRUE)
  result <- classify_yao_ma(model, valence = "positive")
  mplus_object <- structure(list(results = model), class = "mplusObject")
  workflow <- new_rsa_mplus_workflow(results = model, status = "read")
  rsa_result <- RSA_mplus(
    model,
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    plot = FALSE
  )

  expect_s3_class(result, "yao_ma_classification")
  expect_equal(result$summary$type, 2L)
  expect_identical(result$summary$status, "classified")
  expect_identical(result$summary$direction, "none")
  expect_setequal(
    result$decisions$parameter,
    c("CS", "CC", "IS", "IC", "P10", "P11")
  )
  expect_equal(
    classify_yao_ma(mplus_object, valence = "positive")$summary$type,
    2L
  )
  expect_equal(
    classify_yao_ma(workflow, valence = "positive")$summary$type,
    2L
  )
  expect_equal(
    classify_yao_ma(rsa_result, valence = "positive")$summary$type,
    2L
  )
})
