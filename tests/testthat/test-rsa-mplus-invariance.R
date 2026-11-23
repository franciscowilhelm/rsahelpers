test_that("invariance levels share loadings and then intercepts", {
  design <- rsa_mplus_invariance_design(
    list(X = c("x1", "x2", "x3"), Y = c("y1", "y2", "y3")),
    pairs = NULL,
    free_intercepts = "y2"
  )
  expect_identical(design$pairs, list(c("X", "Y")))
  expect_identical(design$free, list(2L))

  configural <- rsa_mplus_invariance_syntax(design, "configural")
  expect_identical(
    configural,
    c(
      "X BY x1;", "X BY x2;", "X BY x3;",
      "Y BY y1;", "Y BY y2;", "Y BY y3;"
    )
  )

  metric <- rsa_mplus_invariance_syntax(design, "metric")
  expect_identical(metric[c(2, 5)], c("X BY x2 (l1_2);", "Y BY y2 (l1_2);"))
  expect_identical(metric[c(1, 4)], c("X BY x1;", "Y BY y1;"))
  expect_length(metric, 6L)

  partial <- rsa_mplus_invariance_syntax(design, "partial")
  expect_identical(partial[1:6], metric)
  expect_identical(
    partial[7:14],
    c(
      "[x1] (i1_1);", "[x2];", "[x3] (i1_3);",
      "[y1] (i1_1);", "[y2];", "[y3] (i1_3);",
      "[X@0];", "[Y];"
    )
  )

  scalar <- rsa_mplus_invariance_syntax(design, "scalar")
  expect_identical(scalar[c(8, 11)], c("[x2] (i1_2);", "[y2] (i1_2);"))
})

test_that("several pairs get their own labels in one model", {
  design <- rsa_mplus_invariance_design(
    list(
      A = c("a1", "a2"),
      B = c("b1", "b2"),
      C = c("c1", "c2"),
      D = c("d1", "d2")
    ),
    pairs = list(c("A", "C"), c("B", "D")),
    free_intercepts = "d2"
  )
  expect_identical(design$free, list(integer(), 2L))
  partial <- rsa_mplus_invariance_syntax(design, "partial")
  expect_true(all(
    c("A BY a2 (l1_2);", "C BY c2 (l1_2);", "B BY b2 (l2_2);", "D BY d2 (l2_2);") %in%
      partial
  ))
  expect_true(all(c("[a2] (i1_2);", "[c2] (i1_2);", "[b2];", "[d2];") %in% partial))
  expect_true(all(c("[A@0];", "[C];", "[B@0];", "[D];") %in% partial))
})

test_that("the invariance design is validated", {
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  expect_error(
    rsa_mplus_invariance_design(measurement, NULL, NULL),
    "`pairs` is needed",
    fixed = TRUE
  )
  expect_error(
    rsa_mplus_invariance_design(measurement, list(c("X", "Y")), NULL),
    "exactly one of the `pairs`",
    fixed = TRUE
  )
  expect_error(
    rsa_mplus_invariance_design(
      list(X = c("x1", "x2", "x3"), Y = c("y1", "y2")),
      NULL,
      NULL
    ),
    "same number of indicators",
    fixed = TRUE
  )
  expect_error(
    rsa_mplus_invariance_design(measurement[1:2], NULL, "q1"),
    "not found: `q1`",
    fixed = TRUE
  )
})

test_that("the scaled chi-square difference follows Satorra and Bentler", {
  # Scaled statistics 87.906 (df 12, c 1.20) and 65.052 (df 10, c 1.25).
  test <- rsa_mplus_scaled_chisq_difference(87.906, 12, 1.20, 65.052, 10, 1.25)
  cd <- (12 * 1.20 - 10 * 1.25) / 2
  expect_equal(test$df, 2)
  expect_equal(test$chisq, (87.906 * 1.20 - 65.052 * 1.25) / cd)
  expect_equal(test$p, stats::pchisq(test$chisq, 2, lower.tail = FALSE))

  # Without scaling factors it is the ordinary difference test.
  ml <- rsa_mplus_scaled_chisq_difference(30, 12, NA, 20, 10, NA)
  expect_equal(ml$chisq, 10)
  # A scaling correction that is not positive leaves the test undefined.
  undefined <- rsa_mplus_scaled_chisq_difference(30, 12, 1, 20, 10, 1.3)
  expect_true(is.na(undefined$chisq) && is.na(undefined$p))
})

test_that("levels are compared with the less restricted model", {
  level <- function(k, chisq, df, cfi) {
    list(results = list(summaries = data.frame(
      Parameters = k,
      ChiSqM_Value = chisq,
      ChiSqM_DF = df,
      ChiSqM_ScalingCorrection = 1,
      CFI = cfi,
      RMSEA_Estimate = 0.05,
      SRMR = 0.03
    )))
  }
  comparison <- rsa_mplus_invariance_comparison(list(
    configural = level(19, 58, 8, 0.970),
    metric = level(17, 65, 10, 0.967),
    partial = level(16, 67, 11, 0.966),
    scalar = level(15, 88, 12, 0.954)
  ))
  expect_identical(
    comparison$against,
    c(NA, "configural", "metric", "metric")
  )
  expect_equal(comparison$diff_chisq, c(NA, 7, 2, 23))
  expect_equal(comparison$diff_df, c(NA, 2, 1, 2))
  expect_equal(comparison$delta_CFI, c(NA, -0.003, -0.001, -0.013))
})

test_that("writing the invariance models does not need Mplus", {
  skip_if_not_installed("MplusAutomation")

  set.seed(6)
  dat <- as.data.frame(matrix(stats::rnorm(400), ncol = 4))
  names(dat) <- c("x1", "x2", "y1", "y2")
  directory <- withr::local_tempdir()
  written <- fit_rsa_mplus_invariance(
    dat,
    list(X = c("x1", "x2"), Y = c("y1", "y2")),
    modelout = file.path(directory, "inv.inp"),
    run = FALSE
  )
  expect_s3_class(written, "rsa_mplus_invariance")
  expect_named(written$models, c("configural", "metric", "scalar"))
  expect_null(written$comparison)
  expect_true(all(file.exists(written$files$input)))
  expect_true(file.exists(written$files$data))
  input <- readLines(written$files$input[["scalar"]])
  expect_true(any(grepl("ESTIMATOR = MLR;", input, fixed = TRUE)))
  expect_true(any(grepl("[x2] (i1_2);", input, fixed = TRUE)))
  expect_error(
    fit_rsa_mplus_invariance(
      dat,
      list(X = c("x1", "x2"), Y = c("y1", "y2")),
      modelout = file.path(directory, "inv.inp"),
      run = FALSE
    ),
    "Refusing to overwrite"
  )
})

test_that("the opt-in Mplus integration tests invariance", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("MplusAutomation")

  set.seed(7)
  n <- 500
  x <- stats::rnorm(n)
  y <- 0.4 * x + stats::rnorm(n) + 0.5
  item <- function(f, loading, intercept) {
    intercept + loading * f + stats::rnorm(n, sd = 0.6)
  }
  dat <- data.frame(
    x1 = item(x, 1, 0), x2 = item(x, 0.8, 0.2), x3 = item(x, 1.1, -0.1),
    y1 = item(y, 1, 0), y2 = item(y, 0.8, 0.9), y3 = item(y, 1.1, -0.1)
  )
  directory <- withr::local_tempdir()
  invariance <- fit_rsa_mplus_invariance(
    dat,
    list(X = c("x1", "x2", "x3"), Y = c("y1", "y2", "y3")),
    free_intercepts = "x2",
    modelout = file.path(directory, "inv.inp")
  )
  comparison <- invariance$comparison
  expect_identical(
    comparison$level,
    c("configural", "metric", "partial", "scalar")
  )
  expect_equal(comparison$df, c(8, 10, 11, 12))
  # The second pair has unequal intercepts: scalar fails, partial holds.
  expect_lt(comparison$diff_p[comparison$level == "scalar"], 0.001)
  expect_gt(comparison$diff_p[comparison$level == "partial"], 0.01)
  expect_gt(comparison$diff_p[comparison$level == "metric"], 0.01)
})
