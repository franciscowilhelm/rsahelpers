workshop_data <- function() {
  skip_if_not_installed("haven")

  path <- system.file(
    "extdata",
    "spline.dta",
    package = "rsahelpers",
    mustWork = TRUE
  )
  haven::read_dta(path)
}

test_that("formula order maps predictors to X and Y", {
  d <- workshop_data()

  prepared <- prepare_congruence_data(JOBSAT ~ ATHWC * ATHHC, d, center = FALSE)

  expect_equal(prepared$x, as.numeric(d$ATHWC))
  expect_equal(prepared$y, as.numeric(d$ATHHC))
  expect_identical(
    attr(prepared, "variables"),
    c(z = "JOBSAT", x = "ATHWC", y = "ATHHC")
  )
})

test_that("formula interface supports non-syntactic column names", {
  d <- data.frame(
    check.names = FALSE,
    "Z score" = c(1, 2, 3, NA),
    "X score" = c(2, 3, 4, 5),
    "Y score" = c(4, 3, 2, 1)
  )

  prepared <- prepare_congruence_data(
    `Z score` ~ `X score` * `Y score`,
    d,
    center = FALSE
  )

  expect_equal(nrow(prepared), 3)
  expect_equal(prepared$.row, 1:3)
  expect_equal(prepared$x, d[["X score"]][1:3])
  expect_equal(prepared$y, d[["Y score"]][1:3])
})

test_that("formula interface rejects unsupported specifications", {
  d <- data.frame(z = 1:4, x = 2:5, y = 5:2, group = letters[1:4])

  expect_snapshot(error = TRUE, prepare_congruence_data(z ~ x + y, d))
  expect_snapshot(error = TRUE, prepare_congruence_data(z ~ x * x, d))
  expect_snapshot(error = TRUE, prepare_congruence_data(z ~ x * missing, d))
  expect_snapshot(error = TRUE, prepare_congruence_data(z ~ x * group, d))
})

test_that("one-seam authority model matches workshop estimates", {
  d <- workshop_data()

  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  expect_equal(
    coef(fit),
    c(
      b0 = 4.2148667,
      b1 = 0.7556690,
      b2 = -0.4927952,
      b3 = 0.8819460,
      c0 = 0.5992022,
      c1 = 0.8302906
    ),
    tolerance = 1e-4
  )
  expect_equal(fit$rss, 677.3972, tolerance = 1e-3)
})

test_that("one-seam variety model matches workshop estimates", {
  d <- workshop_data()

  fit <- fit_spline_congruence(
    JOBSAT ~ VARWC * VARHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  expect_equal(
    coef(fit),
    c(
      b0 = 3.8154259,
      b1 = 0.4012231,
      b2 = -0.1942850,
      b3 = 0.8808640,
      c0 = -0.0770789,
      c1 = 0.6281473
    ),
    tolerance = 1e-4
  )
  expect_equal(fit$rss, 660.0655, tolerance = 1e-3)
})

test_that("two-seam authority model reaches the workshop local solution", {
  d <- workshop_data()

  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 2,
    center = FALSE
  )

  expect_lt(fit$rss, 676.2)
})

test_that("one-seam joint tests include expected terms", {
  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  tests <- suppressWarnings(spline_tests(fit))

  expect_setequal(
    tests$joint$term,
    c(
      "deviation_from_no_seam",
      "absolute_difference_constraints",
      "seam_equals_y_equals_x"
    )
  )
})

test_that("piecewise helper returns all comparison models", {
  d <- workshop_data()

  piecewise <- fit_piecewise_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 2,
    center = FALSE
  )
  piecewise_summary <- tidy_lm_summary(piecewise)

  expect_setequal(
    unique(piecewise_summary$model),
    c(
      "absolute_difference",
      "linear",
      "one_break",
      "constrained_one_seam",
      "fixed_two_seam"
    )
  )
})

test_that("bootstrap helper returns coefficient interval columns", {
  skip_if_not_installed("boot")

  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  boot_smoke <- bootstrap_spline(fit, R = 5)

  expect_named(boot_smoke, c("term", "estimate", "lower", "upper"))
})

test_that("surface plot helper returns plotted surface components", {
  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )
  png_file <- withr::local_tempfile(fileext = ".png")

  grDevices::png(png_file)
  plot_info <- plot_spline_surface(fit, grid_size = 12)
  grDevices::dev.off()

  expect_gt(file.size(png_file), 0)
  expect_named(
    plot_info,
    c("x", "y", "z", "transform", "congruence", "incongruence", "seams")
  )
})

test_that("pooled centering and scaling are the defaults", {
  d <- workshop_data()

  prepared <- prepare_congruence_data(JOBSAT ~ ATHWC * ATHHC, d)

  expect_equal(attr(prepared, "center_method"), "pooled")
  expect_equal(attr(prepared, "scale_method"), "pooled")
  expect_equal(mean(c(prepared$x, prepared$y)), 0, tolerance = 1e-8)
  expect_equal(sd(c(prepared$x, prepared$y)), 1, tolerance = 1e-6)
})

test_that("logical center keeps legacy behaviour (variablewise, no scaling)", {
  d <- workshop_data()

  vw <- prepare_congruence_data(JOBSAT ~ ATHWC * ATHHC, d, center = TRUE)
  expect_equal(attr(vw, "center_method"), "variablewise")
  expect_equal(attr(vw, "scale_method"), "none")
  expect_equal(mean(vw$x), 0, tolerance = 1e-8)
  expect_equal(mean(vw$y), 0, tolerance = 1e-8)

  none <- prepare_congruence_data(JOBSAT ~ ATHWC * ATHHC, d, center = FALSE)
  expect_equal(none$x, as.numeric(d$ATHWC))
})

test_that("hinge_offset parameterizes the fixed two-seam hinges", {
  d <- workshop_data()

  one <- prepare_congruence_data(JOBSAT ~ ATHWC * ATHHC, d, center = FALSE)
  two <- prepare_congruence_data(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    center = FALSE,
    hinge_offset = 2
  )

  expect_false(isTRUE(all.equal(one$hinge_upper, two$hinge_upper)))
  expected_upper <- (two$y - 2 - two$x) * as.numeric(two$y < 2 + two$x)
  expect_equal(two$hinge_upper, expected_upper)
})

test_that("bootstrap supports BCa intervals and reports failure rate", {
  skip_if_not_installed("boot")

  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  set.seed(42)
  boot_bca <- suppressWarnings(bootstrap_spline(fit, R = 40, type = "bca"))

  expect_named(boot_bca, c("term", "estimate", "lower", "upper"))
  expect_false(is.null(attr(boot_bca, "fail_rate")))
  expect_true(all(is.finite(boot_bca$lower)))
})

test_that("logLik, nobs, AIC, and summary methods work", {
  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  expect_equal(nobs(fit), nrow(fit$data))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_equal(attr(logLik(fit), "df"), length(coef(fit)) + 1L)

  smry <- summary(fit, warn_seam = FALSE)
  expect_setequal(
    names(smry$coefficients),
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
})

test_that("spline_tests warns about approximate seam standard errors", {
  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  expect_warning(spline_tests(fit), "seam")
  expect_silent(spline_tests(fit, warn_seam = FALSE))
})

test_that("compare_congruence_models returns nested comparison statistics", {
  d <- workshop_data()
  ols <- fit_piecewise_congruence(JOBSAT ~ ATHWC * ATHHC, d, center = FALSE)

  cmp <- compare_congruence_models(
    linear = ols$linear,
    piecewise = ols$one_break
  )

  expect_setequal(
    names(cmp),
    c(
      "model",
      "npar",
      "df.residual",
      "rss",
      "r.squared",
      "AIC",
      "df",
      "deltaR2",
      "F",
      "p.value"
    )
  )
  expect_true(is.finite(cmp$F[2]))
  expect_true(is.na(cmp$F[1]))
})

test_that("a larger model that fits worse triggers the simpler-model guard", {
  fake_fit <- function(npar, rss, n = 100) {
    structure(
      list(
        n_seams = 1L,
        coefficients = stats::setNames(
          rep(0, npar),
          paste0("p", seq_len(npar))
        ),
        fixed = NULL,
        data = data.frame(x = seq_len(n), y = seq_len(n), z = seq_len(n)),
        rss = rss,
        r.squared = 1 - rss / 1000,
        df.residual = n - npar
      ),
      class = "congruence_spline"
    )
  }

  simpler <- fake_fit(3, 100)
  larger <- fake_fit(5, 120) # more parameters but higher RSS

  expect_warning(
    compare_congruence_models(simpler = simpler, larger = larger),
    "simpler"
  )
})

test_that("fixing seam parameters yields a constrained nested fit", {
  d <- workshop_data()
  free <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )
  fixed <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE,
    fix = c(c0 = 0, c1 = 1)
  )

  expect_equal(unname(coef(fixed)[c("c0", "c1")]), c(0, 1))
  expect_equal(fixed$df.residual, free$df.residual + 2L)
  expect_gte(fixed$rss, free$rss)
})

test_that("two-seam surface features include section slopes and crossing", {
  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 2,
    center = FALSE
  )

  feats <- surface_features(fit)
  expect_true(all(
    c("seam1_x_slope", "both_y_slope", "seams_cross", "n_sections") %in%
      names(feats)
  ))
  expect_true(feats[["n_sections"]] %in% c(3, 4))
})

test_that("tidy_lm_summary includes inferential columns", {
  d <- workshop_data()
  ols <- fit_piecewise_congruence(JOBSAT ~ ATHWC * ATHHC, d, center = FALSE)

  tidied <- tidy_lm_summary(ols)
  expect_true(all(c("std.error", "statistic", "p.value") %in% names(tidied)))
})

test_that("plot_spline_contour returns a ggplot", {
  skip_if_not_installed("ggplot2")

  d <- workshop_data()
  fit <- fit_spline_congruence(
    JOBSAT ~ ATHWC * ATHHC,
    d,
    n_seams = 1,
    center = FALSE
  )

  expect_s3_class(plot_spline_contour(fit), "ggplot")
})
