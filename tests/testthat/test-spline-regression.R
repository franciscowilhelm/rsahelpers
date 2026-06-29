workshop_data <- function() {
  skip_if_not_installed("haven")

  path <- system.file("extdata", "spline.dta",
    package = "splinecongruence",
    mustWork = TRUE
  )
  haven::read_dta(path)
}

test_that("component variables follow field convention", {
  d <- workshop_data()

  prepared <- prepare_congruence_data(d, "ATHWC", "ATHHC", "JOBSAT",
    center = FALSE
  )

  expect_equal(prepared$x, as.numeric(d$ATHWC))
  expect_equal(prepared$y, as.numeric(d$ATHHC))
})

test_that("one-seam authority model matches workshop estimates", {
  d <- workshop_data()

  fit <- fit_spline_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
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

  fit <- fit_spline_congruence(d, "VARWC", "VARHC", "JOBSAT",
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

  fit <- fit_spline_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
    n_seams = 2,
    center = FALSE
  )

  expect_lt(fit$rss, 676.2)
})

test_that("one-seam joint tests include expected terms", {
  d <- workshop_data()
  fit <- fit_spline_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
    n_seams = 1,
    center = FALSE
  )

  tests <- spline_tests(fit)

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

  piecewise <- fit_piecewise_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
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
  fit <- fit_spline_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
    n_seams = 1,
    center = FALSE
  )

  boot_smoke <- bootstrap_spline(fit, R = 5)

  expect_named(boot_smoke, c("term", "estimate", "lower", "upper"))
})

test_that("surface plot helper returns plotted surface components", {
  d <- workshop_data()
  fit <- fit_spline_congruence(d, "ATHWC", "ATHHC", "JOBSAT",
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
