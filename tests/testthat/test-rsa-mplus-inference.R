# Parsed results with the pieces rsa_mplus_draws() reads: the regressions of
# Z, their TECH1 numbers, and a TECH3 covariance matrix (lower triangle).
mock_surface_results <- function(est, numbers, covariance) {
  predictors <- c("X", "Y", "XS", "XY", "YS")
  unstandardized <- data.frame(
    paramHeader = "Z.ON",
    param = predictors,
    est = est,
    se = 0.05,
    est_se = est / 0.05,
    pval = 0.01
  )
  names <- c("X", "Y", "Z", "XY", "XS", "YS")
  beta <- matrix(0L, 6, 6, dimnames = list(names, names))
  beta["Z", predictors] <- numbers
  covariance[upper.tri(covariance)] <- NA
  structure(
    list(
      parameters = list(unstandardized = unstandardized),
      tech1 = list(parameterSpecification = list(X = list(beta = beta))),
      tech3 = list(paramCov = covariance)
    ),
    class = c("mplus.model", "list")
  )
}

draw_surface <- function(results, ...) {
  rsa_mplus_draws(
    results,
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    ...
  )
}

test_that("restrictions and surface parameters are derived from b1 to b5", {
  b <- data.frame(b1 = 0.3, b2 = 0.2, b3 = -0.1, b4 = 0.25, b5 = -0.15)

  full <- rsa_mplus_derive_surface(b, rsa_mplus_model_spec("full"))
  expect_equal(full$a1, 0.5)
  expect_equal(full$a2, 0)
  expect_equal(full$a3, 0.1)
  expect_equal(full$a4, -0.5)
  expect_equal(full$a5, 0.05)
  # The stationary point solves the two first-order conditions.
  expect_equal(0.3 + 2 * -0.1 * full$x0 + 0.25 * full$y0, 0)
  expect_equal(0.2 + 2 * -0.15 * full$y0 + 0.25 * full$x0, 0)
  expect_equal(full$p10, full$y0 - full$p11 * full$x0)

  srr <- rsa_mplus_derive_surface(b, rsa_mplus_model_spec("SRR"))
  expect_equal(srr$b3, -0.15)
  expect_equal(srr$b4, 0.3)
  expect_equal(srr$a2, 0)
  expect_equal(srr$a5, 0)
  expect_equal(srr$c, (0.3 - 0.2) / (4 * -0.15))
  expect_equal(srr$bm, 0.5)

  additive <- rsa_mplus_derive_surface(b, rsa_mplus_model_spec("additive"))
  expect_equal(unlist(additive[c("b3", "b4", "b5", "a2", "a4")]), rep(0, 5), ignore_attr = TRUE)

  srrr <- rsa_mplus_derive_surface(b, rsa_mplus_model_spec("SRRR"))
  expect_equal(srrr$b3, 0.25^2 / (4 * -0.15))
  expect_equal(srrr$s, 0.25 / 0.3)
})

test_that("draws reproduce the TECH3 covariance of the free terms", {
  skip_if_not_installed("MplusAutomation")
  skip_if_not_installed("MASS")

  set.seed(8)
  root <- matrix(stats::rnorm(64, sd = 0.03), 8)
  covariance <- crossprod(root) + diag(0.001, 8)
  est <- c(0.3, 0.2, -0.1, 0.25, -0.15)
  results <- mock_surface_results(est, c(4L, 5L, 6L, 7L, 8L), covariance)

  full <- draw_surface(results, n = 50000, seed = 1)
  expect_s3_class(full, "rsa_mplus_draws")
  expect_identical(full$rsa_model, "full")
  expect_equal(colMeans(full$draws[1:5]), est, tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(
    stats::cov(full$draws[1:5]),
    covariance[4:8, 4:8],
    tolerance = 0.05,
    ignore_attr = TRUE
  )
  a1 <- full$summary[full$summary$parameter == "a1", ]
  expect_equal(a1$est, 0.5)
  expect_equal(
    a1$se,
    sqrt(covariance[4, 4] + covariance[5, 5] + 2 * covariance[5, 4]),
    tolerance = 0.02
  )
  expect_lt(a1$lower, a1$est)
  expect_gt(a1$upper, a1$est)
  expect_true(all(c("x0", "y0", "p11", "p21") %in% names(full$draws)))

  # The seed is local and reproducible.
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  again <- draw_surface(results, n = 100, seed = 1)
  expect_identical(stats::runif(1), before)
  expect_identical(again$draws, draw_surface(results, n = 100, seed = 1)$draws)
})

test_that("draws satisfy the restrictions of the model", {
  skip_if_not_installed("MplusAutomation")
  skip_if_not_installed("MASS")

  covariance <- diag(0.002, 8)
  # Printed estimates of an SRR model do not satisfy b4 = -2*b5 exactly.
  results <- mock_surface_results(
    c(0.284, 0.262, -0.107, 0.213, -0.107),
    c(4L, 5L, 6L, 7L, 8L),
    covariance
  )
  srr <- draw_surface(results, rsa_model = "SRR", n = 500, seed = 2)
  expect_equal(srr$draws$b3, srr$draws$b5)
  expect_equal(srr$draws$b4, -2 * srr$draws$b5)
  expect_equal(srr$draws$a2, numeric(500))
  expect_equal(srr$summary$est[srr$summary$parameter == "b4"], 0.214)
  expect_true("c" %in% srr$summary$parameter)

  # Terms that share a parameter share its draws; terms fixed at 0 stay 0.
  shared <- mock_surface_results(
    c(0.3, 0.3, 0, 0.2, 0),
    c(4L, 4L, 0L, 5L, 0L),
    covariance
  )
  ia <- draw_surface(shared, rsa_model = "IA", n = 200, seed = 3)
  expect_identical(ia$draws$b1, ia$draws$b2)
  expect_identical(ia$draws$b3, numeric(200))
  expect_gt(stats::sd(ia$draws$b4), 0)
})

test_that("draws need TECH1 and TECH3", {
  skip_if_not_installed("MplusAutomation")
  skip_if_not_installed("MASS")

  results <- mock_surface_results(
    c(0.3, 0.2, -0.1, 0.25, -0.15),
    4:8,
    diag(0.002, 8)
  )
  results$tech3 <- NULL
  expect_error(draw_surface(results), "no TECH1 or TECH3 section", fixed = TRUE)
  expect_error(
    draw_surface(results, probs = c(0.9, 0.1)),
    "`probs` must give a lower and an upper quantile",
    fixed = TRUE
  )
})

test_that("the supported region follows from the latent moments", {
  parameters <- data.frame(
    paramHeader = c("Y.WITH", "Means", "Means", "Variances", "Variances"),
    param = c("X", "X", "Y", "X", "Y"),
    est = c(0.4, 0.5, -0.5, 1.0, 1.5)
  )
  results <- structure(
    list(parameters = list(unstandardized = parameters)),
    class = c("mplus.model", "list")
  )
  range <- rsa_mplus_supported_range(results, pred_x = "X", pred_y = "Y")
  expect_equal(
    range$loc,
    stats::qnorm(c(0.05, 0.95), 0, sqrt((2.5 + 0.8) / 4))
  )
  expect_equal(
    range$loic,
    stats::qnorm(c(0.05, 0.95), 0.5, sqrt((2.5 - 0.8) / 4))
  )

  # Means that are not estimated are 0, and several models are averaged.
  other <- results
  other$parameters$unstandardized <- parameters[-(2:3), ]
  other$parameters$unstandardized$est <- c(0.2, 2.0, 0.5)
  averaged <- rsa_mplus_supported_range(
    list(results, other),
    pred_x = "X",
    pred_y = "Y",
    probs = c(0.1, 0.9)
  )
  expect_equal(
    averaged$moments,
    c(mx = 0.25, my = -0.25, vx = 1.5, vy = 1.0, cxy = 0.3)
  )
  expect_error(
    rsa_mplus_supported_range(results, pred_x = "X", pred_y = "Q"),
    "latent variances and covariance of `X` and `Q`",
    fixed = TRUE
  )
  expect_error(rsa_mplus_supported_range(results), "`pred_x` and `pred_y`")
})

test_that("the opt-in Mplus integration draws from a fitted ridge model", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")
  skip_if_not_installed("MASS")

  set.seed(9)
  n <- 400
  x <- stats::rnorm(n)
  y <- 0.4 * x + stats::rnorm(n, sd = 0.9)
  z <- 0.3 * x + 0.2 * y - 0.15 * (x - y)^2 + stats::rnorm(n, sd = 0.7)
  item <- function(f) f + stats::rnorm(n, sd = 0.5)
  dat <- data.frame(
    x1 = item(x), x2 = item(x), x3 = item(x),
    y1 = item(y), y2 = item(y), y3 = item(y),
    z1 = item(z), z2 = item(z), z3 = item(z)
  )
  directory <- withr::local_tempdir()
  fitted <- rsa_mplus_workflow(
    dat,
    measurement = list(
      X = c("x1", "x2", "x3"),
      Y = c("y1", "y2", "y3"),
      Z = c("z1", "z2", "z3")
    ),
    rsa_model = "SRR",
    blocks = list(OUTPUT = "TECH1 TECH3;"),
    modelout = file.path(directory, "srr.inp"),
    run = TRUE
  )
  intervals <- rsa_mplus_draws(fitted, n = 20000, seed = 1)
  new <- fitted$results$parameters$unstandardized
  new <- new[new$paramHeader == "New.Additional.Parameters", ]
  mplus <- function(label) new[new$param == label, ]
  ours <- function(label) intervals$summary[intervals$summary$parameter == label, ]
  # Monte Carlo and delta-method standard errors agree for linear functions.
  expect_equal(ours("a1")$est, mplus("CS")$est, tolerance = 0.01)
  expect_equal(ours("a1")$se, mplus("CS")$se, tolerance = 0.05)
  expect_equal(ours("a4")$se, mplus("IC")$se, tolerance = 0.05)
  expect_equal(ours("c")$est, mplus("C")$est, tolerance = 0.05)

  range <- rsa_mplus_supported_range(fitted)
  expect_lt(range$loc[[1]], 0)
  expect_gt(range$loc[[2]], range$loic[[2]])
})
