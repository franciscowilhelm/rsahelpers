test_that("latent RSA-Mplus models combine tidySEM and XWITH syntax", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(1)
  dat <- as.data.frame(matrix(stats::rnorm(900), ncol = 9))
  names(dat) <- c("x1", "x2", "x3", "y1", "y2", "y3", "z1", "z2", "z3")

  model <- create_rsa_mplus_model(
    data = dat,
    measurement = list(
      X = c("x1", "x2", "x3"),
      Y = c("y1", "y2", "y3"),
      Z = c("z1", "z2", "z3")
    ),
    structural = list(Z = c("X", "Y"))
  )

  expect_s3_class(model, "rsa_mplus_workflow")
  expect_s3_class(model$mplus, "mplusObject")
  expect_identical(model$status, "created")
  expect_identical(model$spec$metadata$pred_xy, "XY")
  expect_snapshot(cat(
    model$mplus$MODEL,
    "\nMODEL CONSTRAINT:\n",
    model$mplus$MODELCONSTRAINT,
    "\n"
  ))
})

test_that("character syntax, covariates, and block customization are retained", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(2)
  dat <- as.data.frame(matrix(stats::rnorm(700), ncol = 7))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2", "age")

  model <- create_rsa_mplus_model(
    data = dat,
    measurement = c(
      "X =~ x1 + x2",
      "Y =~ y1 + y2",
      "Z =~ z1 + z2"
    ),
    structural = "Z ~ X + Y + age",
    constraints = "surface",
    blocks = list(
      analysis = "TYPE = RANDOM;\nALGORITHM = INTEGRATION;\nPROCESSORS = 2;",
      output = "STDYX CINTERVAL;"
    ),
    model_extra = "Z WITH age;",
    constraint_extra = "NEW(total); total = cs + cc;"
  )

  expect_match(model$mplus$MODEL, "Z ON AGE", ignore.case = TRUE)
  expect_match(model$mplus$MODEL, "Z WITH age;", fixed = TRUE)
  expect_match(model$mplus$ANALYSIS, "PROCESSORS = 2", fixed = TRUE)
  expect_match(model$mplus$MODELCONSTRAINT, "NEW(total)", fixed = TRUE)
  expect_no_match(model$mplus$MODELCONSTRAINT, "p11", fixed = TRUE)
})

test_that("SI-LMS fixes residual variances from supplied reliabilities", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(
    xmean = seq_len(10),
    ymean = seq_len(10) * 2,
    zmean = seq_len(10) * 3
  )
  reliability <- c(X = 0.8, Y = 0.9, Z = 0.7)

  model <- create_rsa_mplus_model(
    data = dat,
    measurement = list(X = "xmean", Y = "ymean", Z = "zmean"),
    model_type = "si_lms",
    reliability = reliability,
    constraints = "surface"
  )

  syntax_table <- tidySEM::syntax(model$tidysem)
  for (factor_name in c("X", "Y", "Z")) {
    indicator <- syntax_table$rhs[
      syntax_table$op == "=~" & syntax_table$lhs == factor_name
    ]
    residual <- syntax_table[
      syntax_table$op == "~~" &
        syntax_table$lhs == indicator &
        syntax_table$rhs == indicator,
      ,
      drop = FALSE
    ]
    expect_equal(residual$free, 0)
    expect_equal(
      residual$ustart,
      (1 - reliability[[factor_name]]) * stats::var(dat[[indicator]])
    )
  }
})

test_that("RSA-Mplus validation reports actionable model errors", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(x = 1:5, y = 2:6, z = 3:7)

  expect_snapshot(
    error = TRUE,
    create_rsa_mplus_model(
      dat,
      list(X = "x", Y = "y", Z = "z")
    )
  )
})

test_that("RSA-Mplus validation catches Mplus name collisions", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(
    longname_a = 1:5,
    longname_b = 2:6,
    y1 = 3:7,
    y2 = 4:8,
    z1 = 5:9,
    z2 = 6:10
  )

  expect_snapshot(
    error = TRUE,
    create_rsa_mplus_model(
      dat,
      list(
        X = c("longname_a", "longname_b"),
        Y = c("y1", "y2"),
        Z = c("z1", "z2")
      )
    )
  )
})

test_that("writing creates reproducible input and data files safely", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(3)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  model <- create_rsa_mplus_model(
    dat,
    list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  )
  output_dir <- withr::local_tempdir()
  input_path <- file.path(output_dir, "rsa.inp")

  written <- write_rsa_mplus_model(model, modelout = input_path)

  expect_identical(written$status, "written")
  expect_equal(file.exists(written$files$input), TRUE)
  expect_equal(file.exists(written$files$data), TRUE)
  expect_match(
    paste(readLines(input_path, warn = FALSE), collapse = "\n"),
    "XWITH"
  )
  expect_snapshot(
    error = TRUE,
    transform = function(x) {
      gsub("`[^`]+/rsa\\.(inp|dat)`", "`<file>`", x)
    },
    write_rsa_mplus_model(model, modelout = input_path)
  )
  expect_no_error(write_rsa_mplus_model(
    model,
    modelout = input_path,
    overwrite = TRUE
  ))
})

test_that("read workflows retain metadata for RSA_mplus inference", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  model_path <- system.file(
    "extdata",
    "congruence_sim.out",
    package = "rsahelpers"
  )
  if (!nzchar(model_path)) {
    model_path <- normalizePath(
      file.path("inst", "extdata", "congruence_sim.out"),
      mustWork = TRUE
    )
  }
  parsed <- MplusAutomation::readModels(model_path, quiet = TRUE)

  set.seed(4)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  workflow <- create_rsa_mplus_model(
    dat,
    list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  )
  workflow$results <- parsed
  workflow$status <- "read"

  result <- RSA_mplus(workflow, plot = FALSE)
  expect_equal(
    unname(result$coefficients[c("x", "y", "x2", "xy", "y2", "b0")]),
    c(0.326, 0.166, -0.058, 0.069, -0.085, 0),
    tolerance = 1e-8
  )
  expect_identical(result$workflow, workflow)

  object <- MplusAutomation::mplusObject(MODEL = "Z ON X;", quiet = TRUE)
  object$results <- parsed
  explicit <- RSA_mplus(
    object,
    outcome = "Z",
    pred_x = "X",
    pred_y = "Y",
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    b0 = 0,
    plot = FALSE
  )
  expect_equal(explicit$coefficients, result$coefficients)
})

test_that("the opt-in Mplus integration runs latent and SI-LMS models", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(2026)
  n <- 300
  x <- stats::rnorm(n)
  y <- 0.3 * x + stats::rnorm(n)
  z <- 0.3 *
    x +
    0.2 * y -
    0.1 * x^2 +
    0.15 * x * y -
    0.08 * y^2 +
    stats::rnorm(n, sd = 0.7)
  latent_data <- data.frame(
    x1 = x + stats::rnorm(n, sd = 0.4),
    x2 = x + stats::rnorm(n, sd = 0.4),
    x3 = x + stats::rnorm(n, sd = 0.4),
    y1 = y + stats::rnorm(n, sd = 0.4),
    y2 = y + stats::rnorm(n, sd = 0.4),
    y3 = y + stats::rnorm(n, sd = 0.4),
    z1 = z + stats::rnorm(n, sd = 0.4),
    z2 = z + stats::rnorm(n, sd = 0.4),
    z3 = z + stats::rnorm(n, sd = 0.4)
  )
  output_dir <- withr::local_tempdir()

  latent <- rsa_mplus_workflow(
    data = latent_data,
    measurement = list(
      X = c("x1", "x2", "x3"),
      Y = c("y1", "y2", "y3"),
      Z = c("z1", "z2", "z3")
    ),
    modelout = file.path(output_dir, "latent.inp"),
    run = TRUE,
    Mplus_command = "/Applications/Mplus/mplus"
  )
  expect_identical(latent$status, "run")
  expect_s3_class(RSA_mplus(latent, plot = FALSE), "rsa_mplus")

  si_lms <- rsa_mplus_workflow(
    data = data.frame(xmean = x, ymean = y, zmean = z),
    measurement = list(X = "xmean", Y = "ymean", Z = "zmean"),
    model_type = "si_lms",
    reliability = c(X = 0.8, Y = 0.8, Z = 0.8),
    modelout = file.path(output_dir, "si_lms.inp"),
    run = TRUE,
    Mplus_command = "/Applications/Mplus/mplus"
  )
  expect_identical(si_lms$status, "run")
  expect_s3_class(RSA_mplus(si_lms, plot = FALSE), "rsa_mplus")
})
