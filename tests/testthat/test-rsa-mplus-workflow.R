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
  expect_match(model$mplus$VARIABLE, "USEVARIABLES", fixed = TRUE)
  expect_match(model$mplus$VARIABLE, "age", fixed = TRUE)
})

test_that("USEVARIABLES is explicit and custom VARIABLE directives are retained", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(
    x1 = 1:5,
    x2 = 2:6,
    y1 = 3:7,
    y2 = 4:8,
    z1 = 5:9,
    z2 = 6:10,
    age = 20:24
  )
  analysis_variables <- setdiff(names(dat), "age")

  model <- create_rsa_mplus_model(
    data = dat,
    measurement = list(
      X = c("x1", "x2"),
      Y = c("y1", "y2"),
      Z = c("z1", "z2")
    ),
    usevariables = analysis_variables,
    blocks = list(VARIABLE = "IDVARIABLE = x1;")
  )

  expect_match(model$mplus$VARIABLE, "USEVARIABLES =", fixed = TRUE)
  expect_match(model$mplus$VARIABLE, "x1 x2 y1 y2 z1 z2", fixed = TRUE)
  expect_no_match(model$mplus$VARIABLE, "age", fixed = TRUE)
  expect_match(model$mplus$VARIABLE, "IDVARIABLE = x1;", fixed = TRUE)
  expect_identical(model$spec$usevariables, analysis_variables)

  expect_snapshot(
    error = TRUE,
    create_rsa_mplus_model(
      data = dat,
      measurement = list(
        X = c("x1", "x2"),
        Y = c("y1", "y2"),
        Z = c("z1", "z2")
      ),
      blocks = list(VARIABLE = "USEVARIABLES = x1 x2;")
    )
  )
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

test_that("effects coding frees the marker loadings and constrains their mean", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(3)
  dat <- as.data.frame(matrix(stats::rnorm(900), ncol = 9))
  names(dat) <- c("x1", "x2", "x3", "y1", "y2", "y3", "z1", "z2", "z3")
  measurement <- list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  )

  model <- create_rsa_mplus_model(
    dat,
    measurement,
    identification = "effects"
  )

  for (indicator in c("x1", "y1", "z1")) {
    expect_match(model$mplus$MODEL, paste0("BY ", indicator, "* ("), fixed = TRUE)
  }
  expect_no_match(model$mplus$MODEL, "@1", fixed = TRUE)
  expect_match(model$mplus$MODEL, "[X@0];", fixed = TRUE)
  expect_match(model$mplus$MODELCONSTRAINT, "lx1 = 3 - lx2 - lx3;", fixed = TRUE)
  expect_match(model$mplus$MODELCONSTRAINT, "lz1 = 3 - lz2 - lz3;", fixed = TRUE)
  expect_no_match(model$mplus$MODELCONSTRAINT, "ix1", fixed = TRUE)
  expect_null(model$mplus$DEFINE)
  expect_null(model$spec$centering)
})

test_that("pooled centering shifts X and Y by one constant and frees their means", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(4)
  dat <- as.data.frame(matrix(stats::rnorm(900), ncol = 9))
  names(dat) <- c("x1", "x2", "x3", "y1", "y2", "y3", "z1", "z2", "z3")
  dat[c("x1", "x2", "x3")] <- dat[c("x1", "x2", "x3")] + 4
  dat[c("y1", "y2", "y3")] <- dat[c("y1", "y2", "y3")] + 3

  model <- create_rsa_mplus_model(
    dat,
    list(
      X = c("x1", "x2", "x3"),
      Y = c("y1", "y2", "y3"),
      Z = c("z1", "z2", "z3")
    ),
    identification = "effects",
    center = "pooled",
    blocks = list(DEFINE = "w = x1 + y1;")
  )

  scale_means <- c(
    x = mean(colMeans(dat[c("x1", "x2", "x3")])),
    y = mean(colMeans(dat[c("y1", "y2", "y3")]))
  )
  expect_equal(model$spec$centering$scale_means, scale_means)
  expect_equal(model$spec$centering$constant, mean(scale_means))

  define <- strsplit(model$mplus$DEFINE, "\n", fixed = TRUE)[[1]]
  expect_length(define, 7L)
  expect_match(define[1], "^x1 = x1 - 3\\.")
  expect_match(define[6], "^y3 = y3 - 3\\.")
  expect_identical(define[7], "w = x1 + y1;")
  expect_no_match(model$mplus$DEFINE, "z1", fixed = TRUE)

  expect_match(model$mplus$MODEL, "[X];", fixed = TRUE)
  expect_match(model$mplus$MODEL, "[Y];", fixed = TRUE)
  expect_match(model$mplus$MODEL, "[Z@0];", fixed = TRUE)
  expect_match(model$mplus$MODELCONSTRAINT, "ix1 = -ix2 - ix3;", fixed = TRUE)
  expect_match(model$mplus$MODELCONSTRAINT, "iy1 = -iy2 - iy3;", fixed = TRUE)
  expect_no_match(model$mplus$MODELCONSTRAINT, "iz1", fixed = TRUE)
})

test_that("identification and centering options are validated", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  dat <- data.frame(
    x1 = 1:5,
    x2 = 2:6,
    y1 = 3:7,
    y2 = 4:8,
    z1 = 5:9,
    z2 = 6:10
  )
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))

  expect_error(
    create_rsa_mplus_model(dat, measurement, center = "pooled"),
    "requires `identification = \"effects\"`",
    fixed = TRUE
  )
  expect_error(
    create_rsa_mplus_model(
      data.frame(xm = 1:5, ym = 2:6, zm = 3:7),
      list(X = "xm", Y = "ym", Z = "zm"),
      model_type = "si_lms",
      reliability = c(X = 0.8, Y = 0.8, Z = 0.8),
      identification = "effects"
    ),
    "requires `model_type = \"latent\"`",
    fixed = TRUE
  )
  expect_error(
    create_rsa_mplus_model(
      dat,
      c("X =~ x1 + a*x2", "Y =~ y1 + y2", "Z =~ z1 + z2"),
      identification = "effects"
    ),
    "default (unlabelled, marker-identified) loadings for `X`",
    fixed = TRUE
  )
})

test_that("invariance identification shares loadings and intercepts of X and Y", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(4)
  dat <- as.data.frame(matrix(stats::rnorm(900, mean = 3), ncol = 9))
  names(dat) <- c("x1", "x2", "x3", "y1", "y2", "y3", "z1", "z2", "z3")
  measurement <- list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  )

  model <- create_rsa_mplus_model(
    dat,
    measurement,
    identification = "invariance",
    center = "pooled",
    free_intercepts = "y2"
  )
  syntax <- model$mplus$MODEL
  for (statement in c(
    "X BY x1@1;", "X BY x2 (lxy2);", "Y BY y2 (lxy2);",
    "X BY x3 (lxy3);", "Y BY y3 (lxy3);",
    "[x1@0];", "[y1@0];", "[x2];", "[y2];",
    "[x3](ixy3);", "[y3](ixy3);", "[X];", "[Y];", "[Z@0];", "Z BY z2;"
  )) {
    expect_match(syntax, statement, fixed = TRUE)
  }
  # The measurement model needs no MODEL CONSTRAINT, so it runs under Bayes.
  expect_no_match(model$mplus$MODELCONSTRAINT, "lxy|ixy")
  expect_identical(model$spec$invariance$free_intercepts, c("x2", "y2"))
  expect_identical(model$spec$invariance$indicators$y, c("y1", "y2", "y3"))

  constant <- mean(c(
    mean(colMeans(dat[c("x1", "x2", "x3")])),
    mean(colMeans(dat[c("y1", "y2", "y3")]))
  ))
  expect_equal(model$spec$centering$constant, constant)
  expect_match(model$mplus$DEFINE, "y3 = y3 - ", fixed = TRUE)

  # Without pooled centering the indicators are not shifted.
  uncentred <- create_rsa_mplus_model(
    dat,
    measurement,
    identification = "invariance"
  )
  expect_null(uncentred$mplus$DEFINE)
  expect_null(uncentred$spec$centering)
  expect_match(uncentred$mplus$MODEL, "[x2](ixy2);", fixed = TRUE)

  expect_error(
    create_rsa_mplus_model(
      dat,
      measurement,
      identification = "invariance",
      free_intercepts = "x1"
    ),
    "marker indicators",
    fixed = TRUE
  )
  expect_error(
    create_rsa_mplus_model(
      dat,
      measurement,
      identification = "invariance",
      free_intercepts = "z2"
    ),
    "must name X or Y indicators; not found: `z2`",
    fixed = TRUE
  )
  expect_error(
    create_rsa_mplus_model(dat, measurement, free_intercepts = "x2"),
    "`free_intercepts` requires `identification = \"invariance\"`",
    fixed = TRUE
  )
  expect_error(
    create_rsa_mplus_model(
      dat,
      list(X = c("x1", "x2", "x3"), Y = c("y1", "y2"), Z = c("z1", "z2", "z3")),
      identification = "invariance"
    ),
    "as many `X` as `Y` indicators",
    fixed = TRUE
  )
})

test_that("the models added to the RSA family write their restrictions", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(5)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  generate <- function(rsa_model) {
    create_rsa_mplus_model(dat, measurement, rsa_model = rsa_model)$mplus
  }

  a5zero <- generate("a5zero")
  expect_match(a5zero$MODELCONSTRAINT, "b3 = b5;", fixed = TRUE)
  expect_match(a5zero$MODELCONSTRAINT, "NEW(cs cc is ic);", fixed = TRUE)

  no_ia <- generate("noIA")
  expect_match(no_ia$MODEL, "Z ON XY@0 (b4);", fixed = TRUE)
  expect_match(no_ia$MODELCONSTRAINT, "cc = b3 + b5;", fixed = TRUE)

  addx2 <- generate("addx2")
  expect_match(addx2$MODEL, "Z ON XS (b3);", fixed = TRUE)
  expect_match(addx2$MODEL, "Z ON XY@0 (b4);", fixed = TRUE)
  expect_match(addx2$MODEL, "Z ON YS@0 (b5);", fixed = TRUE)
  expect_match(addx2$MODEL, "Z ON Y (b2);", fixed = TRUE)

  addy2 <- generate("addy2")
  expect_match(addy2$MODEL, "Z ON XS@0 (b3);", fixed = TRUE)
  expect_match(addy2$MODEL, "Z ON YS (b5);", fixed = TRUE)
  expect_match(addy2$MODEL, "Z ON X (b1);", fixed = TRUE)
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
  expect_match(
    paste(readLines(input_path, warn = FALSE), collapse = "\n"),
    "USEVARIABLES ="
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
  workflow <- read_rsa_mplus_model(workflow, output = model_path)

  result <- RSA_mplus(workflow, plot = FALSE)
  expect_equal(
    unname(result$coefficients[c("x", "y", "x2", "xy", "y2", "b0")]),
    c(0.326, 0.166, -0.058, 0.069, -0.085, 0),
    tolerance = 1e-8
  )
  expect_identical(result$workflow, workflow)
  expect_identical(workflow$status, "read")
  expect_identical(workflow$files$output, normalizePath(model_path))

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

  expect_snapshot(
    error = TRUE,
    read_rsa_mplus_model(model_path, output = model_path)
  )
})

test_that("rsa_model fixes zero terms with @0 and never references them", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(5)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))
  regressors <- c(b1 = "X", b2 = "Y", b3 = "XS", b4 = "XY", b5 = "YS")

  for (rsa_model in names(rsa_mplus_model_restrictions)) {
    model <- create_rsa_mplus_model(dat, measurement, rsa_model = rsa_model)
    spec <- rsa_mplus_model_spec(rsa_model)
    expect_identical(model$spec$rsa_model, rsa_model)
    expect_identical(model$spec$metadata$rsa_model, rsa_model)

    fixed <- vapply(
      names(regressors),
      function(label) {
        grepl(
          sprintf("Z ON %s@0 (%s);", regressors[[label]], label),
          model$mplus$MODEL,
          fixed = TRUE
        )
      },
      logical(1)
    )
    expect_identical(names(regressors)[fixed], spec$zero, info = rsa_model)

    # Mplus cannot refer to a parameter fixed with `@0` in MODEL CONSTRAINT.
    constraint <- paste(model$mplus$MODELCONSTRAINT, collapse = "\n")
    for (label in spec$zero) {
      expect_no_match(constraint, paste0("\\b", label, "\\b"), info = rsa_model)
    }
  }

  expect_snapshot({
    for (rsa_model in c("SRRR", "SRSQD", "SRR", "SQD", "IA")) {
      model <- create_rsa_mplus_model(dat, measurement, rsa_model = rsa_model)
      cat("==", rsa_model, "==\n")
      cat(grep(" ON ", strsplit(model$mplus$MODEL, "\n")[[1]], value = TRUE),
        sep = "\n")
      cat(model$mplus$MODELCONSTRAINT, "\n")
    }
  })
})

test_that("restricted surfaces drop exactly the surface parameters fixed at zero", {
  set.seed(6)
  labels <- unname(rsa_mplus_coefficient_labels)
  evaluate <- function(expressions, b) {
    vapply(expressions, function(x) eval(str2lang(x), as.list(b)), numeric(1))
  }

  for (rsa_model in names(rsa_mplus_model_restrictions)) {
    spec <- rsa_mplus_model_spec(rsa_model)
    b <- stats::setNames(numeric(5), labels)
    b[spec$free] <- stats::runif(length(spec$free), -1, 1)
    for (term in names(spec$linear)) {
      weights <- spec$linear[[term]]
      b[[term]] <- sum(weights * b[names(weights)])
    }
    b[names(spec$nonlinear)] <- evaluate(spec$nonlinear, b)

    surface <- vapply(
      rsa_mplus_surface_weights,
      function(weights) sum(weights * b[names(weights)]),
      numeric(1)
    )
    defined <- rsa_mplus_surface_definitions(spec)
    omitted <- setdiff(names(surface), names(defined))
    if (length(defined) > 0L) {
      values <- evaluate(defined, b)
      expect_equal(values, surface[names(defined)], info = rsa_model)
      expect_true(all(abs(values) > 1e-8), info = rsa_model)
    }
    expect_equal(
      unname(surface[omitted]),
      numeric(length(omitted)),
      info = rsa_model
    )
  }
})

test_that("ridge models recover the shift, scaling, and mean effect of the ridge", {
  # Fit-pattern parametrization: Z = bM*M + c1*((S*X + C) - Y)^2 with
  # M = ((S*X + C) + Y)/2, expanded into b1 to b5.
  ridge <- function(c1, S = 1, C = 0, bM = 0) {
    c(
      b1 = bM * S / 2 + 2 * c1 * S * C,
      b2 = bM / 2 - 2 * c1 * C,
      b3 = c1 * S^2,
      b4 = -2 * c1 * S,
      b5 = c1
    )
  }
  evaluate <- function(expressions, b) {
    vapply(expressions, function(x) eval(str2lang(x), as.list(b)), numeric(1))
  }
  cases <- list(
    list("SRRR", c1 = -0.3, S = 1.4, C = 0.25, bM = 0.5),
    list("SRRR", c1 = 0.2, S = 0.7, C = -0.4, bM = -0.3),
    list("SRR", c1 = -0.3, C = 0.25, bM = 0.5),
    list("RR", c1 = -0.3, bM = 0.5),
    list("SRSQD", c1 = -0.3, S = 1.4, C = 0.25),
    list("SSQD", c1 = 0.2, C = -0.4),
    list("SQD", c1 = -0.3)
  )

  for (case in cases) {
    rsa_model <- case[[1]]
    truth <- utils::modifyList(list(S = 1, C = 0, bM = 0), case[-1])
    b <- do.call(ridge, truth)
    spec <- rsa_mplus_model_spec(rsa_model)

    expect_equal(unname(b[spec$zero]), numeric(length(spec$zero)))
    for (term in names(spec$linear)) {
      weights <- spec$linear[[term]]
      expect_equal(b[[term]], sum(weights * b[names(weights)]), info = rsa_model)
    }
    if (length(spec$nonlinear) > 0L) {
      expect_equal(
        evaluate(spec$nonlinear, b),
        b[names(spec$nonlinear)],
        info = rsa_model
      )
    }

    expected <- c(c = truth$C, s = truth$S, bm = truth$bM, a4r = 4 * truth$c1)
    if (length(spec$derived) > 0L) {
      expect_equal(
        evaluate(spec$derived, b),
        expected[names(spec$derived)],
        info = rsa_model
      )
    }
  }
})

test_that("start values are written for the freely estimated terms only", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(7)
  dat <- as.data.frame(matrix(stats::rnorm(700), ncol = 7))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2", "age")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))

  srrr <- create_rsa_mplus_model(
    dat,
    measurement,
    rsa_model = "SRRR",
    start = c(b1 = 0.3, b2 = -0.2, b3 = -0.1, b4 = 0.25, b5 = -0.1)
  )
  expect_match(srrr$mplus$MODEL, "Z ON X*0.3 (b1);", fixed = TRUE)
  expect_match(srrr$mplus$MODEL, "Z ON Y*-0.2 (b2);", fixed = TRUE)
  expect_match(srrr$mplus$MODEL, "Z ON XS (b3);", fixed = TRUE)
  expect_match(srrr$mplus$MODEL, "Z ON XY*0.25 (b4);", fixed = TRUE)
  expect_match(srrr$mplus$MODEL, "Z ON YS*-0.1 (b5);", fixed = TRUE)
  expect_identical(srrr$spec$start, c(b1 = 0.3, b2 = -0.2, b4 = 0.25, b5 = -0.1))

  # Coefficients as RSA_mplus() returns them; fixed terms keep `@0`.
  ia <- create_rsa_mplus_model(
    dat,
    measurement,
    rsa_model = "IA",
    start = c(x = 0.3, y = 0.2, x2 = -0.1, xy = 0.2, y2 = -0.1, b0 = 1)
  )
  expect_match(ia$mplus$MODEL, "Z ON XS@0 (b3);", fixed = TRUE)
  expect_match(ia$mplus$MODEL, "Z ON XY*0.2 (b4);", fixed = TRUE)
  expect_identical(names(ia$spec$start), c("b1", "b2", "b4"))

  expect_snapshot(error = TRUE, {
    create_rsa_mplus_model(dat, measurement, start = c(a1 = 1))
    create_rsa_mplus_model(dat, measurement, start = 1)
  })
})

test_that("unsupported models and reserved labels are rejected", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(8)
  dat <- as.data.frame(matrix(stats::rnorm(700), ncol = 7))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2", "age")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))

  expect_snapshot(error = TRUE, {
    create_rsa_mplus_model(dat, measurement, rsa_model = "weak")
    create_rsa_mplus_model(dat, measurement, rsa_model = "cubic")
    create_rsa_mplus_model(
      dat,
      measurement,
      structural = "Z ~ X + Y + c*age",
      rsa_model = "SRR"
    )
  })
})

test_that("the null model is written without a MODEL CONSTRAINT block", {
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  set.seed(9)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  model <- create_rsa_mplus_model(
    dat,
    list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2")),
    rsa_model = "null"
  )
  expect_null(model$mplus$MODELCONSTRAINT)

  input_path <- file.path(withr::local_tempdir(), "null.inp")
  write_rsa_mplus_model(model, modelout = input_path)
  input <- readLines(input_path)
  expect_false(any(grepl("MODEL CONSTRAINT", input, fixed = TRUE)))
  expect_true(any(grepl("Z ON YS@0 (b5);", input, fixed = TRUE)))
})

test_that("RSA_mplus() takes rsa_model and start values from workflows", {
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
  set.seed(10)
  dat <- as.data.frame(matrix(stats::rnorm(600), ncol = 6))
  names(dat) <- c("x1", "x2", "y1", "y2", "z1", "z2")
  measurement <- list(X = c("x1", "x2"), Y = c("y1", "y2"), Z = c("z1", "z2"))

  # The output is a full model, which contradicts the additive model that the
  # workflow names.
  additive <- create_rsa_mplus_model(dat, measurement, rsa_model = "additive")
  additive <- read_rsa_mplus_model(additive, output = model_path)
  expect_error(
    RSA_mplus(additive, plot = FALSE),
    "`rsa_model = \"additive\"` sets",
    fixed = TRUE
  )
  expect_identical(
    RSA_mplus(additive, rsa_model = "full", plot = FALSE)$rsa_model,
    "full"
  )

  full <- read_rsa_mplus_model(
    create_rsa_mplus_model(dat, measurement),
    output = model_path
  )
  srrr <- create_rsa_mplus_model(
    dat,
    measurement,
    rsa_model = "SRRR",
    start = full
  )
  expect_identical(
    srrr$spec$start,
    c(b1 = 0.326, b2 = 0.166, b4 = 0.069, b5 = -0.085)
  )
  expect_match(srrr$mplus$MODEL, "Z ON X*0.326 (b1);", fixed = TRUE)
})

test_that("Mplus is located by name, including a Windows binary under WSL", {
  skip_on_os("windows")

  expect_identical(resolve_rsa_mplus_command("my-mplus"), "my-mplus")
  expect_error(
    resolve_rsa_mplus_command(c("a", "b")),
    "`Mplus_command` must be a single command or executable path.",
    fixed = TRUE
  )

  bin <- withr::local_tempdir()
  executable <- file.path(bin, "Mplus.exe")
  writeLines("#!/bin/sh", executable)
  Sys.chmod(executable, "755")
  withr::local_envvar(PATH = bin)
  expect_identical(resolve_rsa_mplus_command(), executable)
})

test_that("a lower-cased Mplus output on a case-sensitive file system is explained", {
  skip_on_os("windows")
  skip_if_not_installed("MplusAutomation")
  directory <- withr::local_tempdir()
  input <- file.path(directory, "Upper.inp")
  writeLines("TITLE: a model;", input)
  skip_if(
    file.exists(file.path(directory, "upper.inp")),
    "The file system is not case-sensitive."
  )

  # Stands in for Mplus, which lower-cases the names of the files it writes.
  fake_mplus <- file.path(directory, "fake-mplus")
  writeLines(
    c(
      "#!/bin/sh",
      "out=$(printf '%s' \"$1\" | tr 'A-Z' 'a-z' | sed 's/\\.inp$/.out/')",
      "printf 'Error opening input file\\n' > \"$out\""
    ),
    fake_mplus
  )
  Sys.chmod(fake_mplus, "755")

  error <- expect_error(
    run_rsa_mplus_model(input, Mplus_command = fake_mplus),
    "Mplus did not create the expected output file"
  )
  expect_match(
    conditionMessage(error),
    "give `modelout` a lower-case file name",
    fixed = TRUE
  )
  expect_match(conditionMessage(error), "upper.out", fixed = TRUE)
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
    run = TRUE
  )
  expect_identical(latent$status, "run")
  expect_s3_class(RSA_mplus(latent, plot = FALSE), "rsa_mplus")

  si_lms <- rsa_mplus_workflow(
    data = data.frame(xmean = x, ymean = y, zmean = z),
    measurement = list(X = "xmean", Y = "ymean", Z = "zmean"),
    model_type = "si_lms",
    reliability = c(X = 0.8, Y = 0.8, Z = 0.8),
    modelout = file.path(output_dir, "si_lms.inp"),
    run = TRUE
  )
  expect_identical(si_lms$status, "run")
  expect_s3_class(RSA_mplus(si_lms, plot = FALSE), "rsa_mplus")
})

test_that("pooled centering is a translation that moves a3 and a1 as derived", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  # X sits about 0.8 above Y on the response scale and the surface is curved
  # along the line of incongruence, so the predicted shifts are well above
  # the three-decimal rounding of the Mplus output.
  set.seed(2027)
  n <- 500
  x <- stats::rnorm(n)
  y <- 0.4 * x + stats::rnorm(n, sd = 0.9)
  z <- 0.4 * x + 0.1 * y - 0.15 * x^2 + 0.2 * x * y - 0.1 * y^2 +
    stats::rnorm(n, sd = 0.7)
  item <- function(score, intercept) {
    intercept + score + stats::rnorm(n, sd = 0.5)
  }
  dat <- data.frame(
    x1 = item(x, 4.0), x2 = item(x, 3.9), x3 = item(x, 4.1),
    y1 = item(y, 3.2), y2 = item(y, 3.1), y3 = item(y, 3.3),
    z1 = item(z, 2.5), z2 = item(z, 2.6), z3 = item(z, 2.4)
  )
  measurement <- list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  )
  output_dir <- withr::local_tempdir()

  fit <- function(center) {
    workflow <- rsa_mplus_workflow(
      data = dat,
      measurement = measurement,
      constraints = "surface",
      identification = "effects",
      center = center,
      modelout = file.path(output_dir, paste0(center, ".inp")),
      run = TRUE
    )
    estimates <- workflow$results$parameters$unstandardized
    pick <- function(header, param) {
      estimates$est[estimates$paramHeader == header & estimates$param == param]
    }
    list(
      ll = workflow$results$summaries$LL,
      b = vapply(c("XS", "XY", "YS"), function(p) pick("Z.ON", p), numeric(1)),
      a = vapply(
        c(a1 = "CS", a2 = "CC", a3 = "IS", a4 = "IC", a5 = "A5"),
        function(p) pick("New.Additional.Parameters", p),
        numeric(1)
      ),
      means = c(x = pick("Means", "X"), y = pick("Means", "Y"))
    )
  }

  separate <- fit("none")
  pooled <- fit("pooled")
  delta <- pooled$means[["x"]] - pooled$means[["y"]]

  # Absolute differences: Mplus prints three decimals, and a sum such as a4
  # accumulates the rounding of three printed coefficients.
  expect_close <- function(actual, expected, tolerance) {
    expect_lt(max(abs(actual - expected)), tolerance)
  }

  # Same model, different origin: identical fit and second-order terms.
  expect_close(pooled$ll, separate$ll, 0.01)
  expect_close(sum(pooled$means), 0, 0.01)
  expect_close(pooled$b, separate$b, 0.005)
  expect_close(
    pooled$a[c("a2", "a4", "a5")],
    separate$a[c("a2", "a4", "a5")],
    0.005
  )

  # The shift is large enough to be a real test of the identities.
  expect_gt(abs(delta * separate$a[["a4"]]), 0.1)
  expect_close(
    pooled$a[["a3"]],
    separate$a[["a3"]] - delta * separate$a[["a4"]],
    0.005
  )
  expect_close(
    pooled$a[["a1"]],
    separate$a[["a1"]] - delta * separate$a[["a5"]],
    0.005
  )
})

test_that("the opt-in Mplus integration fits the nested RSA model family", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  # A shifted and rotated rising ridge, Z = bM*M + c1*((S*X + C) - Y)^2,
  # with S = 1.3, C = 0.3, bM = 0.3, and c1 = -0.25.
  set.seed(2028)
  n <- 400
  x <- stats::rnorm(n)
  y <- 0.3 * x + stats::rnorm(n)
  shifted <- 1.3 * x + 0.3
  z <- 0.3 * (shifted + y) / 2 - 0.25 * (shifted - y)^2 +
    stats::rnorm(n, sd = 0.7)
  item <- function(score) score + stats::rnorm(n, sd = 0.4)
  dat <- data.frame(
    x1 = item(x), x2 = item(x), x3 = item(x),
    y1 = item(y), y2 = item(y), y3 = item(y),
    z1 = item(z), z2 = item(z), z3 = item(z)
  )
  measurement <- list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  )
  output_dir <- withr::local_tempdir()
  fit <- function(models, start, file) {
    fit_rsa_mplus_family(
      data = dat,
      measurement = measurement,
      models = models,
      start = start,
      modelout = file.path(output_dir, paste0(file, ".inp"))
    )
  }

  models <- names(rsa_mplus_model_restrictions)
  family <- fit(models, start = NULL, file = "rsa")
  expect_identical(nrow(family$failures), 0L)
  fits <- family$models
  expect_identical(names(fits), models)
  expect_true(all(vapply(fits, function(x) x$status, "") == "run"))
  ll <- vapply(fits, function(x) x$results$summaries$LL, numeric(1))
  npar <- vapply(fits, function(x) x$results$summaries$Parameters, numeric(1))

  # Each restriction removes one free parameter.
  free <- vapply(
    models,
    function(m) length(rsa_mplus_model_spec(m)$free),
    integer(1)
  )
  expect_equal(unname(npar[["full"]] - npar), unname(5 - free))

  # A model never fits better than a model it is nested in. The tolerance
  # absorbs the error of the numerical integration.
  chains <- list(
    c("full", "SRRR", "SRR", "RR", "SQD", "null"),
    c("SRRR", "SRSQD", "SSQD", "SQD"),
    c("SRR", "SSQD", "diff", "null"),
    c("RR", "mean"),
    c("full", "IA", "additive", "mean", "null"),
    c("additive", "diff"),
    c("full", "onlyx2", "onlyx", "null"),
    c("full", "onlyy2", "onlyy", "null"),
    c("additive", "onlyx"),
    c("additive", "onlyy")
  )
  for (chain in chains) {
    for (i in seq_len(length(chain) - 1L)) {
      expect_gt(
        ll[[chain[[i]]]] - ll[[chain[[i + 1L]]]],
        -0.01,
        label = paste(chain[[i]], "vs", chain[[i + 1L]])
      )
    }
  }

  # RSA_mplus() reads every model with the rsa_model of its workflow.
  for (rsa_model in models) {
    result <- RSA_mplus(fits[[rsa_model]], plot = FALSE)
    expect_identical(result$rsa_model, rsa_model)
    expect_true(all(result$coefficients[rsa_mplus_zero_terms[[rsa_model]]] == 0))
  }

  # The ridge parameters follow from the printed weights.
  srrr <- RSA_mplus(fits$SRRR, plot = FALSE)
  b <- stats::setNames(srrr$coefficients[1:5], paste0("b", 1:5))
  derived <- vapply(
    rsa_mplus_model_restrictions$SRRR$derived,
    function(x) eval(str2lang(x), as.list(b)),
    numeric(1)
  )
  printed <- srrr$new_parameters$est[
    match(toupper(names(derived)), srrr$new_parameters$Label)
  ]
  expect_lt(max(abs(printed - derived)), 0.02)

  # Start values from the full model reach the same optimum.
  started <- fit(c("full", "SRRR"), start = "full", file = "started")
  expect_identical(
    started$models$SRRR$spec$start,
    stats::setNames(
      RSA_mplus(fits$full, plot = FALSE)$coefficients[c("x", "y", "xy", "y2")],
      c("b1", "b2", "b4", "b5")
    )
  )
  expect_lt(
    abs(started$models$SRRR$results$summaries$LL - ll[["SRRR"]]),
    0.01
  )

  # The comparison finds no nesting violations and ranks the true model
  # first. AICc agrees with MplusAutomation's, which starts from the AIC.
  expect_no_warning(comparison <- compare_rsa_mplus_models(family))
  expect_setequal(comparison$model, models)
  expect_true(all(comparison$converged))
  expect_identical(comparison$model[[1]], "SRRR")
  aicc <- vapply(fits, function(x) x$results$summaries$AICC, numeric(1))
  expect_lt(max(abs(comparison$AICc - aicc[comparison$model])), 0.01)
  expect_equal(sum(comparison$weight), 1)
  expect_equal(
    comparison$lrt_df,
    ifelse(comparison$model == "full", NA, npar[["full"]] - npar[comparison$model])
  )
  r2 <- vapply(
    fits,
    function(x) with(x$results$parameters$r2, est[param == "Z"]),
    numeric(1)
  )
  expect_equal(comparison$R2, unname(r2[comparison$model]))
  expect_true(all(comparison$flags %in% c("", "cd <= 0")))
})

test_that("the opt-in Mplus integration flags ridge models that drift on a plane", {
  skip_if(Sys.getenv("RSAHELPERS_RUN_MPLUS_TESTS") != "true")
  skip_if_not_installed("tidySEM")
  skip_if_not_installed("MplusAutomation")

  # A plane with no ridge: the rotated models push the ridge out of the data
  # and SRSQD, which has no mean effect, drifts to b5 = 0.
  set.seed(2029)
  n <- 400
  x <- stats::rnorm(n)
  y <- 0.3 * x + stats::rnorm(n)
  z <- 0.4 * x + 0.1 * y + stats::rnorm(n, sd = 0.7)
  item <- function(score) score + stats::rnorm(n, sd = 0.4)
  dat <- data.frame(
    x1 = item(x), x2 = item(x), x3 = item(x),
    y1 = item(y), y2 = item(y), y3 = item(y),
    z1 = item(z), z2 = item(z), z3 = item(z)
  )
  family <- fit_rsa_mplus_family(
    data = dat,
    measurement = list(
      X = c("x1", "x2", "x3"),
      Y = c("y1", "y2", "y3"),
      Z = c("z1", "z2", "z3")
    ),
    models = c("full", "SRRR", "SRSQD", "additive"),
    modelout = file.path(withr::local_tempdir(), "plane.inp")
  )
  # Drifting ridge models stop short of their optimum, so depending on the
  # Mplus version SRSQD can end slightly above SRRR, which is warned about.
  comparison <- suppressWarnings(compare_rsa_mplus_models(family))
  flags <- stats::setNames(comparison$flags, comparison$model)

  expect_identical(comparison$model[[1]], "additive")
  expect_match(flags[["SRSQD"]], "b5 = 0", fixed = TRUE)
  expect_match(flags[["SRSQD"]], "SD from the data", fixed = TRUE)
  expect_match(flags[["SRRR"]], "SD from the data", fixed = TRUE)
  expect_identical(flags[["additive"]], "")
})
