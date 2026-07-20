# rsahelpers

`rsahelpers` bundles response surface analysis tools used across congruence
research projects. It includes an R implementation of Edwards and Parry's
spline regression procedure (Edwards and Parry, 2018) and helpers for working
with polynomial response surfaces estimated in Mplus.

## Installation

Install the development version from GitHub:

```r
pak::pak("franciscowilhelm/rsahelpers")
```

## Edwards-Parry spline regression

Fit one- and two-seam nonlinear spline models, OLS comparison models, Wald
tests, and bootstrap intervals, and plot fitted response surfaces:

```r
library(rsahelpers)

workshop <- haven::read_dta(
  system.file("extdata", "spline.dta", package = "rsahelpers")
)

fit <- fit_spline_congruence(
  workshop,
  wanted = "ATHWC",
  actual = "ATHHC",
  outcome = "JOBSAT",
  center = FALSE
)

coef(fit)
surface_features(fit)
spline_tests(fit)
```

By field convention, `wanted` or ideal values are mapped to `x`, and `actual`
or current values are mapped to `y`.

## Mplus response surfaces

The RSA-Mplus workflow can generate a latent measurement model with tidySEM,
add the LMS interactions and response-surface constraints, and export the
corresponding Mplus input and data files:

```r
model <- rsa_mplus_workflow(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  structural = list(Z = c("X", "Y")),
  modelout = "rsa-latent.inp"
)
```

By default this writes reproducible `.inp` and `.dat` files without running
Mplus. Set `run = TRUE` to execute and read the model. The steps are also
available separately through `create_rsa_mplus_model()`,
`write_rsa_mplus_model()`, `run_rsa_mplus_model()`, and
`read_rsa_mplus_model()`.

```r
fitted <- run_rsa_mplus_model(model)
```

Reliability-corrected single-indicator LMS models use prepared scale scores and
reliability estimates:

```r
si_lms <- create_rsa_mplus_model(
  data = scores,
  measurement = list(X = "xmean", Y = "ymean", Z = "zmean"),
  model_type = "si_lms",
  reliability = c(X = 0.82, Y = 0.79, Z = 0.88)
)
```

`RSA_mplus()` extracts polynomial coefficients and optional `MODEL CONSTRAINT`
parameters from a fitted workflow, an Mplus object, or an output file. Generated
workflows retain the labels needed for plotting, so no variable mapping is
required:

```r
surface <- RSA_mplus(fitted, plot = FALSE)
surface$coefficients
```

For arbitrary output files, specify the labels explicitly:

```r
surface <- RSA_mplus(
  system.file("extdata", "congruence_sim.out", package = "rsahelpers"),
  outcome = "Z",
  pred_x = "X",
  pred_y = "Y",
  pred_x2 = "XS",
  pred_xy = "XY",
  pred_y2 = "YS",
  plot = FALSE
)

surface$coefficients
```
