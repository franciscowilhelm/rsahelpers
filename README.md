# rsahelpers

`rsahelpers` bundles response surface analysis (RSA) tools used across congruence
research projects. It includes helpers for working with polynomial response surfaces estimated
in Mplus and an R implementation of Edwards and Parry's spline regression procedure (Edwards and Parry, 2018). 

Response surface analyses build on the excellent `RSA` package which estimates models via lavaan. The rsahelpers aims to re-implement and enhance the RSA package's functions in Mplus to make use of its capabilities for estimating models in a latent variable framework with latent interactions.

## Installation

Install the development version from GitHub:

```r
pak::pak("franciscowilhelm/rsahelpers")
```


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

The restricted surface models of the `RSA` package (`"SRRR"`, `"SRR"`, `"RR"`,
`"SQD"`, `"IA"`, `"additive"`, and so on) are generated with `rsa_model`.
`fit_rsa_mplus_family()` fits a set of them from one specification, and
`compare_rsa_mplus_models()` ranks them by AICc and tests them against the full
model:

```r
family <- fit_rsa_mplus_family(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  identification = "invariance",
  center = "pooled",
  modelout = "rsa.inp"
)
comparison <- compare_rsa_mplus_models(family)
RSA_mplus(family$models[[comparison$model[[1]]]])
```

`identification` and `center` put the X and Y factors on one metric with one
origin, so that the line X = Y is equality on the response scale;
`fit_rsa_mplus_invariance()` tests the measurement invariance this assumes.
`rsa_mplus_draws()` gives Monte Carlo intervals for the surface and ridge
parameters, and `rsa_mplus_supported_range()` the region of the surface that
the data support.

See `vignette("mplus-response-surfaces", package = "rsahelpers")` for the
complete workflow, including character syntax, custom Mplus blocks, SI-LMS,
execution, result reading, constraints, and plotting.

The same material is available in the [Mplus response surfaces
guide](https://franciscowilhelm.github.io/rsahelpers/articles/mplus-response-surfaces.html).

## Edwards-Parry spline regression

Fit one- and two-seam nonlinear spline models, OLS comparison models, Wald
tests, and bootstrap intervals, and plot fitted response surfaces:

```r
library(rsahelpers)

workshop <- haven::read_dta(
  system.file("extdata", "spline.dta", package = "rsahelpers")
)

fit <- fit_spline_congruence(
  JOBSAT ~ ATHW * ATHH,
  data = workshop
)

coef(fit)
spline_surface_features(fit)
spline_tests(fit)
```

Spline functions use the formula convention `Z ~ X * Y`: the response is Z,
the first predictor is X, and the second predictor is Y. Here, `ATHW` is X and
`ATHH` is Y. The raw predictors are pooled-centered and pooled-scaled by
default. The `*` declares their roles in the spline model; it does not add an
ordinary linear interaction term.

See the [Spline response surfaces
guide](https://franciscowilhelm.github.io/rsahelpers/articles/spline-response-surfaces.html)
for model comparison, interpretation, visualization, and bootstrap inference.
