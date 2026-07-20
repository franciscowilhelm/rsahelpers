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

`RSA_mplus()` extracts polynomial coefficients and optional `MODEL CONSTRAINT`
parameters from an Mplus output file. Set `plot = TRUE` to pass the polynomial
coefficients to `RSA::plotRSA()`.

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
