# splinecongruence

`splinecongruence` provides a small R implementation of Edwards and Parry's
spline regression procedure for congruence research.

The package includes helpers to fit one- and two-seam nonlinear spline models,
fit OLS comparison models, run Wald tests, compute bootstrap intervals, and
plot fitted response surfaces with base R graphics.

```r
library(splinecongruence)

workshop <- haven::read_dta(
  system.file("extdata", "spline.dta", package = "splinecongruence")
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
