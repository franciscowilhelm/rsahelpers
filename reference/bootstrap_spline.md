# Bootstrap coefficients and surface features for a spline model

Refit the spline model to bootstrap resamples and compute confidence
intervals for coefficients and derived surface features via
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html). Edwards
and Parry (2018) recommend the bootstrap over delta-method inference
because the seam parameters sit at a non-differentiable kink (see
[`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)).

## Usage

``` r
bootstrap_spline(
  fit,
  R = 10000,
  conf = 0.95,
  type = c("perc", "basic", "norm", "bca")
)
```

## Arguments

- fit:

  A `congruence_spline` object.

- R:

  Number of bootstrap resamples. Edwards and Parry used 10,000 in the
  example analyses; smaller values are useful for smoke tests.

- conf:

  Confidence level for intervals.

- type:

  Interval type passed to
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html):
  `"perc"` (percentile, the default), `"basic"`, `"norm"` (normal
  approximation), or `"bca"` (bias-corrected and accelerated). If
  `boot.ci()` cannot compute the requested interval for a degenerate
  statistic, the percentile interval is used as a fallback.

## Value

A data frame with columns `term`, `estimate`, `lower`, and `upper`. The
proportion of bootstrap resamples whose refit failed (and were dropped)
is attached as the `fail_rate` attribute; a non-zero failure rate also
triggers a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
bootstrap_spline(fit, R = 100, type = "bca")
} # }
```
