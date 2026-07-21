# Summarize a fitted congruence spline

Assemble a coefficient table with delta-method standard errors, t
statistics, and p-values alongside fit statistics.

## Usage

``` r
# S3 method for class 'congruence_spline'
summary(object, warn_seam = TRUE, ...)
```

## Arguments

- object:

  A `congruence_spline` object.

- warn_seam:

  Logical. If `TRUE` (default), warn that delta-method standard errors
  for seam parameters are approximate. See
  [`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md).

- ...:

  Unused.

## Value

An object of class `summary.congruence_spline`: a list with a
`coefficients` data frame and the fit summary statistics.
