# Summarize OLS comparison models

Convert the list returned by
[`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md)
to a compact data frame with estimates and model R-squared values.

## Usage

``` r
tidy_piecewise_summary(fits)
```

## Arguments

- fits:

  A `congruence_piecewise` object or named list of `lm` objects.

## Value

A data frame with columns `model`, `term`, `estimate`, `std.error`,
`statistic`, `p.value`, and `r.squared`.
