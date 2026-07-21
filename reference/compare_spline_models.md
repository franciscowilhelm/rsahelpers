# Compare nested congruence models

Build a sequential nested-model comparison table for any mix of `lm`
comparison models (from
[`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md))
and spline fits (from
[`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md)).
Each adjacent pair contributes a change in R-squared and an F test on
the change in residual sum of squares, alongside per-model AIC.

## Usage

``` r
compare_spline_models(...)

# S3 method for class 'congruence_spline'
anova(object, ...)
```

## Arguments

- ...:

  Two or more fitted models in nested order. Names, if supplied, are
  used as row labels.

- object:

  A `congruence_spline` object (the first model).

## Value

A data frame with one row per model: `model`, `npar`, `df.residual`,
`rss`, `r.squared`, `AIC`, and the comparison columns `df`, `deltaR2`,
`F`, and `p.value` (the first row's comparison columns are `NA`).

## Details

Pass the models in increasing complexity, for example the Edwards-Parry
chain absolute-difference \\\subset\\ linear \\\subset\\ one-break
piecewise \\\subset\\ spline, or a one-seam spline followed by a
two-seam spline. If a larger model has a higher RSS than the simpler
model nested within it (a sign of a nonlinear local minimum), its F and
p-value are returned as `NA` and a warning recommends keeping the
simpler model.

## Examples

``` r
if (FALSE) { # \dontrun{
ols <- fit_piecewise_congruence(satisfaction ~ x * y, dat)
spline <- fit_spline_congruence(satisfaction ~ x * y, dat)
compare_spline_models(
  absdiff = ols$absolute_difference,
  linear = ols$linear,
  piecewise = ols$one_break,
  spline = spline
)
} # }
```
