# Fit OLS comparison models for congruence analyses

Fit the absolute-difference, linear, one-break piecewise, constrained
one-seam piecewise, and optionally fixed two-seam OLS models used as
comparisons and starting-value sources for spline regression.

## Usage

``` r
fit_piecewise_congruence(
  formula,
  data,
  center = c("pooled", "variablewise", "none"),
  scale = c("pooled", "none"),
  hinge_offset = 1,
  n_seams = c(1L, 2L),
  constrained = TRUE
)
```

## Arguments

- formula:

  A two-sided formula of the form `z ~ x * y`. The response is mapped to
  Z, the first predictor to X, and the second predictor to Y. The `*`
  declares the two predictor roles rather than an ordinary interaction.

- data:

  A data frame containing the Z, X, and Y variables.

- center, scale:

  Centering and scaling of the component variables, passed to
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md).
  Defaults to pooled centering and pooled scaling.

- hinge_offset:

  Numeric half-width, in working-scale units, for the fixed two-seam OLS
  seam lines `Y = X ± hinge_offset`. See
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md).

- n_seams:

  Number of seams to include in comparison models. If `2`, add the fixed
  two-seam OLS model using seam lines `Y = X + hinge_offset` and
  `Y = X - hinge_offset`.

- constrained:

  Logical. If `FALSE`, omit the constrained one-seam OLS model from the
  returned list.

## Value

A named list of `lm` objects with class `congruence_piecewise`.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- fit_piecewise_congruence(satisfaction ~ x * y, dat)
tidy_piecewise_summary(fits)
} # }
```
