# Fit an Edwards-Parry congruence spline surface

Estimate one- or two-seam spline regression surfaces for congruence
hypotheses using nonlinear least squares. The one-seam model matches
Equation 39 in Edwards and Parry (2018); the two-seam model extends it
with a second estimated seam term.

## Usage

``` r
fit_spline_congruence(
  formula,
  data,
  n_seams = c(1L, 2L),
  center = c("pooled", "variablewise", "none"),
  scale = c("pooled", "none"),
  hinge_offset = 1,
  starts = NULL,
  fix = NULL,
  multistart = TRUE,
  control = list(maxit = 10000),
  prefer_nlslm = TRUE
)
```

## Arguments

- formula:

  A two-sided formula of the form `z ~ x * y`. The response is mapped to
  Z, the first predictor to X, and the second predictor to Y. The `*`
  declares the two predictor roles rather than an ordinary interaction.

- data:

  A data frame containing the Z, X, and Y variables.

- n_seams:

  Number of seams to estimate. Supported values are `1` and `2`.

- center, scale:

  Centering and scaling of the component variables, passed to the
  internal model-frame builder. Defaults to pooled centering and pooled
  scaling (the congruence-safe `RSA` convention). See
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md)
  for the accepted values and the logical backward-compatibility
  mapping.

- hinge_offset:

  Numeric half-width used for the fixed-seam starting values, in
  working-scale units. See
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md).

- starts:

  Optional named numeric vector or list of named numeric vectors with
  starting values. One-seam names are `b0`, `b1`, `b2`, `b3`, `c0`, and
  `c1`. Two-seam names are `b0`, `b1`, `b2`, `b3`, `b4`, `c10`, `c11`,
  `c20`, and `c21`.

- fix:

  Optional named numeric vector of parameters to hold constant rather
  than estimate, for example `c(c0 = 0, c1 = 1)` to fix the one-seam
  line of congruence at `Y = X`. Held-constant parameters are dropped
  from the optimization, excluded from the residual degrees of freedom,
  and given zero covariance, yielding a constrained fit suitable for a
  nested comparison against the unconstrained spline (see
  [`compare_spline_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_spline_models.md)).

- multistart:

  Logical. If `TRUE`, fit from Stata-style constrained and unconstrained
  starts plus small jittered variants and keep the lowest-RSS converged
  solution.

- control:

  List of solver control options. `maxit` is accepted for convenience
  and translated to `minpack.lm`'s `maxiter`.

- prefer_nlslm:

  Logical. If `TRUE` and `minpack.lm` is installed, use
  [`minpack.lm::nlsLM()`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html).
  Otherwise use a base-R [`optim()`](https://rdrr.io/r/stats/optim.html)
  fallback.

## Value

A `congruence_spline` object, a list containing coefficients, covariance
matrix, fitted values, residuals, RSS, R-squared, residual degrees of
freedom, solver metadata, and all start attempts. Any held-constant
parameters are recorded in the `fixed` element.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_spline_congruence(satisfaction ~ x * y, dat, n_seams = 1)
coef(fit)
spline_surface_features(fit)
} # }
```
