# Run Edwards-Parry-style Wald tests for a spline model

Compute normal-theory delta and Wald tests for quantities used to
interpret a supported one-seam congruence surface. The joint tests cover
the absolute-difference restrictions and whether the seam equals
`Y = X`. Structural tests of whether any seam exists or whether a second
seam is needed are deliberately omitted because seam locations are
unidentified under those null hypotheses. Use
[`select_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/select_spline_congruence.md)
for the residual-bootstrap seam-existence test.

## Usage

``` r
spline_tests(fit, warn_seam = TRUE)
```

## Arguments

- fit:

  A `congruence_spline` object.

- warn_seam:

  Logical. If `TRUE` (default), warn that the delta-method standard
  errors underlying these tests are approximate for quantities involving
  the seam parameters, because the numeric Jacobian is evaluated at the
  non-differentiable seam. Edwards and Parry (2018) recommend
  [`bootstrap_spline()`](https://franciscowilhelm.github.io/rsahelpers/reference/bootstrap_spline.md)
  for inference on seam features.

## Value

A list with two data frames: `scalar` for one-degree-of-freedom
delta-method tests and `joint` for multi-constraint Wald F tests.

## Examples

``` r
if (FALSE) { # \dontrun{
tests <- spline_tests(fit)
tests$joint
} # }
```
