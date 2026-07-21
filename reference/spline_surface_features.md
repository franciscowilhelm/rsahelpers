# Compute derived surface features for a fitted spline model

Derive interpretable features from the fitted surface, including side
slopes, along-seam slopes, symmetry expressions, and seam shifts along
lines parallel to `Y = -X`.

## Usage

``` r
spline_surface_features(fit, lines = c(-1, 0, 1))
```

## Arguments

- fit:

  A `congruence_spline` object.

- lines:

  Numeric vector giving the `k` values where lines of interest cross
  `Y = X` at `(k, k)`. For one-seam models, each line is `Y = 2k - X`.

## Value

A named numeric vector of derived quantities. One-seam models return the
full Edwards-Parry set currently implemented. Two-seam models return the
seam shifts along `Y = -X`, the section slopes in `x` and `y` for each
active-hinge region (base, seam1-only, seam2-only, both), and a crossing
diagnostic: `seams_cross` (1 if the seams cross within the observed `x`
range, else 0), `crossing_x`, and `n_sections` (3 if the seams cross
inside the data, else 4).

## Examples

``` r
if (FALSE) { # \dontrun{
spline_surface_features(fit, lines = c(-1, 0, 1))
} # }
```
