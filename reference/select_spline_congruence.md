# Select a continuous congruence surface

Fit a linear plane, a one-seam spline fixed to the line of congruence,
and a freely located one-seam spline, then select among them using
residual- bootstrap likelihood-ratio tests. Absolute-difference and
piecewise models are retained as benchmarks, while an optional two-seam
spline is retained as a sensitivity analysis and is never selected as
the focal model.

## Usage

``` r
select_spline_congruence(
  formula,
  data,
  center = c("pooled", "variablewise", "none"),
  scale = c("pooled", "none"),
  hinge_offset = 1,
  R = 1999,
  alpha = 0.05,
  min_arm_n = 30L,
  min_arm_prop = 0.1,
  max_fail = 0.1,
  seed = NULL,
  include_two_seam = TRUE
)
```

## Arguments

- formula, data:

  Model formula and data accepted by
  [`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md).

- center, scale:

  Congruence-preserving transformation options passed to
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md).

- hinge_offset:

  Working-scale offset used for fixed two-seam benchmarks.

- R:

  Number of null-bootstrap resamples for each structural comparison.

- alpha:

  Unadjusted decision threshold.

- min_arm_n, min_arm_prop:

  Minimum count and proportion required on each side of the freely
  estimated seam.

- max_fail:

  Maximum acceptable failed-refit proportion for a bootstrap structural
  test.

- seed:

  Optional integer seed. The LOC test uses `seed + 1`.

- include_two_seam:

  Whether to fit a two-seam sensitivity model.

## Value

An object of class `spline_selection` containing `fits`, `tests`,
`diagnostics`, `selected_model`, `selected_fit`, `status`, and `reason`.

## Details

The seam-existence comparison is non-regular because the seam location
is unidentified under a linear surface. It therefore uses a null
residual bootstrap rather than the ordinary F or Wald reference
distribution. The fixed-versus-free seam comparison is performed only
after the seam is supported and both arms contain enough observations.
