# rsahelpers 0.1.0

- The package has been renamed from `splinecongruence` to `rsahelpers` to
  support a broader collection of response surface analysis utilities.
- `classify_yao_ma()` classifies Mplus or parameter-table response surfaces
  using significance-based rules or user-supplied practical-equivalence
  margins, with diagnostic partial and indeterminate results.
- `create_rsa_mplus_model()` and the related write, run, read, and workflow
  helpers now generate full latent and reliability-corrected SI-LMS response
  surface models, including `XWITH` interactions, surface constraints, and an
  explicit `USEVARIABLES` statement.
- `predictor_coverage()` and `plot_predictor_coverage()` summarize and display
  sparse regions in the joint distribution of commensurate X/Y scale scores.
- `read_rsa_mplus_model()` can attach a cached or relocated Mplus output to an
  existing workflow while preserving its role metadata.
- `RSA_mplus()` is now the canonical implementation for extracting polynomial
  coefficients from Mplus output and optionally plotting them with `RSA`; it
  now also accepts fitted workflow and `mplusObject` inputs and infers generated
  model labels automatically. Three-dimensional plots annotate `a1` through
  `a5` with APA-style significance stars when Mplus constraint p-values are
  available. With `ESTIMATOR = BAYES`, the 95% credibility interval is used for
  the significance decision instead of the one-tailed posterior p-value, so a
  single `*` marks surface parameters whose interval excludes zero. The
  interval bounds are also returned in `new_parameters` as `lower_2.5ci` and
  `upper_2.5ci`.

- rsahelpers no longer declares or installs `franzpak`. The two packages
  previously listed each other in `Suggests` and `Remotes`, which made their
  dependency resolution circular and coupled their CI. The dependency is now
  one-directional: franzpak's deprecated `RSA_mplus()` forwards here, and the
  optional franzpak coefficient table in the Mplus vignette is guarded by
  `requireNamespace()`.

## Centering, scaling, and seams

- `fit_piecewise_congruence()`, `fit_spline_congruence()`, and
  `prepare_congruence_data()` now use the consistent formula interface
  `Z ~ X * Y` instead of separate `wanted`, `actual`, and `outcome` arguments.
- **Pooled centering and scaling are now the default** (`center = "pooled"`,
  `scale = "pooled"`), following the `RSA` package, so the `X = Y` congruence
  interpretation is preserved out of the box. `center`/`scale` accept
  `"pooled"`, `"variablewise"`, and `"none"`; logical `center` is still accepted
  for backward compatibility (`TRUE` = variablewise, `FALSE` = none, with no
  scaling unless `scale` is given explicitly).
- `prepare_congruence_data()`, `fit_spline_congruence()`, and
  `fit_piecewise_congruence()` gain a `hinge_offset` argument that parameterizes
  the previously hardcoded `Y = X ± 1` fixed two-seam hinges. The offset is in
  working-scale units (i.e. SD units under pooled scaling).
- `fit_spline_congruence()` gains a `fix` argument to hold parameters constant
  (e.g. `fix = c(c0 = 0, c1 = 1)`) for a constrained nested fit.

## Inference and model comparison

- `bootstrap_spline()` now routes through `boot::boot.ci()` and supports
  `type = "perc"`, `"basic"`, `"norm"`, and `"bca"`. The bootstrap failure rate
  is reported via a warning and a `fail_rate` attribute.
- New `compare_spline_models()` and `anova()` method build nested-model
  comparison tables (ΔR², F on ΔRSS, df, p, AIC) with a guard that flags when a
  larger model fits worse than the simpler model nested within it.
- `select_spline_congruence()` selects a linear, fixed-LOC, or free one-seam
  surface using residual-bootstrap likelihood-ratio tests and arm-coverage
  diagnostics; absolute-difference, piecewise, and two-seam fits are retained
  as benchmarks or sensitivity analyses.
- New `logLik()`, `nobs()`, and `summary()` methods for `congruence_spline`,
  enabling `AIC()`/`BIC()` and a coefficient table with delta-method standard
  errors, t statistics, and p-values.
- `spline_tests()` omits invalid Wald tests of seam existence and seam count;
  its remaining seam-location and shape tests are conditional on a supported
  one-seam surface. `summary()` and `spline_tests()` warn (controllable via
  `warn_seam`) that
  delta-method standard errors for seam parameters are approximate because the
  Jacobian is evaluated at the non-differentiable seam; the bootstrap is
  recommended for seam inference.
- `spline_surface_features()` now reports two-seam section slopes and a
  seam-crossing diagnostic (`seams_cross`, `crossing_x`, `n_sections`).
- `tidy_piecewise_summary()` now includes `std.error`, `statistic`, and `p.value`.

## Plotting

- `plot_spline_surface()` now shades facets by predicted outcome on a
  high-contrast diverging palette (RSA-style), uses larger tiles by default, and
  adds `show_fit_line` to optionally drop the theoretical `X = Y` line.
- New `plot_spline_contour()` draws a 2D `ggplot2` heatmap with contour lines and
  seam overlays (requires `ggplot2`).

## Initial release notes

- Initial package version with Edwards-Parry one- and two-seam spline
  regression helpers, OLS comparison models, Wald tests, bootstrap intervals,
  and base R surface plotting.
