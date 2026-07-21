# Changelog

## rsahelpers 0.0.0.9000

- The package has been renamed from `splinecongruence` to `rsahelpers`
  to support a broader collection of response surface analysis
  utilities.
- [`classify_yao_ma()`](https://franciscowilhelm.github.io/rsahelpers/reference/classify_yao_ma.md)
  classifies Mplus or parameter-table response surfaces using
  significance-based rules or user-supplied practical-equivalence
  margins, with diagnostic partial and indeterminate results.
- [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  and the related write, run, read, and workflow helpers now generate
  full latent and reliability-corrected SI-LMS response surface models,
  including `XWITH` interactions, surface constraints, and an explicit
  `USEVARIABLES` statement.
- [`predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/predictor_coverage.md)
  and
  [`plot_predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_predictor_coverage.md)
  summarize and display sparse regions in the joint distribution of
  commensurate X/Y scale scores.
- [`read_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  can attach a cached or relocated Mplus output to an existing workflow
  while preserving its role metadata.
- [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  is now the canonical implementation for extracting polynomial
  coefficients from Mplus output and optionally plotting them with
  `RSA`; it now also accepts fitted workflow and `mplusObject` inputs
  and infers generated model labels automatically. Three-dimensional
  plots annotate `a1` through `a5` with APA-style significance stars
  when Mplus constraint p-values are available.

### Centering, scaling, and seams

- [`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md),
  [`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md),
  and
  [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md)
  now use the consistent formula interface `Z ~ X * Y` instead of
  separate `wanted`, `actual`, and `outcome` arguments.
- **Pooled centering and scaling are now the default**
  (`center = "pooled"`, `scale = "pooled"`), following the `RSA`
  package, so the `X = Y` congruence interpretation is preserved out of
  the box. `center`/`scale` accept `"pooled"`, `"variablewise"`, and
  `"none"`; logical `center` is still accepted for backward
  compatibility (`TRUE` = variablewise, `FALSE` = none, with no scaling
  unless `scale` is given explicitly).
- [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md),
  [`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md),
  and
  [`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md)
  gain a `hinge_offset` argument that parameterizes the previously
  hardcoded `Y = X ± 1` fixed two-seam hinges. The offset is in
  working-scale units (i.e. SD units under pooled scaling).
- [`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md)
  gains a `fix` argument to hold parameters constant
  (e.g. `fix = c(c0 = 0, c1 = 1)`) for a constrained nested fit.

### Inference and model comparison

- [`bootstrap_spline()`](https://franciscowilhelm.github.io/rsahelpers/reference/bootstrap_spline.md)
  now routes through
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) and
  supports `type = "perc"`, `"basic"`, `"norm"`, and `"bca"`. The
  bootstrap failure rate is reported via a warning and a `fail_rate`
  attribute.
- New
  [`compare_spline_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_spline_models.md)
  and [`anova()`](https://rdrr.io/r/stats/anova.html) method build
  nested-model comparison tables (ΔR², F on ΔRSS, df, p, AIC) with a
  guard that flags when a larger model fits worse than the simpler model
  nested within it.
- New [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for
  `congruence_spline`, enabling
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  and a coefficient table with delta-method standard errors, t
  statistics, and p-values.
- [`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)
  and [`summary()`](https://rdrr.io/r/base/summary.html) now warn
  (controllable via `warn_seam`) that delta-method standard errors for
  seam parameters are approximate because the Jacobian is evaluated at
  the non-differentiable seam; the bootstrap is recommended for seam
  inference.
- [`spline_surface_features()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_surface_features.md)
  now reports two-seam section slopes and a seam-crossing diagnostic
  (`seams_cross`, `crossing_x`, `n_sections`).
- [`tidy_piecewise_summary()`](https://franciscowilhelm.github.io/rsahelpers/reference/tidy_piecewise_summary.md)
  now includes `std.error`, `statistic`, and `p.value`.

### Plotting

- [`plot_spline_surface()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_spline_surface.md)
  now shades facets by predicted outcome on a high-contrast diverging
  palette (RSA-style), uses larger tiles by default, and adds
  `show_fit_line` to optionally drop the theoretical `X = Y` line.
- New
  [`plot_spline_contour()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_spline_contour.md)
  draws a 2D `ggplot2` heatmap with contour lines and seam overlays
  (requires `ggplot2`).

### Initial release notes

- Initial package version with Edwards-Parry one- and two-seam spline
  regression helpers, OLS comparison models, Wald tests, bootstrap
  intervals, and base R surface plotting.
