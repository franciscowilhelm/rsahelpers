# Package index

## RSA with Mplus

Build Mplus response-surface models, extract their polynomial surfaces,
and classify fitted results.

- [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  [`write_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  [`run_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  [`read_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  : Build, write, run, and read RSA models for Mplus
- [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  : Extract and plot RSA surfaces from Mplus output
- [`classify_yao_ma()`](https://franciscowilhelm.github.io/rsahelpers/reference/classify_yao_ma.md)
  : Classify a response surface using the Yao-Ma typology

## Predictor diagnostics

Check whether combinations of the two predictors are supported by the
observed data.

- [`predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/predictor_coverage.md)
  : Diagnose joint predictor coverage for response-surface analysis
- [`plot_predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_predictor_coverage.md)
  : Plot joint predictor coverage

## Spline response surfaces

Prepare data and fit, compare, summarize, test, bootstrap, and plot
Edwards-Parry spline models.

- [`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md)
  : Prepare variables for congruence regression models
- [`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md)
  : Fit an Edwards-Parry congruence spline surface
- [`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md)
  : Fit OLS comparison models for congruence analyses
- [`compare_spline_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_spline_models.md)
  [`anova(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_spline_models.md)
  : Compare nested congruence models
- [`spline_surface_features()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_surface_features.md)
  : Compute derived surface features for a fitted spline model
- [`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)
  : Run Edwards-Parry-style Wald tests for a spline model
- [`bootstrap_spline()`](https://franciscowilhelm.github.io/rsahelpers/reference/bootstrap_spline.md)
  : Bootstrap coefficients and surface features for a spline model
- [`plot_spline_surface()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_spline_surface.md)
  : Plot a fitted congruence spline surface in 3D
- [`plot_spline_contour()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_spline_contour.md)
  : Plot a fitted congruence spline surface as a 2D contour map
- [`tidy_piecewise_summary()`](https://franciscowilhelm.github.io/rsahelpers/reference/tidy_piecewise_summary.md)
  : Summarize OLS comparison models
- [`coef(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/coef.congruence_spline.md)
  : Extract congruence spline coefficients
- [`logLik(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/logLik.congruence_spline.md)
  : Log-likelihood of a fitted congruence spline
- [`nobs(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/nobs.congruence_spline.md)
  : Number of observations in a fitted congruence spline
- [`print(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/print.congruence_spline.md)
  : Print a congruence spline model
- [`summary(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/summary.congruence_spline.md)
  : Summarize a fitted congruence spline
- [`vcov(`*`<congruence_spline>`*`)`](https://franciscowilhelm.github.io/rsahelpers/reference/vcov.congruence_spline.md)
  : Extract a congruence spline covariance matrix

## Package

Package overview.

- [`rsahelpers`](https://franciscowilhelm.github.io/rsahelpers/reference/rsahelpers-package.md)
  [`rsahelpers-package`](https://franciscowilhelm.github.io/rsahelpers/reference/rsahelpers-package.md)
  : rsahelpers: Helpers for Response Surface Analysis
