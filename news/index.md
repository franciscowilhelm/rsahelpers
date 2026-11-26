# Changelog

## rsahelpers 0.2.0

This release brings the model family of the `RSA` package to the latent
Mplus workflow: restricted response surface models, a runner and a
comparison table for the family, a shared metric for the X and Y
factors, and intervals for derived surface parameters. The defaults
reproduce the output of 0.1.0.

### Restricted response surface models

- [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  and
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  gain `rsa_model`, which generates the restricted second-order models
  of [`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html) as latent
  Mplus models: the fit-pattern models `"SRRR"`, `"SRR"`, `"RR"`,
  `"SRSQD"`, `"SSQD"`, and `"SQD"`, as well as `"IA"`, `"additive"`,
  `"mean"`, `"diff"`, `"onlyx"`, `"onlyy"`, `"onlyx2"`, `"onlyy2"`, and
  `"null"`. Terms set to zero are fixed with `@0`; the other
  restrictions and the ridge parameters `c`, `s`, `bm`, and `a4r` (RSA’s
  C, S, `meaneffect`, and `a4.rescaled`) are written in
  `MODEL CONSTRAINT`. The rotated ridge models express `b3` as a
  function of the other weights, so a single run fits ridges and
  troughs. Surface and principal-axis parameters that a model fixes or
  leaves undefined are not generated. `start` takes start values, for
  example from the fitted full model.

- `rsa_model` also accepts four models that are not in the `RSA`
  package: `"addx2"` and `"addy2"` (the plane with one squared term),
  `"noIA"` (no product term), and `"a5zero"` (equal squared terms, which
  nests `"SRR"` and `"IA"`).

- [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  gains `rsa_model`, which names the fitted model in the `RSA` package’s
  family (`"additive"`, `"IA"`, `"SQD"`, `"onlyx"`, `"null"`, and so
  on). Terms the model sets to zero may be dropped from the Mplus syntax
  or left without a label, and are then returned as 0. They are listed
  in the new `omitted_terms` element. Terms the model estimates must
  still be found, and a nonzero estimate for a term the model sets to
  zero is an error. By default the model is taken from a generated
  workflow and is otherwise `"full"`, as before. Surface parameters that
  the model fixes at zero are shown as 0 in the plot annotation.

### Fitting and comparing the model family

- New
  [`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md)
  writes, runs, and reads a set of latent RSA models from one
  specification, by default the 14 models that
  [`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html) fits. The full
  model is fitted first and its estimates are the start values of the
  others. A model that fails is recorded with the stage at which it
  failed, and the other models are still fitted. The four additional
  models are fitted when they are named in `models`. Input files are
  named in lower case (`rsa_srrr.inp`), because Mplus lower-cases file
  names and would not find `rsa_SRRR.inp` on a case-sensitive file
  system.

- New
  [`compare_rsa_mplus_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_rsa_mplus_models.md)
  compares the fitted models, modelled on
  [`RSA::aictab()`](https://rdrr.io/pkg/RSA/man/aictab.html) and
  [`RSA::compare()`](https://rdrr.io/pkg/RSA/man/compare.html). Models
  are ranked by AICc (counting all free parameters, including the
  measurement model) with Akaike weights and evidence ratios. Each
  restricted model is tested against the full model with the
  Satorra-Bentler scaled likelihood-ratio test from the MLR
  log-likelihoods and scaling correction factors; `reference` names
  another model to test against, and the models nested in it are then
  tested. The table also has the R-square of the latent outcome and
  flags failed or non-converged runs, Mplus errors and warnings, a model
  that fits better than a model it is nested in (with a warning), ridge
  lines of the SRRR, SRR, SRSQD, and SSQD models that lie far outside
  the data or rotated ridge models at `b5 = 0`, and a scaled test with
  `cd <= 0`.

### A shared metric for X and Y

- [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md),
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md),
  and
  [`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md)
  gain `identification` and `center`. `identification = "effects"` puts
  each factor on the metric of its average indicator (effects coding).
  `center = "pooled"` gives the X and Y factors one shared origin, so
  the latent line X = Y is equality on the response scale, as in
  `RSA::RSA(center = "pooled")`. It requires
  `identification = "effects"` or `"invariance"`.

- `identification = "invariance"` is a second way to put the X and Y
  factors on one metric. Matched X and Y indicators share their loadings
  and intercepts, the marker intercepts are fixed at 0, and both factor
  means are freed. Unlike effects coding, this writes no
  `MODEL CONSTRAINT` statements for the measurement model, so it runs
  under `ESTIMATOR = BAYES`, and it is a testable restriction.
  `free_intercepts` frees the intercepts of single indicator pairs
  (partial scalar invariance).

- New
  [`fit_rsa_mplus_invariance()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_invariance.md)
  writes, runs, and compares configural, metric, partial scalar, and
  scalar invariance models of one or more X/Y factor pairs, with
  Satorra-Bentler scaled chi-square difference tests.

### Intervals and the supported region

- New
  [`rsa_mplus_draws()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_draws.md)
  gives Monte Carlo intervals for the surface parameters of a model
  estimated by maximum likelihood. It draws `b1` to `b5` from their
  asymptotic distribution (`TECH3`), imposes the restrictions of
  `rsa_model`, and computes `a1` to `a5` and the ridge or principal-axis
  parameters on each draw.

- New
  [`rsa_mplus_supported_range()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_supported_range.md)
  gives the range of the lines of congruence and incongruence that the
  data support, from the model-implied distribution of the latent X and
  Y factors.

### Minor improvements and fixes

- Mplus is now found from WSL when only the Windows executable
  (`Mplus.exe`) is on the `PATH`.

- [`run_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  explains the failure when Mplus wrote a lower-cased output file next
  to an input file whose name has capital letters, which happens on
  case-sensitive file systems (Linux, WSL). Give `modelout` a lower-case
  file name there.

## rsahelpers 0.1.0

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
  when Mplus constraint p-values are available. With
  `ESTIMATOR = BAYES`, the 95% credibility interval is used for the
  significance decision instead of the one-tailed posterior p-value, so
  a single `*` marks surface parameters whose interval excludes zero.
  The interval bounds are also returned in `new_parameters` as
  `lower_2.5ci` and `upper_2.5ci`.

- rsahelpers no longer declares or installs `franzpak`. The two packages
  previously listed each other in `Suggests` and `Remotes`, which made
  their dependency resolution circular and coupled their CI. The
  dependency is now one-directional: franzpak’s deprecated
  [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  forwards here, and the optional franzpak coefficient table in the
  Mplus vignette is guarded by
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html).

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
- [`select_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/select_spline_congruence.md)
  selects a linear, fixed-LOC, or free one-seam surface using
  residual-bootstrap likelihood-ratio tests and arm-coverage
  diagnostics; absolute-difference, piecewise, and two-seam fits are
  retained as benchmarks or sensitivity analyses.
- New [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for
  `congruence_spline`, enabling
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  and a coefficient table with delta-method standard errors, t
  statistics, and p-values.
- [`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)
  omits invalid Wald tests of seam existence and seam count; its
  remaining seam-location and shape tests are conditional on a supported
  one-seam surface. [`summary()`](https://rdrr.io/r/base/summary.html)
  and
  [`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)
  warn (controllable via `warn_seam`) that delta-method standard errors
  for seam parameters are approximate because the Jacobian is evaluated
  at the non-differentiable seam; the bootstrap is recommended for seam
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
