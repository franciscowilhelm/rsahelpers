# `select_spline_congruence()` review

## Scope

- Review the model-selection sequence and the null residual-bootstrap LRT.
- Expand `vignettes/spline-response-surfaces.qmd` with the inputs, algorithm,
  interpretation, and comparison to delta R-squared, AIC/BIC variants, and
  Wald tests.
- Preserve the user's existing in-progress changes and use targeted validation.

## Status

- 2026-07-22: Started implementation and documentation trace. The worktree was
  already modified in the spline implementation, tests, generated docs, NEWS,
  NAMESPACE, pkgdown config, and the target vignette.

## Decisions and findings

- The seam-existence null is the linear plane (`b3 = 0`); `c0` and `c1` are
  unidentified under that null. The observed statistic is the Gaussian profile
  LRT `n * log(RSS0 / RSS1)`.
- Each null replicate fixes X/Y, samples centered residuals from the fitted
  null, and refits both null and free-seam models. The free fit repeats its
  multistart search. The p-value uses the `(1 + exceedances) / (1 + successes)`
  correction, and failed paired refits are audited against `max_fail`.
- Successful refits are analyzed after dropping failures. The `max_fail` gate
  is useful but cannot rule out selection bias when failure probability is
  related to the bootstrap statistic.
- The residual bootstrap assumes independent, exchangeable/homoskedastic
  errors. It is not robust to heteroskedasticity, clustering, or serial
  dependence; a wild, cluster, or block null bootstrap would be needed.
- Partial R-squared is a monotone transformation of this LRT:
  `T = -n * log(1 - partial_R2)`. It is useful as an effect size, but its
  ordinary F calibration is invalid for seam existence.
- Standard AIC/AICc/BIC parameter-count penalties are descriptive sensitivity
  checks in the non-regular linear-vs-free-seam comparison, not substitute
  p-values. Wald seam-existence tests fail because seam location disappears
  under the null.
- The fixed-LOC comparison is performed only after seam support and arm
  coverage. Its bootstrap is a prudent default because the fitted surface is
  non-differentiable at the seam; a conditional Wald test is only approximate.
- The sequence uses unadjusted alpha and does not claim family-wise error
  control. Sparse or failed structural evidence triggers conservative fallback
  statuses; the two-seam fit is sensitivity-only.
- Expanded the vignette with equations, the full bootstrap algorithm,
  assumptions, failure/coverage gates, Monte Carlo cautions, an alternatives
  table, and code for R-squared/AIC/BIC sensitivity checks.

## Validation

- `air format vignettes/spline-response-surfaces.qmd` completed successfully.
- `quarto render vignettes/spline-response-surfaces.qmd --to html` completed
  successfully after the final edits.
- `Rscript -e "devtools::test(filter = '^spline-regression$')"` passed:
  65 tests, 0 failures, 0 warnings, 0 skips.
- `git diff --check` passed.
