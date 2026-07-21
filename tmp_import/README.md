# tmp_export — RSA helper functions for reuse

Copied and lightly repackaged from the `rsasim` project `lib/` for dropping
into a general RSA helper R package. Every external call is now
namespace-qualified (`pkg::fn`) so the files work as package sources with only
`Imports:` entries in `DESCRIPTION` — no `library()` calls, no `source()`
chains. Roxygen headers with `@export`/`@param`/`@return` are included; run
`devtools::document()` after copying in.

## What's here

| File | Provides | Origin |
|------|----------|--------|
| `surface-math.R` | `ld_to_b`, `b_to_ld`, `rsa_a_parameters`, `rsa_p_parameters`, `lateral_shift_quantity` | `lib/surface_math.R` (subset) |
| `classify-unconstrained.R` | Yao & Ma type **from the full model**: `classify_yao_ma` (tolerance), `classify_yao_ma_rsa` (significance) + `fit_quadratic_lm`, `estimates_from_b`, `estimates_from_lm`, `yao_ma_type_from_flags`, `fit_rsa`, `tidy_rsa_params`, `flip_valence_b` | `lib/estimate.R` |
| `classify-constrained.R` | Yao & Ma type **from constrained models**: `surface_model_specs`, `fit_constrained_ls`, `compare_surface_models_ls`/`_lavaan`, `test_surface_theory`, and the new `classify_yao_ma_constrained` | `lib/model_comparison.R` (+ 1 new fn) |
| `predictor-coverage.R` | Cell-proportions coverage + heatmap: `predictor_coverage`, `plot_predictor_coverage`, `coverage_cells_long` + the Likert generator (`simulate_pe_predictors`, `likert_thresholds`, `likert_margin_presets`) | `lib/simulate.R` + report 07 inline ggplot |

## The two classification routes (as requested)

Both take centered predictors `P`, `E` and an `outcome` column and return a
Yao & Ma (2023) Table-4 type.

- **Non-constrained (one full fit, read the a-/p-parameters):**
  ```r
  est <- estimates_from_lm(fit_quadratic_lm(data))
  classify_yao_ma(est, valence = "positive")          # tolerance tier
  classify_yao_ma_rsa(fit_rsa(data), boot_ci = FALSE)  # significance tier
  ```
- **Constrained (fit the family of type-constrained models, select one):**
  ```r
  classify_yao_ma_constrained(data, method = "ls")     # fast, dependency-light
  classify_yao_ma_constrained(data, method = "lavaan") # robust MLR fits
  ```
  `classify_yao_ma_constrained` is **new** (not in the original `lib/`): it
  composes `compare_surface_models_*` over all 16 type specs and returns the
  declared `yao_ma_type` of the AICc-best tenable model. The original
  `test_surface_theory` (tests one named theory) is also included.

## The coverage heatmap (report 07)

```r
cov <- predictor_coverage(data)          # data needs raw P_raw / E_raw
plot_predictor_coverage(cov)             # single heatmap
plot_predictor_coverage(list(a = cov1, b = cov2))  # faceted by name
```
`plot_predictor_coverage` packages the inline ggplot from
`qmd/07-predictor-complications.qmd`. It accepts a `predictor_coverage()`
result, a bare `cell_props` table, or a named list (faceted). The report's
own facet grid (`margins ~ rho`) can be reproduced by binding several
`coverage_cells_long()` frames and adding your own `facet_grid()`.

## Dependencies

Add to `DESCRIPTION` under `Imports:`:

- **Core / always:** `tibble`, `dplyr`, `purrr`, `stats`, `MASS`
- **Coverage heatmap:** `ggplot2` (and `rlang` for the `.data` pronoun used in
  `plot_predictor_coverage`)
- **Constrained lavaan path & the RSA-based classifier:** `RSA`, `lavaan`
  (only needed for `method = "lavaan"`, `classify_yao_ma_rsa`, `fit_rsa`; the
  `ls` path and the tolerance-tier classifier need neither — consider
  `Suggests:` + `requireNamespace()` guards if you want them optional)

## Notes / things to adjust on integration

- `%||%` is defined locally in `classify-constrained.R`. R ≥ 4.4 has a base
  `%||%`; `rlang`/`purrr` also export one. Delete the local copy if your
  package already imports one, to avoid a masking NOTE.
- Functions prefixed `.` (e.g. `.b_constraint_lexicon`, `.null_space`,
  `.lavaan_aicc`) are internal helpers — leave them unexported.
- `plot_predictor_coverage` uses `ggplot2::aes(.data$E, ...)`; either
  `@importFrom rlang .data` or add `utils::globalVariables(".data")` /
  `importFrom(rlang, .data)` to NAMESPACE.
- The tolerance defaults in `classify_yao_ma` (`eps_a = 0.10`, etc.) are
  calibrated to the **centered 1-5 Likert scale** (predictors in `[-2, 2]`).
  Rescale the tolerances if your predictors live on a different range.
- `fit_rsa` assumes predictors are **pre-centered** (`center = "none"`) and
  the data is clean (`out.rm = FALSE`). Change these if you feed raw scores.
- Sign/valence convention throughout: deficiency = P > E, excess = E > P;
  negative-valence outcomes are reverse-coded before classification.

## Notes for our integration

- we dont use P / E usually but X / Y. It matters for classifying excess and deficiency, but is not generalizable naturally to all applications of RSA functions.
- the fucntions were made with manifest scores in mind, implemented via lm or lavaan. instead we use mplus and latent models, so we need to consider their differences.
- centering / pre-centering isnt relevant to SEM for now.
- tolerance-based classification into Yao and Ma: we need to think about this.
- classify-constrained: not possible right now, defer, as very computationally expensive if requiring different model constraitns and model (re)fits of a latent variable model with latent interactions.
- surface math differs from our nomenclature, stick to ours for now.