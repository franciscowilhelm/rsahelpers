# Compare the models of the latent RSA family

Ranks fitted models of the RSA family by AICc, as
[`RSA::aictab()`](https://rdrr.io/pkg/RSA/man/aictab.html) does, and
tests each restricted model against the full model with a scaled
likelihood-ratio test, as
[`RSA::compare()`](https://rdrr.io/pkg/RSA/man/compare.html) does. Mplus
reports no chi-square, CFI, RMSEA, or SRMR for models with `XWITH`, so
the comparison rests on the log-likelihoods.

## Usage

``` r
compare_rsa_mplus_models(x, tolerance = 0.01, boundary = 3, reference = "full")
```

## Arguments

- x:

  An `rsa_mplus_family` from
  [`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md),
  or a named list of fitted workflows (see
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)).
  Each workflow's model is taken from the workflow, or otherwise from
  its name.

- tolerance:

  Largest log-likelihood by which a model may exceed a model it is
  nested in before this is flagged. The default absorbs the error of the
  numerical integration.

- boundary:

  Distance, in standard deviations of the latent predictors, beyond
  which a ridge or valley line is flagged as lying outside the data.

- reference:

  The model against which the likelihood-ratio tests are run. The
  default tests every restricted model against the full model. With
  another model, for example `"a5zero"` or `"additive"`, only the models
  nested in it are tested.

## Value

A data frame of class `rsa_mplus_comparison` with one row per model,
ordered by AICc: `model`, `k`, `LL`, `AICc`, `delta_AICc`, `weight`
(Akaike weight), `cum_weight`, `evidence_ratio`, the test against the
`reference` model (`lrt_df`, `lrt_cd`, `lrt_chisq`, `lrt_p`), `R2`,
`converged`, and `flags`. The attribute `messages` holds the Mplus
errors and warnings and the failures of
[`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md),
with columns `model`, `type`, and `message`.

## Details

Model selection follows Schönbrodt's *Testing fit patterns with
polynomial regression models*: AICc is primary and the likelihood-ratio
tests are secondary. AICc is `-2*LL + 2*k + 2*k*(k + 1)/(N - k - 1)`,
with `k` all free parameters, including the measurement model, and `N`
the number of observations. Akaike weights are computed over the
converged models, and the evidence ratio is the weight of the best model
divided by the weight of each model, as in
[`RSA::aictab()`](https://rdrr.io/pkg/RSA/man/aictab.html).

The likelihood-ratio test uses the Satorra-Bentler scaled difference of
the MLR log-likelihoods. With `p0` parameters and scaling correction
factor `c0` in the restricted model and `p1` and `c1` in the full model,
the difference test scaling correction is
`cd = (p0*c0 - p1*c1)/(p0 - p1)` and the test statistic is
`-2*(LL0 - LL1)/cd` on `p1 - p0` degrees of freedom. A `cd` that is not
positive leaves the test undefined; it is flagged and its statistic and
p-value are `NA`. Without a scaling correction factor in the output, as
under `ESTIMATOR = ML`, `c = 1` and the test is the ordinary
likelihood-ratio test.

`flags` collects what needs attention:

- The run failed, did not converge, or Mplus printed errors or warnings.
  Their text is in the `messages` attribute. Models without a
  log-likelihood are listed last and excluded from the ranking.

- A model fits better, by more than `tolerance`, than a model it is
  nested in. Then one of them did not reach its optimum, usually the
  larger one, and a warning is raised. The nesting follows the
  restrictions of
  [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md).

- The ridge or valley line of `"SRRR"`, `"SRR"`, `"SRSQD"`, or `"SSQD"`,
  `Y = S*X + C`, lies more than `boundary` standard deviations from the
  centroid of the latent X and Y, or the rotated models `"SRRR"` and
  `"SRSQD"` have `b5` at 0 (to three decimals). These models are not
  defined at `b5 = 0`. When the surface has no such line, they drift
  towards it and mimic a plane with a far-away line, and `C` and `S`
  cannot be interpreted.

- `cd <= 0` in the likelihood-ratio test.

`R2` is the R-square of the outcome factor as printed by Mplus, and `NA`
when Mplus prints none (it is part of the standardized output).

## See also

[`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md)

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- compare_rsa_mplus_models(family)
comparison
attr(comparison, "messages")
} # }
```
