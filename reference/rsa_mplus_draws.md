# Monte Carlo intervals for response surface parameters

Draws the polynomial weights `b1` to `b5` of a fitted Mplus response
surface from their asymptotic sampling distribution (a multivariate
normal with the estimates as means and the `TECH3` covariance matrix)
and computes the surface parameters on each draw. This gives intervals
for parameters whose delta-method standard errors are poor, such as
ratios, without refitting the model, and draws on which any further
function of the surface can be evaluated.

## Usage

``` r
rsa_mplus_draws(
  model,
  outcome = NULL,
  pred_x = NULL,
  pred_y = NULL,
  pred_x2 = NULL,
  pred_xy = NULL,
  pred_y2 = NULL,
  rsa_model = NULL,
  n = 10000L,
  probs = c(0.025, 0.975),
  seed = NULL
)
```

## Arguments

- model:

  Path to an Mplus `.out` file, an object returned by
  [`MplusAutomation::readModels()`](https://michaelhallquist.github.io/MplusAutomation/reference/readModels.html),
  an `mplusObject` containing results, or an
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  object.

- outcome:

  Dependent variable label from the Mplus regression table. Inferred for
  generated RSA-Mplus workflows.

- pred_x:

  Label of the linear X predictor in the Mplus output. Inferred for
  generated RSA-Mplus workflows.

- pred_y:

  Label of the linear Y predictor in the Mplus output. Inferred for
  generated RSA-Mplus workflows.

- pred_x2:

  Label of the squared X term in the Mplus output. Inferred for
  generated RSA-Mplus workflows.

- pred_xy:

  Label of the XY interaction term in the Mplus output. Inferred for
  generated RSA-Mplus workflows.

- pred_y2:

  Label of the squared Y term in the Mplus output. Inferred for
  generated RSA-Mplus workflows.

- rsa_model:

  Name of the fitted model in the `RSA` package's model family. Terms
  that this model sets to zero may be missing from the Mplus output, and
  their labels are optional. One of `"full"`, `"additive"`, `"IA"`,
  `"SQD"`, `"SRRR"`, `"SRR"`, `"RR"`, `"SSQD"`, `"SRSQD"`, `"null"`,
  `"onlyx2"`, `"onlyy2"`, `"onlyx"`, `"onlyy"`, `"mean"`, `"diff"`,
  `"weak"`, or `"strong"`. If `NULL` (the default), the model that
  [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  generated for a workflow, and otherwise `"full"`. Also passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html) as
  `model`. See Details.

- n:

  Number of draws.

- probs:

  Lower and upper quantile of the reported intervals.

- seed:

  Optional seed, set locally for the draws.

## Value

An object of class `rsa_mplus_draws`, a list with

- `draws`: a data frame with one row per draw.

- `summary`: a data frame with the `parameter`, its estimate `est`
  computed from the point estimates, the standard deviation `se` of the
  draws, and the `lower` and `upper` quantiles.

- `rsa_model`, `n`, and `probs`.

## Details

The model must have been estimated by maximum likelihood with `TECH1`
and `TECH3` in the `OUTPUT` command, for example
`blocks = list(OUTPUT = "STDYX CINTERVAL TECH1 TECH3;")` in
[`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md).
Bayesian models have no `TECH3`; use their posterior draws
(`BPARAMETERS`) instead.

Terms that share a parameter in the model, for example through a shared
label, share its draws. The restrictions of `rsa_model` are imposed on
the draws and on the estimates: terms the model sets to zero are 0, and
terms it defines from others (such as `b4 = -2*b5` in `"SRR"`) are
computed from them. This also removes the rounding of the printed
estimates, which do not satisfy the restrictions exactly.

Besides `b1` to `b5`, the draws hold the surface parameters `a1` to `a5`
(`cs`, `cc`, `is`, `ic`, and `a5` in the generated `MODEL CONSTRAINT`),
and the parameters specific to `rsa_model` in
[`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md):
the ridge parameters `c`, `s`, `bm`, and `a4r`, or the stationary point
and principal axes. Ratios can be undefined or extreme in single draws,
so their intervals are quantiles and the standard deviation in `se` may
be large.

Mplus prints estimates and covariances with limited precision. The draws
are therefore centred on estimates rounded to three decimals.

## See also

[`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md),
[`rsa_mplus_supported_range()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_supported_range.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fitted <- rsa_mplus_workflow(
  data = dat,
  measurement = measurement,
  rsa_model = "SRR",
  blocks = list(OUTPUT = "TECH1 TECH3;"),
  modelout = "rsa_srr.inp",
  run = TRUE
)
intervals <- rsa_mplus_draws(fitted, seed = 1)
intervals$summary
} # }
```
