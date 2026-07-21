# Extract and plot RSA surfaces from Mplus output

Extract and plot RSA surfaces from Mplus output

## Usage

``` r
RSA_mplus(
  model,
  outcome = NULL,
  pred_x = NULL,
  pred_y = NULL,
  pred_x2 = NULL,
  pred_xy = NULL,
  pred_y2 = NULL,
  b0 = NULL,
  coef_type = c("un", "std", "stdy", "stdyx"),
  new_labels = NULL,
  include_new = TRUE,
  plot = TRUE,
  xlab = NULL,
  ylab = NULL,
  zlab = NULL,
  ...
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

- b0:

  Optional intercept passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). If
  `NULL`, the function looks for `<outcome><-Intercepts>` in the Mplus
  expectation parameters and otherwise falls back to `0`.

- coef_type:

  Which Mplus coefficient table to use. Passed to
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html) for `mplus.model`
  objects.

- new_labels:

  Optional character vector of `NEW` parameter labels to return.
  Matching is case-insensitive.

- include_new:

  If `TRUE`, include `NEW` parameters from `MODEL CONSTRAINT` in the
  returned object.

- plot:

  If `TRUE`, call
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). If
  `FALSE`, only return the extracted coefficients and metadata.

- xlab:

  Optional x-axis label passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). Defaults
  to `pred_x`.

- ylab:

  Optional y-axis label passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). Defaults
  to `pred_y`.

- zlab:

  Optional z-axis label passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). Defaults
  to `outcome`.

- ...:

  Additional arguments passed to
  [`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html).

## Value

A list with the extracted RSA coefficients, selected Mplus parameters,
optional `NEW` parameters, and the plot object when `plot = TRUE`.

## Details

[`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html) accepts the
polynomial coefficients directly (`x`, `y`, `x2`, `xy`, `y2`, and
optionally `b0`). `NEW` parameters from `MODEL CONSTRAINT` are returned
for inspection. When all five generated surface parameters (`CS`, `CC`,
`IS`, `IC`, and `A5`) and their p-values are available, the default
three-dimensional plot annotation uses these Mplus estimates and adds
`*`, `**`, and `***` at p-values no greater than .05, .01, and .001,
respectively. Other plot types and incomplete output retain the
annotation behavior of
[`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html).

## Examples

``` r
if (FALSE) { # \dontrun{
rsa_pos <- RSA_mplus(
  system.file("extdata", "congruence_sim.out", package = "rsahelpers"),
  outcome = "Z",
  pred_x = "X",
  pred_y = "Y",
  pred_x2 = "XS",
  pred_xy = "XY",
  pred_y2 = "YS",
  plot = TRUE
)
} # }
```
