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
  rsa_model = NULL,
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
`coefficients` always contains all five polynomial terms and `b0`.
`regression_parameters` contains only the terms found in the Mplus
output, and `omitted_terms` names the terms that were absent and set to
zero.

## Details

[`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html) accepts the
polynomial coefficients directly (`x`, `y`, `x2`, `xy`, `y2`, and
optionally `b0`).

The `RSA` package compares a family of restricted second-order surfaces
against the full polynomial. In Mplus a restriction is written in one of
three ways: a term fixed at zero (`Z ON XS@0;`), shared labels with
`MODEL CONSTRAINT` (as in the ridge models), or a term dropped from the
model, which leaves no row in the output. Terms fixed at zero and
constrained terms are read as usual. Dropped terms need `rsa_model`,
which names the model and thereby the terms it sets to zero:

|  |  |
|----|----|
| `rsa_model` | Terms set to zero |
| `"full"`, `"SRRR"`, `"SRR"`, `"RR"`, `"SSQD"`, `"SRSQD"`, `"weak"`, `"strong"`, `"a5zero"` | none |
| `"noIA"` | `xy` |
| `"addx2"` | `xy`, `y2` |
| `"addy2"` | `x2`, `xy` |
| `"SQD"` | `x`, `y` |
| `"IA"` | `x2`, `y2` |
| `"additive"`, `"mean"`, `"diff"` | `x2`, `xy`, `y2` |
| `"onlyx2"` | `y`, `xy`, `y2` |
| `"onlyy2"` | `x`, `x2`, `xy` |
| `"onlyx"` | `y`, `x2`, `xy`, `y2` |
| `"onlyy"` | `x`, `x2`, `xy`, `y2` |
| `"null"` | all five |

A term set to zero is taken as 0 when its label is `NULL` or its row is
absent from the output. If the row is present, its estimate must be 0; a
nonzero estimate means the output is not the named model and is an
error. The remaining terms must be found, so a mistyped label still
fails. One set of labels can therefore be reused across all variants.
Equality constraints such as `b1 = b2` are not checked. The cubic and
absolute-difference models of `RSA` are not polynomial surfaces of
second order and are not supported.

`NEW` parameters from `MODEL CONSTRAINT` are returned for inspection.
When the generated surface parameters (`CS`, `CC`, `IS`, `IC`, and `A5`)
are available, the default three-dimensional plot annotation uses these
Mplus estimates and adds significance markers. Surface parameters that
`rsa_model` fixes at zero, such as `CC` in the ridge models, may be
missing and are shown as 0 without a marker. For frequentist estimators
the markers are `*`, `**`, and `***` at p-values no greater than .05,
.01, and .001, respectively. For `ESTIMATOR = BAYES` Mplus reports a 95%
credibility interval next to a one-tailed posterior p-value; the
credibility interval is the preferred decision rule, so a single `*`
marks parameters whose interval excludes zero and the posterior p-value
is ignored. The credibility interval bounds are also added to
`new_parameters` as `lower_2.5ci` and `upper_2.5ci`. Other plot types
and incomplete output retain the annotation behavior of
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

# An additive model without quadratic terms in the Mplus syntax
additive <- RSA_mplus(
  "additive.out",
  outcome = "Z",
  pred_x = "X",
  pred_y = "Y",
  rsa_model = "additive"
)
} # }
```
