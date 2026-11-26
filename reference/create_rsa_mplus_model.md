# Build, write, run, and read RSA models for Mplus

These functions provide a composable workflow for response surface
models estimated in Mplus. `create_rsa_mplus_model()` uses tidySEM for
ordinary measurement and structural paths, then adds the `XWITH`
interactions and response-surface constraints required for latent
polynomial models.

## Usage

``` r
create_rsa_mplus_model(
  data,
  measurement,
  structural = NULL,
  roles = c(x = "X", y = "Y", outcome = "Z"),
  model_type = c("latent", "si_lms"),
  reliability = NULL,
  constraints = c("surface", "principal"),
  blocks = list(),
  model_extra = NULL,
  constraint_extra = NULL,
  usevariables = NULL,
  identification = c("marker", "effects", "invariance"),
  center = c("none", "pooled"),
  free_intercepts = NULL,
  rsa_model = "full",
  start = NULL,
  quiet = TRUE
)

write_rsa_mplus_model(
  x,
  modelout,
  dataout = NULL,
  overwrite = FALSE,
  quiet = TRUE
)

run_rsa_mplus_model(
  x,
  Mplus_command = NULL,
  replace_outfile = c("always", "never", "modifiedDate"),
  quiet = TRUE
)

read_rsa_mplus_model(target, output = NULL, quiet = TRUE)

rsa_mplus_workflow(
  data,
  measurement,
  structural = NULL,
  roles = c(x = "X", y = "Y", outcome = "Z"),
  model_type = c("latent", "si_lms"),
  reliability = NULL,
  constraints = c("surface", "principal"),
  blocks = list(),
  model_extra = NULL,
  constraint_extra = NULL,
  usevariables = NULL,
  identification = c("marker", "effects", "invariance"),
  center = c("none", "pooled"),
  free_intercepts = NULL,
  rsa_model = "full",
  start = NULL,
  modelout,
  dataout = NULL,
  overwrite = FALSE,
  run = FALSE,
  Mplus_command = NULL,
  quiet = TRUE
)
```

## Arguments

- data:

  A data frame containing all variables exported to Mplus.

- measurement:

  A named list mapping latent variables to indicators, or a character
  vector of lavaan/tidySEM measurement statements.

- structural:

  Optional structural model. Supply a named list mapping outcomes to
  predictors, or a character vector of lavaan/tidySEM paths. The core X
  and Y paths to the RSA outcome are added when absent.

- roles:

  Named character vector with elements `x`, `y`, and `outcome`
  identifying the latent variables used in the RSA model.

- model_type:

  Either `"latent"` for multi-indicator factors or `"si_lms"` for
  reliability-corrected single-indicator factors.

- reliability:

  For `model_type = "si_lms"`, a named numeric vector of reliability
  estimates for the three latent variables in `roles`.

- constraints:

  Which built-in constraints to include. `"surface"` produces `cs`,
  `cc`, `is`, `ic`, and `a5`; `"principal"` additionally produces the
  stationary point and principal axes for the models that have them (see
  Details).

- blocks:

  Named list of Mplus input blocks. Values replace defaults for the
  corresponding block. `MODEL` and `MODELCONSTRAINT` are reserved; use
  `model_extra` and `constraint_extra` to extend them. The `VARIABLE`
  block may contain additional directives, but `USEVARIABLES` is
  generated from `usevariables` and must not be supplied through
  `blocks`.

- model_extra:

  Additional raw statements appended to the Mplus `MODEL` block.

- constraint_extra:

  Additional raw statements appended to the Mplus `MODEL CONSTRAINT`
  block.

- usevariables:

  Character vector of data columns to export and include in the Mplus
  `USEVARIABLES` statement. Defaults to all columns in `data`.

- identification:

  How the X, Y, and outcome factors are scaled, for
  `model_type = "latent"`. `"marker"` (default) fixes the first loading
  of each factor to 1. `"effects"` uses effects coding (Little, Slegers,
  & Card, 2006): loadings are constrained to average 1, so each factor
  is on the metric of its average indicator. This changes the units of
  b1-b5 slightly but not the fit of the model. `"invariance"` keeps the
  marker loadings and puts the X and Y factors on one metric through
  scalar measurement invariance: the k-th X and Y indicators share their
  loading and their intercept, the two marker intercepts are fixed at 0,
  and both factor means are freed. It requires as many X as Y
  indicators, given in matching order. Unlike `"effects"`, it restricts
  the model, and it writes no `MODEL CONSTRAINT` statements for the
  measurement model, so it also runs under `ESTIMATOR = BAYES`. See
  Details.

- center:

  Origin of the X/Y plane. `"none"` (default) fixes both factor means at
  0, so each factor is centred at its own mean and the line X = Y means
  "both at their average" rather than equality. `"pooled"`, which
  requires `identification = "effects"` or `"invariance"`, subtracts one
  constant (the mean of the X and Y scale means) from every X and Y
  indicator in a `DEFINE` statement and frees the two factor means;
  under effects coding it also constrains the intercepts to sum to zero.
  X = Y is then equality on the shared response scale, the latent
  counterpart of `RSA::RSA(center = "pooled")`. This assumes the X and Y
  items share a response metric. The move is a translation by the latent
  mean difference Delta = E(X) - E(Y): b3-b5, a2, a4, and a5 are
  unchanged, while a3 shifts by `-Delta * a4` and a1 by `-Delta * a5`.
  User `DEFINE` statements are placed after the shift.

- free_intercepts:

  For `identification = "invariance"`, names of indicators whose
  intercept is not shared between X and Y (partial scalar invariance).
  Naming either member of a matched pair frees the pair. The marker pair
  cannot be freed.

- rsa_model:

  Response surface model to generate, named as in the model family of
  the `RSA` package. `"full"` (the default) estimates all five
  polynomial terms. The restricted models are `"SRRR"`, `"SRR"`, `"RR"`,
  `"SRSQD"`, `"SSQD"`, `"SQD"`, `"IA"`, `"additive"`, `"mean"`,
  `"diff"`, `"onlyx"`, `"onlyy"`, `"onlyx2"`, `"onlyy2"`, and `"null"`.
  Four further models are not part of the `RSA` package: `"addx2"`,
  `"addy2"`, `"noIA"`, and `"a5zero"`. See Details.
  [`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md)
  fits several of them from one specification.

- start:

  Optional start values for the polynomial terms, typically the
  estimates of the full model. Supply a named numeric vector with names
  `b1` to `b5` (or `x`, `y`, `x2`, `xy`, and `y2`), a result of
  [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md),
  or a fitted workflow. Only the terms that `rsa_model` estimates freely
  receive a start value.

- quiet:

  If `TRUE`, suppress status messages from MplusAutomation.

- x:

  An `rsa_mplus_workflow` object.

- modelout:

  Path to the Mplus `.inp` file to create. Mplus lower-cases file names,
  so on a case-sensitive file system (Linux, WSL) the file name must be
  lower case for Mplus to find it.

- dataout:

  Path to the Mplus data file. By default, it uses the same stem as
  `modelout` with a `.dat` extension.

- overwrite:

  If `TRUE`, replace existing input and data files.

- Mplus_command:

  Optional path or command used to invoke Mplus.

- replace_outfile:

  Passed to
  [`MplusAutomation::runModels()`](https://michaelhallquist.github.io/MplusAutomation/reference/runModels.html).

- target:

  An `rsa_mplus_workflow`, Mplus `.inp` or `.out` path, an
  `mplus.model`, or an `mplusObject` containing parsed results.

- output:

  Optional Mplus `.out` path to read into an existing workflow. This is
  useful for attaching a cached or relocated output while retaining the
  workflow's role metadata. It can only be used when `target` is an
  `rsa_mplus_workflow`.

- run:

  If `TRUE`, run and read the generated Mplus model after writing its
  files. The default only writes reproducible input and data files.

## Value

An object of class `rsa_mplus_workflow`. The object is updated and
returned by the write, run, and read functions. Under pooled centering,
`spec$centering` holds the subtracted constant and the two scale means
it was computed from. Under invariance identification, `spec$invariance`
holds the matched indicators and the freed pairs. `spec$rsa_model` names
the generated model and `spec$start` holds the start values that were
used.

## Details

`rsa_model` writes the restrictions of the second-order models in
[`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html) (Schönbrodt,
*Testing fit patterns with polynomial regression models*) on the
regression weights `b1` to `b5`. Terms set to zero are fixed with `@0`.
The other restrictions are placed in `MODEL CONSTRAINT`, together with
the `NEW` parameters that describe the restricted surface:

|  |  |  |
|----|----|----|
| `rsa_model` | Restrictions | Model-specific `NEW` parameters |
| `"full"` | none | `x0`, `y0`, `p10`, `p11`, `p20`, `p21` (with `"principal"`) |
| `"SRRR"` | `b3 = b4*b4/(4*b5)` | `c`, `s`, `bm`, `a4r` |
| `"SRR"` | `b3 = b5`, `b4 = -2*b5` | `c`, `bm` |
| `"RR"` | `b1 = b2`, `b3 = b5`, `b4 = -2*b5` | `bm` |
| `"SRSQD"` | `b1 = b2*b4/(2*b5)`, `b3 = b4*b4/(4*b5)` | `c`, `s`, `a4r` |
| `"SSQD"` | `b1 = -b2`, `b3 = b5`, `b4 = -2*b5` | `c` |
| `"SQD"` | `b1 = b2 = 0`, `b3 = b5`, `b4 = -2*b5` | none |
| `"IA"` | `b3 = b5 = 0` | `x0`, `y0` (with `"principal"`) |
| `"additive"` | `b3 = b4 = b5 = 0` | none |
| `"mean"`, `"diff"` | `b3 = b4 = b5 = 0`, and `b1 = b2` or `b1 = -b2` | none |
| `"onlyx2"`, `"onlyy2"` | `b2 = b4 = b5 = 0` or `b1 = b3 = b4 = 0` | none |
| `"onlyx"`, `"onlyy"` | all terms but `b1` or `b2` are 0 | none |
| `"null"` | all five terms are 0 | none |
| `"addx2"`, `"addy2"` | `b4 = b5 = 0` or `b3 = b4 = 0` | none |
| `"noIA"` | `b4 = 0` | none |
| `"a5zero"` | `b3 = b5` | none |

The last four models are additions to the family of
[`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html). `"addx2"` and
`"addy2"` add one squared term to the plane; unlike `"onlyx2"` and
`"onlyy2"`, they keep both linear terms. `"noIA"` keeps both squared
terms and drops the product. `"a5zero"` sets the two squared terms
equal, so that `a5 = 0`; it is the smallest model that nests both
`"SRR"` and `"IA"`, which makes it a reference model for testing either
of them (see `reference` in
[`compare_rsa_mplus_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_rsa_mplus_models.md)).

The ridge parameters correspond to those of
[`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html): `c` is the lateral
shift C of the ridge in units of X, `s` the scaling S of X that rotates
the ridge, `bm` the effect of the mean level along the ridge
(`meaneffect`), and `a4r` the curvature across the rotated ridge
(`a4.rescaled`, a4'). Each is written in its simplest closed form in
`b1` to `b5`; for example, a4' is `4*b5`. The rotated models `"SRRR"`
and `"SRSQD"` express `b3` (and `b1`) as functions of the other weights.
One run therefore covers both ridges and troughs, where
[`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html) fits an
inequality-constrained model twice. These models are not defined at
`b5 = 0`: a flat-ridge model that does not fit the data can drift
towards `b5 = 0` with a large `c`.

Surface parameters (`cs`, `cc`, `is`, `ic`, `a5`) that a model fixes at
zero are not generated, and the others omit the terms fixed with `@0`,
because Mplus cannot refer to such a term in `MODEL CONSTRAINT`. With
`"principal"`, the full model gets its stationary point and principal
axes and `"IA"` its stationary point; the principal axes of `"IA"` have
slopes of 1 and -1 by construction, and the other models have no
stationary point. Restricted models keep the `XWITH` products, with
fixed weights where needed, so all models are nested in the full model
and their log-likelihoods can be compared. The `"weak"` and `"strong"`
fit patterns of [`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html) are
not generated.

With `identification = "invariance"`, the latent line X = Y is equality
on the response scale because matched items measure X and Y with the
same loadings and intercepts. This is an assumption about the items: it
holds only if the X and Y items are parallel in wording and response
format, and it can be tested with
[`fit_rsa_mplus_invariance()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_invariance.md).
Each factor mean is then the mean of its marker item. With
`center = "none"`, the origin X = Y = 0 is the value 0 of the response
scale, which usually lies outside the data; `center = "pooled"` moves it
to the mean of the two scale means. To centre elsewhere, for example at
the scale midpoint, use `center = "none"` and subtract that value from
the X and Y indicators in `blocks$DEFINE`. The outcome factor keeps its
marker identification.

## Examples

``` r
if (FALSE) { # \dontrun{
latent_model <- create_rsa_mplus_model(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  structural = list(Z = c("X", "Y"))
)
} # }
```
