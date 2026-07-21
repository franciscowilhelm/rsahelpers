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
  stationary point and principal axes.

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

- quiet:

  If `TRUE`, suppress status messages from MplusAutomation.

- x:

  An `rsa_mplus_workflow` object.

- modelout:

  Path to the Mplus `.inp` file to create.

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
returned by the write, run, and read functions.

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
