# Test measurement invariance between the X and Y factors

Writes, runs, and reads a sequence of Mplus confirmatory factor models
that test whether matched X and Y indicators measure their factors in
the same way: configural, metric (equal loadings), and scalar (equal
loadings and intercepts). Scalar invariance is what
`identification = "invariance"` in
[`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
assumes when it reads the latent line X = Y as equality on the response
scale.

## Usage

``` r
fit_rsa_mplus_invariance(
  data,
  measurement,
  pairs = NULL,
  free_intercepts = NULL,
  blocks = list(),
  usevariables = NULL,
  modelout,
  overwrite = FALSE,
  run = TRUE,
  Mplus_command = NULL,
  quiet = TRUE
)
```

## Arguments

- data:

  A data frame containing the indicators.

- measurement:

  A named list mapping the factors to their indicators. The two factors
  of a pair need the same number of indicators, in matching order.

- pairs:

  A list of character vectors of length two, each naming an X and a Y
  factor of `measurement`. With several pairs, all factors are estimated
  in one model with free covariances. The default pairs the two factors
  of a two-factor `measurement`.

- free_intercepts:

  Names of indicators whose intercept is not shared within its pair.
  Naming either member frees the pair. When given, a `"partial"` level
  is fitted between `"metric"` and `"scalar"`. At least one matched pair
  of each X/Y factor pair must keep its shared intercept.

- blocks:

  Named list of Mplus input blocks, as in
  [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md).
  The defaults are `ESTIMATOR = MLR;` and `OUTPUT: STDYX;`.

- usevariables:

  Columns of `data` to export. Defaults to the indicators in
  `measurement`.

- modelout:

  Path of an Mplus `.inp` file that names the set. Each level is written
  next to it with its name appended, for example `inv_metric.inp` for
  `modelout = "inv.inp"`, and all levels read the data file `inv.dat`.
  On a case-sensitive file system the file name must be lower case (see
  [`write_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)).

- overwrite:

  If `TRUE`, replace existing input and data files.

- run:

  If `TRUE` (the default), run the models. If `FALSE`, only write the
  input and data files.

- Mplus_command:

  Optional path or command used to invoke Mplus.

- quiet:

  If `TRUE`, suppress status messages from MplusAutomation.

## Value

An object of class `rsa_mplus_invariance`, a list with

- `models`: the named list of workflows, one per level.

- `comparison`: a data frame with one row per level: `level`, `k`,
  `chisq`, `df`, `scaling`, `CFI`, `RMSEA`, `SRMR`, the level it is
  tested `against`, `delta_CFI`, and the scaled difference test
  (`diff_chisq`, `diff_df`, `diff_p`). It is `NULL` when `run = FALSE`.

- `files`: the `input` paths, named by level, and the shared `data`
  path.

## Details

The levels are:

- `"configural"`: the same factor structure, all loadings and intercepts
  free (marker loadings fixed at 1, factor means at 0).

- `"metric"`: matched indicators share their loading.

- `"partial"`: metric, and matched indicators share their intercept
  except those in `free_intercepts`. The Y factor mean is freed.

- `"scalar"`: metric, and all matched indicators share their intercept.

`comparison` tests each level against a less restricted one with the
Satorra-Bentler scaled chi-square difference: metric against configural,
and partial and scalar each against metric. `delta_CFI` is the change in
CFI relative to the same model.

## See also

[`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
invariance <- fit_rsa_mplus_invariance(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3")
  ),
  free_intercepts = "x2",
  modelout = "invariance/inv.inp"
)
invariance$comparison
} # }
```
