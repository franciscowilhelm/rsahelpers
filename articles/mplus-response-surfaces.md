# Response Surface Models with Mplus

## 1 Overview

`rsahelpers` provides an end-to-end workflow for polynomial response
surface models estimated in Mplus. It can:

1.  generate conventional measurement and structural syntax with
    tidySEM;
2.  add the latent `XWITH` interactions required for quadratic terms;
3.  generate response-surface and principal-axis constraints;
4.  export an Mplus input file and its data;
5.  run and read the model through MplusAutomation; and
6.  extract and plot the fitted surface with
    [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md).

The workflow supports two measurement strategies:

- **latent LMS**, with multiple indicators for each RSA variable; and
- **SI-LMS**, with one prepared score per construct and a fixed
  measurement error based on a supplied reliability estimate.

> **Note:** Mplus is needed only to estimate a generated input file.
> Model construction, syntax inspection, and reading an existing `.out`
> file do not invoke the Mplus executable.

## 2 Model and parameter conventions

The package uses the common RSA convention in which `Z` is the outcome
and `X` and `Y` are commensurate predictors:

``` math
Z = b_0 + b_1X + b_2Y + b_3X^2 + b_4XY + b_5Y^2 + e.
```

The generated latent interactions are named `XS`, `XY`, and `YS`. The
five structural paths are labeled `b1` through `b5`, which lets Mplus
calculate the surface parameters in `MODEL CONSTRAINT`.

| Mplus label | RSA parameter | Formula | Interpretation |
|----|---:|----|----|
| `cs` | $`a_1`$ | $`b_1 + b_2`$ | Slope along the line of congruence |
| `cc` | $`a_2`$ | $`b_3 + b_4 + b_5`$ | Curvature along the line of congruence |
| `is` | $`a_3`$ | $`b_1 - b_2`$ | Slope along the line of incongruence |
| `ic` | $`a_4`$ | $`b_3 - b_4 + b_5`$ | Curvature along the line of incongruence |
| `a5` | $`a_5`$ | $`b_3 - b_5`$ | Difference between predictor curvatures |

The default constraints also calculate the stationary point and the
intercepts and slopes of both principal axes: `x0`, `y0`, `p10`, `p11`,
`p20`, and `p21`.

## 3 Checking predictor coverage

Response surfaces should not be interpreted in combinations of X and Y
that are unsupported by the data. For scale scores derived from ordinal
items,
[`predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/predictor_coverage.md)
uses common equal-width bins across the pooled X/Y range. This preserves
the geometry of `X = Y` without creating a cell for every possible
decimal score.

``` r

scale_scores <- data.frame(
  X = rowMeans(dat[c("x1", "x2", "x3")]),
  Y = rowMeans(dat[c("y1", "y2", "y3")])
)
coverage <- predictor_coverage(
  scale_scores,
  X,
  Y,
  breaks = 10,
  discrepancy = c(0.5, 1)
)
coverage$summary
```

      x y   n n_missing correlation x_greater_than_y y_greater_than_x ties
    1 X Y 300         0   0.3954739              0.5              0.5    0
      empty_cells sparse_cells min_cell_count min_count
    1          44           75              0         5

``` r

coverage$discrepancy
```

      cutoff        direction count proportion
    1    0.5 x_greater_than_y    92  0.3066667
    2    0.5 y_greater_than_x   102  0.3400000
    3    1.0 x_greater_than_y    58  0.1933333
    4    1.0 y_greater_than_x    67  0.2233333

``` r

plot_predictor_coverage(coverage)
```

![Heatmap of joint X and Y scale-score coverage with sparse bins marked
by red
crosses.](mplus-response-surfaces_files/figure-html/fig-predictor-coverage-1.png)

Figure 1: Joint coverage of the two scale-score predictors. Red crosses
mark bins containing fewer than five observations.

The summary uses neutral X/Y direction labels. Substantive labels such
as deficiency or excess depend on the constructs assigned to X and Y and
should be applied during interpretation rather than inferred by the
function.

## 4 Multi-indicator latent LMS

### 4.1 List input

The most concise specification is a named list mapping each factor to
its indicators. `roles` defaults to
`c(x = "X", y = "Y", outcome = "Z")`, so it does not need to be supplied
when those names are used.

``` r

latent_model <- create_rsa_mplus_model(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  structural = list(Z = c("X", "Y")),
  usevariables = analysis_variables
)

latent_model
```

    <rsa_mplus_workflow>
    Status: created
    Model type: latent 

The core `Z ON X` and `Z ON Y` paths are added when they are absent, so
the `structural` argument can be omitted for the basic RSA model. A
structural list is useful when the outcome has additional predictors:

``` r

latent_with_age <- create_rsa_mplus_model(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  structural = list(Z = c("X", "Y", "age")),
  usevariables = c(analysis_variables, "age")
)
```

### 4.2 Character syntax

Measurement and structural models can instead be supplied as
lavaan/tidySEM-style character vectors. This form is useful for labeled
loadings, covariances, intercepts, and larger structural models.

``` r

character_model <- create_rsa_mplus_model(
  data = dat,
  measurement = c(
    "X =~ x1 + x2 + x3",
    "Y =~ y1 + y2 + y3",
    "Z =~ z1 + z2 + z3"
  ),
  structural = c(
    "Z ~ X + Y + age",
    "X ~~ Y"
  ),
  usevariables = c(analysis_variables, "age")
)
```

### 4.3 Inspecting generated syntax

The returned `rsa_mplus_workflow` stores both the tidySEM representation
and the assembled `mplusObject`.

``` r

cat(latent_model$mplus$VARIABLE)
```

    USEVARIABLES =
      x1 x2 x3 y1 y2 y3 z1 z2 z3
    ;

``` r

model_lines <- strsplit(latent_model$mplus$MODEL, "\n", fixed = TRUE)[[1]]
cat(grep(" BY | ON |XWITH", model_lines, value = TRUE), sep = "\n")
```

    X BY x1@1;
    X BY x2;
    X BY x3;
    Y BY y1@1;
    Y BY y2;
    Y BY y3;
    Z BY z1@1;
    Z BY z2;
    Z BY z3;
    Z ON X (b1);
    Z ON Y (b2);
    XY | X XWITH Y;
    XS | X XWITH X;
    YS | Y XWITH Y;
    Z ON XS (b3);
    Z ON XY (b4);
    Z ON YS (b5);

``` r

cat(latent_model$mplus$MODELCONSTRAINT)
```

    NEW(cs cc is ic a5 x0 y0 p10 p11 p20 p21);
    cs = b1 + b2;
    cc = b3 + b4 + b5;
    is = b1 - b2;
    ic = b3 - b4 + b5;
    a5 = b3 - b5;
    x0 = (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4*b4);
    y0 = (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4*b4);
    p11 = (b5 - b3 + SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;
    p21 = (b5 - b3 - SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;
    p10 = y0 - p11*x0;
    p20 = y0 - p21*x0;

Use `tidySEM::syntax(latent_model$tidysem)` when the underlying
parameter table is more useful than the rendered Mplus statements.

## 5 Single-indicator LMS

SI-LMS uses one prepared score for each construct. If $`r`$ is the
supplied reliability and $`s^2`$ is the sample variance of the score,
its residual variance is fixed to

``` math
\theta = (1-r)s^2.
```

The score columns must be numeric, and reliabilities must be greater
than zero and no greater than one.

``` r

scores <- data.frame(
  xmean = rowMeans(dat[c("x1", "x2", "x3")]),
  ymean = rowMeans(dat[c("y1", "y2", "y3")]),
  zmean = rowMeans(dat[c("z1", "z2", "z3")])
)

si_lms_model <- create_rsa_mplus_model(
  data = scores,
  measurement = list(
    X = "xmean",
    Y = "ymean",
    Z = "zmean"
  ),
  model_type = "si_lms",
  reliability = c(X = 0.82, Y = 0.79, Z = 0.88)
)

si_lines <- strsplit(si_lms_model$mplus$MODEL, "\n", fixed = TRUE)[[1]]
cat(grep(" BY |mean@", si_lines, value = TRUE), sep = "\n")
```

    X BY xmean@1;
    Y BY ymean@1;
    Z BY zmean@1;
    xmean@0.164651744826098;
    ymean@0.280183186480885;
    zmean@0.0924327057256224;

`rsahelpers` does not calculate the scores or their reliability
estimates. This keeps the scoring decisions explicit and allows alpha,
omega, or another appropriate reliability estimate to be supplied.

## 6 Customizing Mplus blocks

`blocks` replaces generated defaults by Mplus section name. The default
`ANALYSIS` block uses `TYPE = RANDOM` and integration, and the default
`OUTPUT` block requests `STDYX`, confidence intervals, and `TECH1`.

``` r

custom_model <- create_rsa_mplus_model(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  usevariables = c(analysis_variables, "age"),
  blocks = list(
    TITLE = "Latent RSA with four processors;",
    ANALYSIS = paste(
      "TYPE = RANDOM;",
      "ALGORITHM = INTEGRATION;",
      "PROCESSORS = 4;",
      sep = "\n"
    ),
    OUTPUT = "STDYX CINTERVAL TECH1 TECH4;"
  ),
  model_extra = "Z WITH age;",
  constraint_extra = "NEW(total_loc); total_loc = cs + cc;"
)
```

`TYPE = RANDOM` remains mandatory because the generated model contains
`XWITH`. `MODEL` and `MODEL CONSTRAINT` cannot be replaced through
`blocks`; use `model_extra` and `constraint_extra` to append statements
without removing the required RSA foundation.

To omit the stationary point and principal-axis formulas, request only
the surface constraints:

``` r

surface_only <- create_rsa_mplus_model(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  constraints = "surface",
  usevariables = analysis_variables
)
```

This is useful when a principal-axis expression is undefined or
unstable, such as when $`b_4`$ or $`4b_3b_5-b_4^2`$ is near zero.

## 7 Writing input and data files

The composable writer creates the `.inp` and `.dat` files without
running Mplus. Existing targets are protected unless `overwrite = TRUE`.

``` r

written <- write_rsa_mplus_model(
  latent_model,
  modelout = "models/latent-rsa.inp"
)

written$files
```

The convenience workflow performs model construction and file writing in
one call. Its default `run = FALSE` is deliberate.

``` r

written <- rsa_mplus_workflow(
  data = dat,
  measurement = list(
    X = c("x1", "x2", "x3"),
    Y = c("y1", "y2", "y3"),
    Z = c("z1", "z2", "z3")
  ),
  usevariables = analysis_variables,
  modelout = "models/latent-rsa.inp"
)
```

MplusAutomation writes every exported column under `NAMES`. rsahelpers
also writes an explicit `USEVARIABLES` statement containing the columns
selected by `usevariables`. This distinction matters: when
`USEVARIABLES` is absent, Mplus analyzes every variable under `NAMES`,
including columns that were present only for data management. The basic
example therefore selects the nine indicators and excludes `age`; the
covariate examples include it intentionally.

By default `usevariables` contains all columns in `data`. Selected names
must be valid Mplus identifiers and unique in their first eight
characters so parsed output remains unambiguous. Other `VARIABLE`
directives can be supplied through `blocks$VARIABLE`, but `USEVARIABLES`
itself is controlled by the dedicated argument.

## 8 Running and reading Mplus

Run a written workflow with
[`run_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md):

``` r

fitted <- run_rsa_mplus_model(written)
```

The runner first checks an explicitly supplied executable and otherwise
looks for `mplus` or `Mplus` on the command path. An explicit path is
useful when the application is installed but not discoverable by the
shell:

``` r

fitted <- run_rsa_mplus_model(
  written,
  Mplus_command = "/Applications/Mplus/mplus"
)
```

Set `run = TRUE` in
[`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
to build, write, run, and parse in a single call. Model files and Mplus
output are retained if estimation fails.

An existing `.inp` or `.out` file can be read independently:

``` r

fitted <- read_rsa_mplus_model("models/latent-rsa.out")
```

This vignette reads a packaged output cache generated from the corrected
model above. Supplying it as `output` attaches the estimates to
`latent_model` and retains the role metadata. Thus the vignette
evaluates result extraction and plotting without requiring Mplus during
package installation or rendering.

``` r

cached_output <- system.file(
  "extdata",
  "mplus-response-surfaces.out",
  package = "rsahelpers"
)
if (!nzchar(cached_output)) {
  cache_candidates <- c(
    file.path("inst", "extdata", "mplus-response-surfaces.out"),
    file.path("..", "inst", "extdata", "mplus-response-surfaces.out")
  )
  cached_output <- cache_candidates[file.exists(cache_candidates)][[1]]
}

fitted <- read_rsa_mplus_model(latent_model, output = cached_output)
fitted
```

    <rsa_mplus_workflow>
    Status: read
    Model type: latent
    Output: /home/runner/work/_temp/Library/rsahelpers/extdata/mplus-response-surfaces.out 

Standalone files do not contain the role metadata created by
[`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md).
Consequently, plotting an arbitrary output requires explicit variable
labels.

## 9 Extracting and plotting the response surface

[`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
accepts a fitted workflow directly. It infers `X`, `Y`, `Z`, and the
three interaction labels from workflow metadata.

``` r

surface <- RSA_mplus(fitted, plot = FALSE)
surface$coefficients
```

         x      y     x2     xy     y2     b0
     0.287  0.260 -0.180  0.220 -0.086  0.000 

``` r

surface$new_parameters[c("Label", "est", "se", "pval")]
```

       Label    est     se  pval
    1     CS  0.548  0.052 0.000
    2     CC -0.046  0.033 0.165
    3     IS  0.027  0.084 0.747
    4     IC -0.486  0.110 0.000
    5     A5 -0.095  0.046 0.039
    6     X0  8.011  8.955 0.371
    7     Y0 11.816 13.514 0.382
    8    P10 -0.348  0.306 0.255
    9    P11  1.519  0.307 0.000
    10   P20 17.091 19.379 0.378
    11   P21 -0.659  0.133 0.000

## 10 Classifying the surface

[`classify_yao_ma()`](https://franciscowilhelm.github.io/rsahelpers/reference/classify_yao_ma.md)
reads the surface and principal-axis constraints directly from the
fitted Mplus model. Outcome valence is required because positively
valenced outcomes use the first principal axis (`P10`, `P11`), whereas
negatively valenced outcomes use the second (`P20`, `P21`).

``` r

classification <- classify_yao_ma(fitted, valence = "positive")
classification
```

    <yao_ma_classification>
    Method: significance
    Valence: positive
    Status: classified
    Type: 2 - Exact correspondence & LLE 

``` r

classification$decisions
```

           condition parameter estimate conf.low conf.high null margin
    1             CS        CS    0.548    0.447     0.649    0     NA
    2             CC        CC   -0.046   -0.111     0.019    0     NA
    3             IS        IS    0.027   -0.137     0.191    0     NA
    4             IC        IC   -0.486   -0.701    -0.271    0     NA
    5 axis_intercept       P10   -0.348   -0.948     0.251    0     NA
    6     axis_slope       P11    1.519    0.918     2.119    1     NA
              state direction
    1     different     above
    2 not_different  overlaps
    3 not_different  overlaps
    4     different     below
    5 not_different  overlaps
    6 not_different  overlaps

By default, an interval excluding its null value indicates a difference
and an interval containing the null is treated as absence of a
difference.

Alternatively, a more conservative equivalence analysis can be done by
supplying margins of equivalence. Supply separate margins for `CS`,
`CC`, `IS`, `IC`, the principal-axis intercept, and its slope. Margins
do not have defaults because their units depend on the model’s working
scale. Intervals that neither fit fully inside nor fully outside the
requested equivalence region produce an indeterminate classification.

``` r

classify_yao_ma(
  fitted,
  valence = "positive",
  equivalence = c(
    CS = 0.10,
    CC = 0.10,
    IS = 0.10,
    IC = 0.10,
    axis_intercept = 0.15,
    axis_slope = 0.15
  )
)
```

Principal-axis constraints are nonlinear functions of the polynomial
coefficients. Their Mplus intervals use delta-method standard errors;
unstable or unavailable axes are therefore reported as a partial result
rather than being assigned a subtype.

### 10.1 Publication-style coefficient table

The optional [franzpak](https://github.com/franciscowilhelm/franzpak)
package can format both ordinary Mplus parameters and `MODEL CONSTRAINT`
results. Its function expects the parsed `mplus.model`, which is stored
in `fitted$results`. franzpak is not on CRAN and is not required by
rsahelpers, so the table below is only rendered when it is installed.

``` r

if (requireNamespace("franzpak", quietly = TRUE)) {
  franzpak::coef_table_mplus(fitted$results, constraints = TRUE)
} else {
  knitr::asis_output(
    "*Install `franzpak` with `pak::pak(\"franciscowilhelm/franzpak\")` to render this table.*"
  )
}
```

*Install `franzpak` with `pak::pak("franciscowilhelm/franzpak")` to
render this table.*

Table 1: Unstandardized Mplus estimates and model constraints.

Set `plot = TRUE` to pass the polynomial coefficients to
[`RSA::plotRSA()`](https://rdrr.io/pkg/RSA/man/plotRSA.html). Additional
arguments are forwarded to that function.

``` r

plotted <- RSA_mplus(
  fitted,
  plot = TRUE,
  xlab = "Actual (X)",
  ylab = "Desired (Y)",
  zlab = "Outcome (Z)"
)

plotted$plot
```

![Three-dimensional response surface with lines of congruence and
incongruence and starred a1 through a5
parameters.](mplus-response-surfaces_files/figure-html/fig-mplus-rsa-1.png)

Figure 2: Latent response surface fitted in Mplus. Stars mark Mplus
constraint p-values (p ≤ .05, .01, and .001).

The surface itself is defined by `b0` and the five polynomial
coefficients. For the default three-dimensional plot, rsahelpers
replaces RSA’s calculated parameter text with the `CS`, `CC`, `IS`,
`IC`, and `A5` estimates from Mplus. It adds `*`, `**`, and `***` for
p-values no greater than .05, .01, and .001.

With `ESTIMATOR = BAYES`, Mplus reports a 95% credibility interval next
to a one-tailed posterior p-value. The credibility interval is the
preferred decision rule, so rsahelpers ignores the posterior p-value and
adds a single `*` to surface parameters whose interval excludes zero —
matching the `*` column of the Mplus output itself. The interval bounds
are also returned in `new_parameters` as `lower_2.5ci` and
`upper_2.5ci`.

If these constraint estimates are unavailable, RSA’s original unstarred
text is retained. `param = FALSE` hides the annotation, and contour
plots do not display a parameter block.

Generated latent and SI-LMS models fix the latent outcome mean to zero,
so their plotting intercept defaults to `b0 = 0`.

`coef_type` selects the Mplus coefficient table:

``` r

RSA_mplus(fitted, coef_type = "un", plot = FALSE)
```

    $plot
    NULL

    $coefficients
         x      y     x2     xy     y2     b0
     0.287  0.260 -0.180  0.220 -0.086  0.000

    $regression_parameters
      Label    est    se  pval term
    1  Z<-X  0.287 0.054 0.000    x
    2  Z<-Y  0.260 0.044 0.000    y
    3 Z<-XS -0.180 0.043 0.000   x2
    4 Z<-XY  0.220 0.057 0.000   xy
    5 Z<-YS -0.086 0.029 0.003   y2

    $new_parameters
       Label    est     se  pval
    1     CS  0.548  0.052 0.000
    2     CC -0.046  0.033 0.165
    3     IS  0.027  0.084 0.747
    4     IC -0.486  0.110 0.000
    5     A5 -0.095  0.046 0.039
    6     X0  8.011  8.955 0.371
    7     Y0 11.816 13.514 0.382
    8    P10 -0.348  0.306 0.255
    9    P11  1.519  0.307 0.000
    10   P20 17.091 19.379 0.378
    11   P21 -0.659  0.133 0.000

    $model
    RSA model generated by rsahelpers

    Estimated using MLR
    Number of obs: 300, number of (free) parameters: 33

    Fit Indices:

    CFI = NA, TLI = NA, SRMR = NA
    RMSEA = NA, 90% CI [NA, NA], p < .05 = NA
    AIC = 5165.521, BIC = 5287.746
    NULL

    $workflow
    <rsa_mplus_workflow>
    Status: read
    Model type: latent
    Output: /home/runner/work/_temp/Library/rsahelpers/extdata/mplus-response-surfaces.out

    $outcome
    [1] "Z"

    $coefficient_type
    [1] "un"

    attr(,"class")
    [1] "rsa_mplus"

``` r

RSA_mplus(fitted, coef_type = "stdyx", plot = FALSE)
```

    $plot
    NULL

    $coefficients
         x      y     x2     xy     y2     b0
     0.314  0.354 -0.179  0.272 -0.132  0.000

    $regression_parameters
      Label    est    se  pval term
    1  Z<-X  0.314 0.054 0.000    x
    2  Z<-Y  0.354 0.059 0.000    y
    3 Z<-XS -0.179 0.043 0.000   x2
    4 Z<-XY  0.272 0.068 0.000   xy
    5 Z<-YS -0.132 0.046 0.004   y2

    $new_parameters
    NULL

    $model
    RSA model generated by rsahelpers

    Estimated using MLR
    Number of obs: 300, number of (free) parameters: 33

    Fit Indices:

    CFI = NA, TLI = NA, SRMR = NA
    RMSEA = NA, 90% CI [NA, NA], p < .05 = NA
    AIC = 5165.521, BIC = 5287.746
    NULL

    $workflow
    <rsa_mplus_workflow>
    Status: read
    Model type: latent
    Output: /home/runner/work/_temp/Library/rsahelpers/extdata/mplus-response-surfaces.out

    $outcome
    [1] "Z"

    $coefficient_type
    [1] "stdyx"

    attr(,"class")
    [1] "rsa_mplus"

Unstandardized coefficients are the default and are generally the
clearest choice when the response surface is interpreted in the model’s
working scale.

## 11 Choosing between latent LMS and SI-LMS

| Consideration | Latent LMS | SI-LMS |
|----|----|----|
| Indicators | Two or more per RSA factor | One prepared score per factor |
| Measurement error | Estimated from the measurement model | Fixed from reliability and score variance |
| Input required | Item-level columns | Scores plus reliability estimates |
| Flexibility | Factor-specific loadings and residuals | Simpler, fewer measurement parameters |
| Typical use | Primary latent-variable analysis | Too few indicators or convergence difficulties |

The two approaches estimate the same five-term polynomial structure.
They differ in how measurement error is represented, so switching
between them is a modeling decision rather than only a computational
shortcut. Su et al. (2019) recommend using LMS as the first choice.

## 12 Troubleshooting

- **Mplus is not found:** pass the executable through `Mplus_command` or
  add it to the system path.
- **The model has not been written:** call
  [`write_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md)
  before the standalone runner, or use
  [`rsa_mplus_workflow()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md).
- **A role cannot be inferred:** workflow metadata is available only for
  models created by
  [`create_rsa_mplus_model()`](https://franciscowilhelm.github.io/rsahelpers/reference/create_rsa_mplus_model.md);
  supply all six labels for an arbitrary output.
- **Mplus names collide:** rename variables that share their first eight
  characters.
- **Principal-axis constraints fail:** use `constraints = "surface"` and
  inspect whether the fitted polynomial has singular principal-axis
  formulas.
- **SI-LMS construction fails:** confirm that each factor has exactly
  one numeric score and that reliability is named by the latent factors.
- **A model file already exists:** choose a new path or explicitly set
  `overwrite = TRUE`.

The generated object deliberately exposes `mplus`, `tidysem`, `files`,
and `results`, so the exact syntax, paths, parsed warnings, and
estimates remain available when diagnosing a model.

## 13 References

Su, R., Zhang, Q., Liu, Y., & Tay, L. (2019). Modeling congruence in
organizational research with latent moderated structural equations.
*Journal of Applied Psychology, 104*(11), 1404–1433.
<https://doi.org/10.1037/apl0000411>
