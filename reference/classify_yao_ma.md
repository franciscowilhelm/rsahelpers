# Classify a response surface using the Yao-Ma typology

Classify a response surface using the Yao-Ma typology

## Usage

``` r
classify_yao_ma(x, valence, equivalence = NULL, level = 0.95)
```

## Arguments

- x:

  An
  [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  result, an `rsa_mplus_workflow`, an `mplus.model`, an `mplusObject`,
  or a data frame containing response-surface parameters. A data frame
  must have columns `parameter`, `estimate`, `conf.low`, and
  `conf.high`.

- valence:

  Whether higher (`"positive"`) or lower (`"negative"`) values of the
  outcome are desirable. This argument is required because it selects
  the relevant principal axis.

- equivalence:

  Optional named numeric vector of practical-equivalence margins. Supply
  one nonnegative margin for each of `CS`, `CC`, `IS`, `IC`,
  `axis_intercept`, and `axis_slope`. When omitted, confidence intervals
  that include the null are treated as absence of an effect. When
  supplied, each condition is classified as equivalent, different, or
  indeterminate.

- level:

  Confidence level used when extracting intervals from Mplus. Intervals
  already present in a parameter data frame are used unchanged.

## Value

An object of class `yao_ma_classification`. It contains a one-row
`summary`, the normalized `parameters`, and a condition-level
`decisions` data frame.

## Details

The function uses `CS`, `CC`, `IS`, and `IC`, corresponding to the usual
RSA parameters `a1` through `a4`. For positively valenced outcomes it
uses `P10` and `P11`; for negatively valenced outcomes it uses `P20` and
`P21`.

Principal-axis parameters are nonlinear functions of the polynomial
coefficients. Mplus confidence intervals for these `MODEL CONSTRAINT`
parameters use its delta-method standard errors. A missing or degenerate
principal axis therefore produces a partial result rather than an
invented Yao-Ma subtype.

## Examples

``` r
params <- data.frame(
  parameter = c("CS", "CC", "IS", "IC", "P10", "P11"),
  estimate = c(0.4, 0.01, 0.02, -0.5, 0.03, 1.02),
  conf.low = c(0.2, -0.1, -0.1, -0.7, -0.1, 0.9),
  conf.high = c(0.6, 0.1, 0.1, -0.3, 0.1, 1.1)
)
classify_yao_ma(params, valence = "positive")
#> <yao_ma_classification>
#> Method: significance 
#> Valence: positive 
#> Status: classified 
#> Type: 2 - Exact correspondence & LLE 
```
