# Diagnose joint predictor coverage for response-surface analysis

Diagnose joint predictor coverage for response-surface analysis

## Usage

``` r
predictor_coverage(
  data,
  x,
  y,
  breaks = 10L,
  discrepancy = NULL,
  min_count = 5L
)
```

## Arguments

- data:

  A data frame containing the two commensurate predictors.

- x, y:

  Predictor columns, supplied as unquoted names or strings.

- breaks:

  Either the number of shared equal-width bins (default `10`) or a
  numeric vector of common breakpoints for both predictors.

- discrepancy:

  Optional positive cutoffs at which to summarize the proportions of
  observations with `X - Y` or `Y - X` at least that large.

- min_count:

  Minimum cell count considered adequately covered.

## Value

An object of class `predictor_coverage`. It contains `cells`, a one-row
`summary`, an optional `discrepancy` table, the shared `breaks`, and
plotting metadata.

## Details

Shared breaks preserve the geometry of the line `X = Y`. This is
especially useful for scale scores derived from Likert items: using
every observed decimal value creates an unnecessarily sparse grid, while
separate quantile bins distort congruence and discrepancy regions.

## Examples

``` r
scores <- data.frame(
  self = c(1.0, 1.5, 2.2, 3.1, 3.8, 4.4, 4.8),
  other = c(1.2, 1.6, 2.7, 2.9, 3.5, 4.1, 4.7)
)
coverage <- predictor_coverage(
  scores,
  self,
  other,
  breaks = 4,
  discrepancy = c(0.5, 1)
)
coverage$summary
#>      x     y n n_missing correlation x_greater_than_y y_greater_than_x ties
#> 1 self other 7         0   0.9858663        0.5714286        0.4285714    0
#>   empty_cells sparse_cells min_cell_count min_count
#> 1          11           16              0         5
```
