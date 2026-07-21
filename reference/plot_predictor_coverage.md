# Plot joint predictor coverage

Plot joint predictor coverage

## Usage

``` r
plot_predictor_coverage(coverage, mark_sparse = TRUE, show_diagonal = TRUE)
```

## Arguments

- coverage:

  A
  [`predictor_coverage()`](https://franciscowilhelm.github.io/rsahelpers/reference/predictor_coverage.md)
  result or a named list of results.

- mark_sparse:

  If `TRUE`, mark cells below the result's minimum count.

- show_diagonal:

  If `TRUE`, draw the line `X = Y`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  scores <- data.frame(
    x = c(1, 1.5, 2, 3, 4, 4.5),
    y = c(1.2, 1.4, 2.5, 2.8, 3.7, 4.6)
  )
  plot_predictor_coverage(predictor_coverage(scores, x, y, breaks = 4))
}
```
