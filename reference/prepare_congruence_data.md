# Prepare variables for congruence regression models

Create the centered component variables and auxiliary terms used by
Edwards-Parry absolute-difference, piecewise, and spline congruence
models.

## Usage

``` r
prepare_congruence_data(
  formula,
  data,
  center = c("pooled", "variablewise", "none"),
  scale = c("pooled", "none"),
  hinge_offset = 1
)
```

## Arguments

- formula:

  A two-sided formula of the form `z ~ x * y`. The response is mapped to
  Z, the first predictor to X, and the second predictor to Y. The `*`
  declares the two predictor roles; it does not add an ordinary linear
  interaction term to the spline model.

- data:

  A data frame containing the Z, X, and Y variables.

- center:

  Centering of the component variables, following the `RSA` package. One
  of `"pooled"` (subtract the shared mean of `x` and `y`; default),
  `"variablewise"` (subtract each component's own mean), or `"none"`.
  Logical values are accepted for backward compatibility: `TRUE` maps to
  `"variablewise"`, `FALSE` to `"none"`. Pooled centering preserves the
  `X = Y` congruence interpretation; use `"none"` for columns that are
  already centered, such as the downloaded `*HC` and `*WC` columns.

- scale:

  Scaling of the component variables. One of `"pooled"` (divide both `x`
  and `y` by their pooled standard deviation; default) or `"none"`.

- hinge_offset:

  Numeric half-width, in **working-scale** units, of the fixed two-seam
  hinges `hinge_upper`/`hinge_lower` (`Y = X ± hinge_offset`). Because
  the offset is applied after centering and scaling, under pooled
  `scale = "pooled"` an offset of `1` corresponds to one pooled standard
  deviation rather than one raw scale unit.

## Value

A data frame with standardized columns `x`, `y`, `z`, `.row`, and helper
columns for absolute difference, one-break piecewise regression,
constrained piecewise regression, and fixed two-seam comparisons. The
centering and scaling constants are recorded as attributes (`x_center`,
`y_center`, `scale`, `center_method`, `scale_method`).

## Examples

``` r
if (FALSE) { # \dontrun{
prepared <- prepare_congruence_data(satisfaction ~ x * y, my_data)
} # }
```
