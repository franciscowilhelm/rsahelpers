# Plot a fitted congruence spline surface in 3D

Draw a three-dimensional response surface for a fitted Edwards-Parry
spline model. The surface is evaluated on a regular grid over the
observed `x` and `y` ranges and plotted with base R's
[`persp()`](https://rdrr.io/r/graphics/persp.html). By default the
facets are shaded by predicted outcome on a high-contrast diverging
palette (as the `RSA` package does) so the surface shape is legible, and
the estimated seam line(s) are projected onto it.

## Usage

``` r
plot_spline_surface(
  fit,
  grid_size = 25,
  xlim = NULL,
  ylim = NULL,
  equal_limits = TRUE,
  theta = -35,
  phi = 25,
  expand = 0.65,
  color_by = c("outcome", "none"),
  palette = c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b", "#ffffbf",
    "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"),
  n_color = 16,
  col = "lightblue",
  border = NA,
  ticktype = "detailed",
  xlab = "X (centered)",
  ylab = "Y (centered)",
  zlab = "Outcome",
  main = "Congruence spline surface",
  show_seams = TRUE,
  show_fit_line = TRUE,
  show_congruence = NULL,
  show_incongruence = FALSE,
  seam_col = "#2166ac",
  seam_lwd = 2.5,
  congruence_col = "grey25",
  congruence_lty = 2,
  congruence_lwd = 1.2,
  incongruence_col = "grey45",
  incongruence_lty = 3,
  incongruence_lwd = 1.1,
  ...
)
```

## Arguments

- fit:

  A `congruence_spline` object.

- grid_size:

  Integer number of grid points per axis. The default of `25` yields
  larger, higher-contrast facet tiles than a dense grid.

- xlim, ylim:

  Optional axis limits. If either is `NULL` and `equal_limits = TRUE`,
  both axes use one shared symmetric range.

- equal_limits:

  Logical. If `TRUE`, use equal limits for the `x` and `y` axes so the
  line of congruence is visually anchored by points such as `(-2, -2)`
  and `(2, 2)`.

- theta, phi:

  Viewing angles passed to
  [`persp()`](https://rdrr.io/r/graphics/persp.html).

- expand:

  Expansion factor passed to
  [`persp()`](https://rdrr.io/r/graphics/persp.html).

- color_by:

  Facet shading. `"outcome"` (default) colors each facet by its
  predicted outcome; `"none"` uses the flat `col` fill.

- palette:

  Either a vector of colors to interpolate or the name of a
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
  palette used when `color_by = "outcome"`.

- n_color:

  Number of color bins for outcome shading.

- col:

  Flat surface fill color used when `color_by = "none"`.

- border:

  Facet border color. Defaults to `NA` (no border) so colored tiles read
  as solid blocks.

- ticktype:

  Tick type passed to
  [`persp()`](https://rdrr.io/r/graphics/persp.html).

- xlab, ylab, zlab:

  Axis labels.

- main:

  Plot title.

- show_seams:

  Logical. If `TRUE`, draw the estimated seam line(s) (the fitted ridge
  of the surface).

- show_fit_line:

  Logical. If `TRUE` (default), draw the purely theoretical line of
  congruence, `X = Y`, on the surface. Set to `FALSE` to drop it for
  models where congruence along `X = Y` is not the hypothesis.

- show_congruence:

  Deprecated alias for `show_fit_line`, retained for backward
  compatibility.

- show_incongruence:

  Logical. If `TRUE`, draw the line of incongruence, `X = -Y`, on the
  fitted surface.

- seam_col, seam_lwd:

  Seam line color and width.

- congruence_col, congruence_lty, congruence_lwd:

  Line style for `X = Y`.

- incongruence_col, incongruence_lty, incongruence_lwd:

  Line style for `X = -Y`.

- ...:

  Additional arguments passed to
  [`persp()`](https://rdrr.io/r/graphics/persp.html).

## Value

Invisibly returns a list containing the grid vectors, fitted surface
matrix, [`persp()`](https://rdrr.io/r/graphics/persp.html)
transformation matrix, and plotted seam coordinates.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_spline_surface(fit)
plot_spline_surface(fit, show_fit_line = FALSE)
} # }
```
