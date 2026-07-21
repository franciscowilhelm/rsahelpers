# Plot a fitted congruence spline surface as a 2D contour map

Draw the fitted spline surface as a filled `ggplot2` heatmap with
contour lines, the estimated seam line(s), and (optionally) the line of
congruence. A 2D view is often clearer than
[`plot_spline_surface()`](https://franciscowilhelm.github.io/rsahelpers/reference/plot_spline_surface.md)
for reading where the seams fall. Requires the `ggplot2` package.

## Usage

``` r
plot_spline_contour(
  fit,
  grid_size = 80,
  xlim = NULL,
  ylim = NULL,
  equal_limits = TRUE,
  palette = c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee08b", "#ffffbf",
    "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#006837"),
  bins = 12,
  show_seams = TRUE,
  show_fit_line = TRUE,
  seam_col = "#2166ac",
  seam_lwd = 1.1,
  fit_line_col = "grey25",
  fit_line_lty = "dashed",
  xlab = "X (centered)",
  ylab = "Y (centered)",
  fill_lab = "Outcome",
  main = "Congruence spline surface"
)
```

## Arguments

- fit:

  A `congruence_spline` object.

- grid_size:

  Integer number of grid points per axis.

- xlim, ylim:

  Optional axis limits. If either is `NULL` and `equal_limits = TRUE`,
  both axes use one shared symmetric range.

- equal_limits:

  Logical. If `TRUE`, use one shared symmetric range for both axes.

- palette:

  Either a vector of colors to interpolate or the name of a
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
  palette for the fill gradient.

- bins:

  Number of contour bins.

- show_seams:

  Logical. If `TRUE`, overlay the estimated seam line(s).

- show_fit_line:

  Logical. If `TRUE` (default), overlay the theoretical line of
  congruence `X = Y`.

- seam_col, seam_lwd:

  Seam line color and width.

- fit_line_col, fit_line_lty:

  Line style for `X = Y`.

- xlab, ylab, fill_lab, main:

  Axis, legend, and title labels.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_spline_contour(fit)
} # }
```
