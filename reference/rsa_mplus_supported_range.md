# Supported region of a latent response surface

Gives the range of the line of congruence and the line of incongruence
that the data support, from the model-implied distribution of the latent
X and Y factors. Factor scores are shrunk towards the mean, so their
range is narrower than the latent distribution to which the surface
refers; this function uses the estimated latent means, variances, and
covariance instead and assumes bivariate normality.

## Usage

``` r
rsa_mplus_supported_range(
  model,
  pred_x = NULL,
  pred_y = NULL,
  probs = c(0.05, 0.95)
)
```

## Arguments

- model:

  A fitted model as in
  [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  (a workflow, an Mplus output path, or a parsed Mplus model), or a list
  of them. The moments of several models are averaged, so that surfaces
  of the same X and Y factors for different outcomes share one region.

- pred_x, pred_y:

  Names of the X and Y factors. Taken from a workflow when `NULL`.

- probs:

  Lower and upper quantile that bound the region.

## Value

A list with the quantiles `loc` and `loic`, and `moments`, the
(averaged) latent means `mx` and `my`, variances `vx` and `vy`, and
covariance `cxy`.

## Details

A point `(X, Y)` is at position `(X + Y)/2` along the line of congruence
and `(X - Y)/2` along the line of incongruence, the coordinates in which
`a1` to `a4` describe the surface. Their means and variances follow from
the latent moments. A factor mean that the model does not estimate is 0,
the Mplus default.

## See also

[`rsa_mplus_draws()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_draws.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rsa_mplus_supported_range(fitted)
rsa_mplus_supported_range("model.out", pred_x = "X", pred_y = "Y")
} # }
```
