# Log-likelihood of a fitted congruence spline

Gaussian log-likelihood implied by the residual sum of squares, enabling
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) for `congruence_spline`
objects.

## Usage

``` r
# S3 method for class 'congruence_spline'
logLik(object, ...)
```

## Arguments

- object:

  A `congruence_spline` object.

- ...:

  Unused.

## Value

An object of class `logLik` with `df` and `nobs` attributes. The degrees
of freedom count the free (estimated) parameters plus the residual
variance.
