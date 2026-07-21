# Extract a congruence spline covariance matrix

Extract a congruence spline covariance matrix

## Usage

``` r
# S3 method for class 'congruence_spline'
vcov(object, ...)
```

## Arguments

- object:

  A `congruence_spline` object.

- ...:

  Unused.

## Value

A numeric covariance matrix. Values are `NA` if `numDeriv` is not
installed.
