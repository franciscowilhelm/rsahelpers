# Spline Response Surfaces

## 1 Overview

Congruence spline regression represents a response surface with linear
regions joined at one or two estimated seam lines. Unlike polynomial
response surface analysis, the fitted surface can change slope abruptly
where it crosses a seam.

This vignette uses the workshop data distributed with `rsahelpers` to
show how to:

1.  prepare commensurate predictors;
2.  fit the OLS models used for comparison and starting values;
3.  estimate and interpret one- and two-seam spline surfaces;
4.  compare constrained and unconstrained models;
5.  visualize the fitted surface; and
6.  obtain bootstrap confidence intervals.

## 2 Formula convention

All spline entry points use the same formula structure:

``` r

Z ~ X * Y
```

The response maps to Z, the first predictor to X, and the second
predictor to Y. The `*` declares the X and Y roles. It does not request
the ordinary `X:Y` interaction used by
[`lm()`](https://rdrr.io/r/stats/lm.html).

> **Important:** Reversing the predictors reverses their X and Y roles.
> This matters because the estimated seam is written as a function of X,
> $`Y = c_0 + c_1X`$.

The main example models job satisfaction (`JOBSAT`) from the uncentered
authority variables `ATHW` (X) and `ATHH` (Y):

``` r

authority_formula <- JOBSAT ~ ATHW * ATHH
```

The workshop data also contain precomputed variables with `*WC` and
`*HC` suffixes. Here we deliberately use the raw columns so the package
applies its default pooled centering and scaling, which preserve the
meaning of the line $`X = Y`$.

## 3 Prepare the model data

[`prepare_congruence_data()`](https://franciscowilhelm.github.io/rsahelpers/reference/prepare_congruence_data.md)
applies complete-case selection, centering, and scaling, then creates
the hinge variables needed by the comparison models.

``` r

authority_data <- prepare_congruence_data(
  authority_formula,
  data = workshop
)

head(authority_data[c("x", "y", "z", "abs_diff", "zd")])
```

               x         y        z abs_diff zd
    1 -0.5086229 0.6445793 3.666667 1.153202  0
    2 -0.1242221 1.0289801 3.333333 1.153202  0
    3 -0.5086229 1.7977815 3.000000 2.306404  0
    4  0.2601786 1.7977815 3.666667 1.537603  0
    5 -0.5086229 0.6445793 4.000000 1.153202  0
    6  0.6445793 1.7977815 4.000000 1.153202  0

The standardized output always uses `x`, `y`, and `z` internally. The
source variables and transformation constants remain attached as
attributes.

``` r

attr(authority_data, "variables")
```

           z        x        y
    "JOBSAT"   "ATHW"   "ATHH" 

``` r

c(
  x_center = attr(authority_data, "x_center"),
  y_center = attr(authority_data, "y_center"),
  scale = attr(authority_data, "scale")
)
```

     x_center  y_center     scale
    3.4410526 3.4410526 0.8671506 

## 4 Fit comparison models

[`fit_piecewise_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_piecewise_congruence.md)
fits the absolute-difference, linear, one-break, and constrained
one-seam OLS models. Requesting two seams also adds a fixed two-seam
comparison.

``` r

authority_ols <- fit_piecewise_congruence(
  authority_formula,
  data = workshop,
  n_seams = 2
)

tidy_piecewise_summary(authority_ols)
```

                      model        term    estimate  std.error   statistic
    1   absolute_difference (Intercept)  3.96663155 0.04143884  95.7225500
    2   absolute_difference    abs_diff -0.27434173 0.02921700  -9.3897991
    3                linear (Intercept)  3.68213604 0.02940738 125.2112919
    4                linear           x  0.19331465 0.03467558   5.5749509
    5                linear           y  0.17739435 0.02543773   6.9736716
    6             one_break (Intercept)  3.95830290 0.04894892  80.8659874
    7             one_break           x  0.47750398 0.05329113   8.9602894
    8             one_break           y -0.19635081 0.05204842  -3.7724648
    9             one_break           w -0.06465692 0.08900720  -0.7264235
    10            one_break          xw -0.43588391 0.08011661  -5.4406188
    11            one_break          yw  0.56465105 0.07116440   7.9344595
    12 constrained_one_seam (Intercept)  3.93194618 0.04063923  96.7524857
    13 constrained_one_seam           x  0.49544099 0.04855548  10.2036051
    14 constrained_one_seam           y -0.16089818 0.04644013  -3.4646367
    15 constrained_one_seam          zd  0.52124670 0.06077270   8.5769880
    16       fixed_two_seam (Intercept)  4.48988102 0.13394865  33.5194206
    17       fixed_two_seam           x  0.79364338 0.08773341   9.0460790
    18       fixed_two_seam           y -0.45199947 0.08854476  -5.1047566
    19       fixed_two_seam hinge_upper  0.66189240 0.12014453   5.5091349
    20       fixed_two_seam hinge_lower  0.18167810 0.09016325   2.0149905
             p.value  r.squared
    1   0.000000e+00 0.08509074
    2   4.333909e-20 0.08509074
    3   0.000000e+00 0.07892639
    4   3.228903e-08 0.07892639
    5   5.789407e-12 0.07892639
    6   0.000000e+00 0.14784308
    7   1.697415e-18 0.14784308
    8   1.717102e-04 0.14784308
    9   4.677592e-01 0.14784308
    10  6.766895e-08 0.14784308
    11  5.964800e-15 0.14784308
    12  0.000000e+00 0.14538468
    13  2.906905e-23 0.14538468
    14  5.547810e-04 0.14538468
    15  3.938720e-17 0.14538468
    16 6.101550e-163 0.15307043
    17  8.249184e-19 0.15307043
    18  4.004505e-07 0.15307043
    19  4.649675e-08 0.15307043
    20  4.418843e-02 0.15307043

These models provide interpretable baselines and starting values for
nonlinear estimation. They are not substitutes for an estimated-seam
spline.

## 5 Estimate a one-seam surface

``` r

authority_one <- fit_spline_congruence(
  authority_formula,
  data = workshop,
  n_seams = 1
)

summary(authority_one, warn_seam = FALSE)
```

    Congruence spline regression
      seams: 1
      centering: pooled  scaling: pooled
      solver: minpack.lm::nlsLM

       estimate std.error statistic   p.value
    b0  4.34492   0.13952   31.1416 < 2.2e-16 ***
    b1  0.67543   0.09806    6.8879 1.033e-11 ***
    b2 -0.42723   0.10011   -4.2674 2.177e-05 ***
    b3  0.76195   0.10522    7.2418 9.181e-13 ***
    c0  0.62838   0.13858    4.5343 6.521e-06 ***
    c1  0.86158   0.14142    6.0926 1.616e-09 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual SE on 944 df | R-squared: 0.1554 | RSS: 677.4 | AIC: 2388.66 | n: 950

``` r

coef(authority_one)
```

            b0         b1         b2         b3         c0         c1
     4.3449150  0.6754307 -0.4272323  0.7619538  0.6283846  0.8615831 

For a one-seam model, `c0` and `c1` define the seam $`Y = c_0 + c_1X`$.
The `b` coefficients define the plane on one side of the seam and the
change in slope after the seam is crossed.

[`spline_surface_features()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_surface_features.md)
translates these coefficients into region slopes, seam features,
symmetry expressions, and shifts along selected lines parallel to the
line of incongruence.

``` r

spline_surface_features(authority_one)
```

         right_intercept        right_x_slope        right_y_slope
              3.86611490           0.01894407           0.33472155
          seam_intercept           seam_slope  equal_opposite_left
              4.07644877           0.30733452           0.24819836
    equal_opposite_right           symmetry_x           symmetry_y
              0.35366562           0.69437472          -0.09251074
            shift_y=-2-x          shift_y=0-x          shift_y=2-x
             -0.58252625          -0.47737328          -0.37222030 

[`spline_tests()`](https://franciscowilhelm.github.io/rsahelpers/reference/spline_tests.md)
supplies the corresponding scalar and joint Wald tests.

``` r

authority_tests <- spline_tests(authority_one, warn_seam = FALSE)
authority_tests$scalar
```

                       term    estimate         se  statistic      p.value
    1       right_intercept  3.86611490 0.03944758 98.0064000 0.000000e+00
    2         right_x_slope  0.01894407 0.04203462  0.4506778 6.523253e-01
    3         right_y_slope  0.33472155 0.03236368 10.3425057 8.002019e-24
    4        seam_intercept  4.07644877 0.04989713 81.6970607 0.000000e+00
    5            seam_slope  0.30733452 0.04962670  6.1929261 8.801467e-10
    6            shift_y=-x  0.47737328 0.10731399  4.4483788 9.681899e-06
    7           shift_y=2-x  0.37222030 0.13883675  2.6809926 7.468422e-03
    8          shift_y=-2-x  0.58252625 0.17435545  3.3410270 8.674253e-04
    9   equal_opposite_left  0.24819836 0.10580028  2.3459141 1.918658e-02
    10 equal_opposite_right  0.35366562 0.04456064  7.9367263 5.863535e-15
    11           symmetry_x  0.69437472 0.10668967  6.5083594 1.231699e-10
    12           symmetry_y -0.09251074 0.10521606 -0.8792454 3.794920e-01

``` r

authority_tests$joint
```

                                 term df statistic      p.value
    1          deviation_from_no_seam  3  53.44753 6.534274e-32
    2 absolute_difference_constraints  4  25.40216 5.192819e-20
    3          seam_equals_y_equals_x  2  11.42475 1.251392e-05

> **Warning:** Delta-method standard errors involving an estimated seam
> are approximate because the fitted surface is not differentiable at
> the seam. Use the bootstrap workflow below for final inference about
> seam features.

## 6 Test a theoretically fixed seam

A seam fixed to the line $`Y = X`$ has `c0 = 0` and `c1 = 1`. The `fix`
argument removes those values from optimization, producing a constrained
model that is nested in the freely estimated one-seam spline.

``` r

authority_fixed <- fit_spline_congruence(
  authority_formula,
  data = workshop,
  n_seams = 1,
  fix = c(c0 = 0, c1 = 1)
)

compare_spline_models(
  fixed_seam = authority_fixed,
  estimated_seam = authority_one
)
```

               model npar df.residual      rss r.squared      AIC df    deltaR2
    1     fixed_seam    4         946 685.4264 0.1453847 2395.883 NA         NA
    2 estimated_seam    6         944 677.3760 0.1554222 2388.660  2 0.01003752
             F     p.value
    1       NA          NA
    2 5.609558 0.003785848

The comparison table reports fit statistics and the change in residual
sum of squares. Models supplied to
[`compare_spline_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_spline_models.md)
must be ordered from simpler to more complex.

## 7 Estimate a two-seam surface

``` r

authority_two <- fit_spline_congruence(
  authority_formula,
  data = workshop,
  n_seams = 2
)

summary(authority_two, warn_seam = FALSE)
```

    Congruence spline regression
      seams: 2
      centering: pooled  scaling: pooled
      solver: minpack.lm::nlsLM

        estimate std.error statistic   p.value
    b0   4.33467   0.17402   24.9095 < 2.2e-16 ***
    b1   0.67148   0.10334    6.4975 1.321e-10 ***
    b2  -0.41956   0.11667   -3.5962 0.0003396 ***
    b3   0.61715   0.17288    3.5698 0.0003753 ***
    b4   0.19955   0.14328    1.3927 0.1640472
    c10  0.75726   0.22723    3.3326 0.0008939 ***
    c11  0.90475   0.19737    4.5841 5.173e-06 ***
    c20 -0.42760   0.54385   -0.7862 0.4319230
    c21  0.47065   0.44086    1.0676 0.2859834
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual SE on 941 df | R-squared: 0.1573 | RSS: 675.8 | AIC: 2392.51 | n: 950

``` r

spline_surface_features(authority_two)
```

    seam1_shift_y_neg_x seam2_shift_y_neg_x        base_x_slope        base_y_slope
             0.56224318         -0.41119112          0.67148043         -0.41956029
          seam1_x_slope       seam1_y_slope       seam2_x_slope       seam2_y_slope
             0.11310996          0.19759461          0.57756258         -0.22001119
           both_x_slope        both_y_slope          crossing_x         seams_cross
             0.01919211          0.39714371         -2.72947855          1.00000000
             n_sections
             3.00000000 

``` r

anova(authority_one, authority_two)
```

       model npar df.residual      rss r.squared      AIC df     deltaR2         F
    1 model1    6         944 677.3760 0.1554222 2388.660 NA          NA        NA
    2 model2    9         941 675.8475 0.1573280 2392.513  3 0.001905844 0.7094099
        p.value
    1        NA
    2 0.5464788

Two-seam models can have local solutions or label-switched seams. By
default,
[`fit_spline_congruence()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_spline_congruence.md)
tries several starting values and retains the converged solution with
the lowest residual sum of squares. A more complicated model should not
be retained merely because it contains another seam.

## 8 Plot the fitted surface

A contour plot makes the seam location and supported regions easy to
read.

``` r

plot_spline_contour(
  authority_one,
  xlab = "Authority X",
  ylab = "Authority Y",
  main = "Authority: one-seam spline"
)
```

![Filled contour plot of predicted job satisfaction over X and Y, with
an estimated seam and a dashed diagonal
line.](spline-response-surfaces_files/figure-html/fig-authority-contour-1.png)

Figure 1: One-seam authority surface with the estimated seam and the
line X = Y.

The three-dimensional view emphasizes changes in the surface slope.

``` r

plot_spline_surface(
  authority_one,
  xlab = "Authority X",
  ylab = "Authority Y",
  zlab = "Job satisfaction",
  main = "Authority: one-seam spline"
)
```

![Three-dimensional tiled response surface for job satisfaction with an
estimated seam
line.](spline-response-surfaces_files/figure-html/fig-authority-surface-1.png)

Figure 2: Three-dimensional view of the one-seam authority surface.

## 9 Bootstrap seam inference

Edwards and Parry recommend bootstrap intervals for the nonlinear seam
parameters and derived surface features. A final analysis should use a
large number of resamples; the code is not executed during vignette
rendering.

``` r

set.seed(2026)
authority_boot <- bootstrap_spline(
  authority_one,
  R = 10000,
  type = "bca"
)
authority_boot
attr(authority_boot, "fail_rate")
```

Inspect the reported failure rate. A substantial rate indicates that
bootstrap samples frequently do not support a stable nonlinear solution.

## 10 Repeat the workflow for another predictor pair

The variety variables use exactly the same Z ~ X \* Y logic.

``` r

variety_formula <- JOBSAT ~ VARW * VARH

variety_one <- fit_spline_congruence(
  variety_formula,
  data = workshop
)

summary(variety_one, warn_seam = FALSE)
```

    Congruence spline regression
      seams: 1
      centering: pooled  scaling: pooled
      solver: minpack.lm::nlsLM

        estimate std.error statistic   p.value
    b0  3.896197  0.046371   84.0231 < 2.2e-16 ***
    b1  0.333611  0.073557    4.5354 6.490e-06 ***
    b2 -0.156410  0.076673   -2.0400   0.04163 *
    b3  0.743771  0.095365    7.7992 1.644e-14 ***
    c0 -0.297078  0.117450   -2.5294   0.01159 *
    c1  0.624971  0.091152    6.8564 1.275e-11 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual SE on 944 df | R-squared: 0.1770 | RSS: 660.1 | AIC: 2364.06 | n: 950

``` r

spline_surface_features(variety_one)
```

         right_intercept        right_x_slope        right_y_slope
              4.11715510          -0.13122404           0.58736084
          seam_intercept           seam_slope  equal_opposite_left
              3.94266307           0.23585933           0.17720109
    equal_opposite_right           symmetry_x           symmetry_y
              0.45613680           0.20238685           0.43095104
            shift_y=-2-x          shift_y=0-x          shift_y=2-x
             -0.06784093           0.25854733           0.58493558 

Keeping the formula visible in each analysis makes the X/Y orientation
auditable while the remaining spline workflow stays unchanged.
