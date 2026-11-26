# Roadmap

## Latent RSA model family

Goal: the procedure of the `RSA` package and of Schönbrodt, *Testing fit
patterns with polynomial regression models*, as closely as possible in
Mplus (MLR with `XWITH`).

Done (0.2.0):

- Effects coding and pooled centering of the latent X and Y factors.
- [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
  reads restricted models (`rsa_model`).
- `create_rsa_mplus_model(rsa_model =, start =)` generates the 16
  second-order models, including SRRR and SRSQD with the ridge
  parameters C, S, bM, and a4’.
- [`fit_rsa_mplus_family()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_family.md)
  fits the family from one specification (by default
  [`RSA::RSA()`](https://rdrr.io/pkg/RSA/man/RSA.html)’s 14 models, with
  start values from the full model) and keeps going when a model fails.
  [`compare_rsa_mplus_models()`](https://franciscowilhelm.github.io/rsahelpers/reference/compare_rsa_mplus_models.md)
  ranks the models by AICc (Akaike weights, evidence ratios) and adds
  the Satorra-Bentler scaled LRT against the full model from the MLR
  log-likelihoods and scaling factors, the latent R-square, and flags
  for failed or non-converged runs, Mplus messages, nested models that
  fit better, ridge lines that drift out of the data, and `cd <= 0`.
  Mplus reports no CFI or SRMR with `XWITH`. All 16 models fit in about
  40 s at N = 400. `reference` runs the tests against another model.
- Four models beyond the RSA package (`addx2`, `addy2`, `noIA`,
  `a5zero`).
- `identification = "invariance"`: a shared X/Y metric through scalar
  invariance (shared labels, no `MODEL CONSTRAINT`), which runs under
  Bayes, with `free_intercepts` for partial invariance, and
  [`fit_rsa_mplus_invariance()`](https://franciscowilhelm.github.io/rsahelpers/reference/fit_rsa_mplus_invariance.md)
  to test the assumption.
- [`rsa_mplus_draws()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_draws.md):
  Monte Carlo intervals for the surface, ridge, and principal-axis
  parameters from the `TECH3` covariance matrix.
- [`rsa_mplus_supported_range()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_supported_range.md):
  the supported range of the LOC and LOIC from the model-implied latent
  distribution.

Next:

1.  Intervals for the derived parameters under Bayes, where the
    restricted models’ `MODEL CONSTRAINT` statements may not run (see
    Inference below).
2.  Plotting the observed region (the supported range, or points or a
    hull from `FSCORES`, which work with `XWITH`, or from scale means)
    and outlier screening (RSA’s Bollen-Jackman rule on scale means;
    Mplus Cook’s D and influence).
3.  A worked, executed example of the model family in the vignette
    (0.2.0 has the workflow with code only), and naming
    `model_type = "si_lms"` as RCSLMS (Cheung et al., 2021) in the
    documentation.
4.  Later: cubic models (CA, CL, RRCA, RRCL); an `XWITH` of a product
    works.

Not planned: the absolute-difference models (a step function of latent
variables cannot be expressed in LMS) and the weak and strong fit
patterns.

## Inference for derived parameters

Mplus gives delta-method standard errors for `NEW` parameters. These are
poor for ratios such as C, S, and the stationary point, where the `RSA`
package bootstraps. Mplus refuses `BOOTSTRAP` with MLR. With
`ESTIMATOR = ML` it works, but for the full latent model it takes about
0.7 s per draw at N = 300 (about 2 h for 10,000 draws per model). That
is too compute-intensive and is not pursued.

Fallbacks, in order of preference:

1.  **Monte Carlo CIs.** Done in
    [`rsa_mplus_draws()`](https://franciscowilhelm.github.io/rsahelpers/reference/rsa_mplus_draws.md):
    draw parameters from the asymptotic sampling distribution of the
    estimates (TECH3 covariance) and compute the derived parameters
    in R. No refitting is needed, so this is fast for any model.

2.  **Bayes.** `ESTIMATOR = BAYES` supports `XWITH`, and
    [`RSA_mplus()`](https://franciscowilhelm.github.io/rsahelpers/reference/RSA_mplus.md)
    already reads its credibility intervals. To check: whether the
    restricted models’ `MODEL CONSTRAINT` statements run under Bayes,
    and how models are compared without AICc.

3.  **Bootstrap with RCSLMS**, only if the draws can run in parallel.
    The reliability-corrected single-indicator LMS of Cheung et
    al. (2021) is `model_type = "si_lms"`, which always fixes the
    residual variances from the supplied reliabilities. Its measurement
    model is much smaller, but the bootstrap gains little. Timings on
    `models/latent-rsa.dat` (N = 300, Mplus 8.10):

    |                        | Latent | RCSLMS |
    |------------------------|--------|--------|
    | Free parameters        | 33     | 12     |
    | One MLR fit            | 2.0 s  | 0.4 s  |
    | ML bootstrap, per draw | 0.7 s  | 0.4 s  |

    The two-dimensional numerical integration over X and Y, which both
    share, dominates each draw, so 10,000 draws still take about 70 min
    per model. Mplus is not known to run bootstrap draws in parallel
    (not verified). Resampling in R with one Mplus process per draw
    would spread the draws over cores. Each process then starts Mplus
    once, but one RCSLMS fit including that start took 0.4 s. Eight
    processes might manage 10,000 draws in roughly 10 min; this is an
    untested estimate. The restricted models run with `si_lms`: SRRR
    gave S = 1.48, C = 0.22, and a4’ = -0.31, against 1.51, 0.23, and
    -0.31 in the latent model.
