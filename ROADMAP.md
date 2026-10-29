# Roadmap

## Latent RSA model family

Goal: the procedure of the `RSA` package and of Schönbrodt, *Testing fit
patterns with polynomial regression models*, as closely as possible in Mplus
(MLR with `XWITH`).

Done (development version):

- Effects coding and pooled centering of the latent X and Y factors.
- `RSA_mplus()` reads restricted models (`rsa_model`).
- `create_rsa_mplus_model(rsa_model =, start =)` generates the 16 second-order
  models, including SRRR and SRSQD with the ridge parameters C, S, bM, and a4'.

Next:

1. Batch runner over the family and a comparison table: AICc, Akaike weights,
   evidence ratios, Satorra-Bentler scaled LRT against the full model (from
   the log-likelihoods and MLR scaling factors), and checks that nested models
   never fit better. Mplus reports no CFI or SRMR with `XWITH`.
2. Confidence intervals for the derived parameters (see Inference below).
3. Plotting the observed region (points or hull from `FSCORES`, which work
   with `XWITH`, or from scale means) and outlier screening (RSA's
   Bollen-Jackman rule on scale means; Mplus Cook's D and influence).
4. Vignette section on the model family.
5. Later: cubic models (CA, CL, RRCA, RRCL); an `XWITH` of a product works.

Not planned: the absolute-difference models (a step function of latent
variables cannot be expressed in LMS) and the weak and strong fit patterns.

## Inference for derived parameters

Mplus gives delta-method standard errors for `NEW` parameters. These are poor
for ratios such as C, S, and the stationary point, where the `RSA` package
bootstraps. Mplus refuses `BOOTSTRAP` with MLR. With `ESTIMATOR = ML` it works,
but for the full latent model it takes about 0.7 s per draw at N = 300 (about
2 h for 10,000 draws per model). That is too compute-intensive and is not
pursued.

Fallbacks, from least to most effort:

1. **RCSLMS.** The reliability-corrected single-indicator LMS of Cheung et
   al. (2021) is `model_type = "si_lms"`, which always fixes the residual
   variances from the supplied reliabilities. Its measurement model is much
   smaller. Timings on `models/latent-rsa.dat` (N = 300, Mplus 8.10):

   | | Latent | RCSLMS |
   |---|---|---|
   | Free parameters | 33 | 12 |
   | One MLR fit | 2.0 s | 0.4 s |
   | ML bootstrap, per draw | 0.7 s | 0.4 s |

   The bootstrap gains less than a single fit because the two-dimensional
   numerical integration over X and Y, which both share, dominates each draw.
   At 0.4 s per draw, 10,000 draws take about 70 min. That is workable for a
   selected model but not for the whole family. The restricted models run
   with `si_lms`: SRRR gave S = 1.48, C = 0.22, and a4' = -0.31, against
   1.51, 0.23, and -0.31 in the latent model.

   To do:
   - An option to rerun a selected model with `ESTIMATOR = ML` and
     `BOOTSTRAP`, and read the bootstrap CIs of the `NEW` parameters.
   - Test the family with `si_lms`.
   - Try `PROCESSORS` and fewer integration points to cut the time per draw
     further; this is untested.
   - Name `si_lms` as RCSLMS (Cheung et al., 2021) in the documentation.
2. **Monte Carlo CIs.** Draw parameters from the asymptotic sampling
   distribution of the estimates (TECH3 covariance) and compute the derived
   parameters in R. No refitting is needed, so this is fast for any model.
3. **Bayes.** `ESTIMATOR = BAYES` supports `XWITH`, and `RSA_mplus()` already
   reads its credibility intervals. To check: whether the restricted models'
   `MODEL CONSTRAINT` statements run under Bayes, and how models are compared
   without AICc.
