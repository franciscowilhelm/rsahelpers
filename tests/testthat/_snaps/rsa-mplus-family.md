# models are ranked by AICc and tested against the full model

    Code
      print(comparison)
    Output
      <rsa_mplus_comparison>
       model  k       LL    AICc delta_AICc weight cum_weight evidence_ratio lrt_df
        SRRR 32 -3534.26 7138.27       0.00  0.668      0.668                     1
        full 33 -3533.77 7139.67       1.40  0.332      0.999           2.01     NA
         SRR 31 -3542.26 7151.91      13.63  0.001      1.000            913      2
       SRSQD 31 -3550.40 7168.18      29.91  0.000      1.000       3.12e+06      2
          IA 31 -3639.77 7346.94     208.66  0.000      1.000       2.04e+45      2
        null 28 -3670.07 7400.51     262.23  0.000      1.000       8.78e+56      5
       lrt_cd lrt_chisq  lrt_p    R2
        0.863      1.13  0.287 0.596
           NA        NA        0.600
        0.837     20.28 <0.001 0.585
        0.857     38.78 <0.001 0.551
       -0.745        NA        0.156
        0.789    345.62 <0.001 0.000
      
      Flags:
        IA: cd <= 0

# nesting violations and lines outside the data are flagged

    Code
      comparison <- rsa_mplus_comparison_table(fits, no_messages(), 0.01, 3)
    Condition
      Warning:
      A model fits better than a model it is nested in:
      x `SRSQD` fits better than `SRRR`, which nests it, by 0.050.
      i One of them did not reach its optimum, usually the larger model. Refit it, for example with start values from the smaller one or more random starts.

# the family runner continues after failures and starts from the full model

    Code
      family <- fit_rsa_mplus_family(dat, measurement, structural = "Z ~ X + Y + c*age",
        models = c("IA", "additive", "full", "SRR"), modelout = modelout,
        Mplus_command = "mplus")
    Condition
      Warning:
      2 of 4 models failed:
      x `additive` (run): Mplus reported an error.
      x `SRR` (create): Parameter label `c` is reserved for the RSA model.

# the family runner validates its arguments before writing

    Code
      family(models = c("full", "weak", "cubic"))
    Condition
      Error in `fit_rsa_mplus_family()`:
      ! Unknown models in `models`: `weak`, `cubic`.
      i The family has `full`, `SRRR`, `SRR`, `RR`, `SRSQD`, `SSQD`, `SQD`, `IA`, `additive`, `mean`, `diff`, `onlyx2`, `onlyy2`, `onlyx`, `onlyy`, `null`, `addx2`, `addy2`, `noIA`, `a5zero`.
    Code
      family(models = c("full", "full"))
    Condition
      Error in `fit_rsa_mplus_family()`:
      ! `models` must name each model once.
    Code
      family(models = "SRRR")
    Condition
      Error in `fit_rsa_mplus_family()`:
      ! `start = "full"` needs "full" in `models`.
    Code
      family(start = "SRRR")
    Condition
      Error in `fit_rsa_mplus_family()`:
      ! `start` must be "full", `NULL`, or start values for `create_rsa_mplus_model()`.
    Code
      fit_rsa_mplus_family(dat, measurement, modelout = "rsa")
    Condition
      Error in `fit_rsa_mplus_family()`:
      ! `modelout` must be a single path with an `.inp` extension.

