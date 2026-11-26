# RSA_mplus extracts unstandardized coefficients from Mplus output

    Code
      result <- RSA_mplus(model = model_path, outcome = "Z", pred_x = "X", pred_y = "Y",
        pred_x2 = "XS", pred_xy = "XY", pred_y2 = "YS", new_labels = c("CS", "CC",
          "IS", "IC", "A5"), plot = FALSE)
    Condition
      Warning:
      No intercept for outcome `Z` was found in the Mplus expectation parameters; using `b0 = 0`.

# RSA_mplus accepts an mplus.model object and standardized coefficients

    Code
      result <- RSA_mplus(model = mplus_model, outcome = "Z", pred_x = "X", pred_y = "Y",
        pred_x2 = "XS", pred_xy = "XY", pred_y2 = "YS", coef_type = "stdyx",
        include_new = FALSE, plot = FALSE)
    Condition
      Warning:
      No intercept for outcome `Z` was found in the Mplus expectation parameters; using `b0 = 0`.

# rsa_model reads terms fixed at zero and rejects estimated ones

    Code
      RSA_mplus(find_rsa_mplus_example(), outcome = "Z", pred_x = "X", pred_y = "Y",
      pred_x2 = "XS", pred_xy = "XY", pred_y2 = "YS", rsa_model = "IA", b0 = 0, plot = FALSE)
    Condition
      Error in `RSA_mplus()`:
      ! `rsa_model = "IA"` sets `x2`, `y2` to zero, but the Mplus output has nonzero estimates:
      x `Z<-XS` = -0.058
      x `Z<-YS` = -0.085

# RSA_mplus requires the labels of the terms the model estimates

    Code
      RSA_mplus(model_path, outcome = "Z", pred_x = "X", pred_y = "Y")
    Condition
      Error in `RSA_mplus()`:
      ! Missing required arguments: `pred_x2`, `pred_xy`, `pred_y2`.
      i For a surface without some terms, name it with `rsa_model`, e.g. `rsa_model = "additive"`.

