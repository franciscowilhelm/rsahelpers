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

