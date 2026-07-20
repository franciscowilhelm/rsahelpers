# latent RSA-Mplus models combine tidySEM and XWITH syntax

    Code
      cat(model$mplus$MODEL, "\nMODEL CONSTRAINT:\n", model$mplus$MODELCONSTRAINT,
      "\n")
    Output
      X BY x1@1;
      X BY x2;
      X BY x3;
      Y BY y1@1;
      Y BY y2;
      Y BY y3;
      Z BY z1@1;
      Z BY z2;
      Z BY z3;
      x1;
      x2;
      x3;
      y1;
      y2;
      y3;
      z1;
      z2;
      z3;
      X;
      Y;
      Z;
      X WITH Y;
      [x1];
      [x2];
      [x3];
      [y1];
      [y2];
      [y3];
      [z1];
      [z2];
      [z3];
      [X@0];
      [Y@0];
      [Z@0];
      Z ON X (b1);
      Z ON Y (b2);
      XY | X XWITH Y;
      XS | X XWITH X;
      YS | Y XWITH Y;
      Z ON XS (b3);
      Z ON XY (b4);
      Z ON YS (b5); 
      MODEL CONSTRAINT:
       NEW(cs cc is ic a5 x0 y0 p10 p11 p20 p21);
      cs = b1 + b2;
      cc = b3 + b4 + b5;
      is = b1 - b2;
      ic = b3 - b4 + b5;
      a5 = b3 - b5;
      x0 = (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4*b4);
      y0 = (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4*b4);
      p11 = (b5 - b3 + SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;
      p21 = (b5 - b3 - SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;
      p10 = y0 - p11*x0;
      p20 = y0 - p21*x0; 

# RSA-Mplus validation reports actionable model errors

    Code
      create_rsa_mplus_model(dat, list(X = "x", Y = "y", Z = "z"))
    Condition
      Error in `apply_rsa_mplus_measurement_type()`:
      ! Latent RSA factors `X`, `Y`, `Z` must have at least two indicators; use `model_type = "si_lms"` for single indicators.

# RSA-Mplus validation catches Mplus name collisions

    Code
      create_rsa_mplus_model(dat, list(X = c("longname_a", "longname_b"), Y = c("y1",
        "y2"), Z = c("z1", "z2")))
    Condition
      Error in `validate_rsa_mplus_variables()`:
      ! Variable names must be unique in their first eight characters for Mplus output; collision: LONGNAME.

# writing creates reproducible input and data files safely

    Code
      write_rsa_mplus_model(model, modelout = input_path)
    Condition
      Error in `write_rsa_mplus_model()`:
      ! Refusing to overwrite existing files: `<file>`, `<file>`. Use `overwrite = TRUE` to replace them.

