# coverage reports invalid inputs clearly

    Code
      predictor_coverage(data.frame(x = letters[1:3], y = 1:3), x, y)
    Condition
      Error in `predictor_coverage()`:
      ! `x` and `y` must select numeric columns.

---

    Code
      predictor_coverage(data.frame(x = 1:3, y = 1:3), x, y, breaks = c(0, 2))
    Condition
      Error in `coverage_breaks()`:
      ! Explicit `breaks` must span all finite X and Y values.

---

    Code
      predictor_coverage(data.frame(x = 1:3, y = 1:3), x, y, discrepancy = 0)
    Condition
      Error in `validate_discrepancy_cutoffs()`:
      ! `discrepancy` must contain finite positive numbers.

