# formula interface rejects unsupported specifications

    Code
      prepare_congruence_data(z ~ x + y, d)
    Condition
      Error:
      ! `formula` must have one response and two untransformed predictors: `z ~ x * y`.

---

    Code
      prepare_congruence_data(z ~ x * x, d)
    Condition
      Error:
      ! The Z, X, and Y variables in `formula` must be different columns.

---

    Code
      prepare_congruence_data(z ~ x * missing, d)
    Condition
      Error:
      ! Missing column(s): missing

---

    Code
      prepare_congruence_data(z ~ x * group, d)
    Condition
      Error:
      ! Spline variables must be numeric: group.

