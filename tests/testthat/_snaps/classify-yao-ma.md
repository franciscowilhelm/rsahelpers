# classification validates parameter tables and margins

    Code
      classify_yao_ma(data.frame(parameter = "CS"), valence = "positive")
    Condition
      Error in `normalize_yao_ma_table()`:
      ! Parameter data is missing columns: `estimate`, `conf.low`, `conf.high`.

---

    Code
      classify_yao_ma(yao_ma_test_parameters(), valence = "positive", equivalence = c(
        CS = 0.1))
    Condition
      Error in `validate_yao_ma_margins()`:
      ! `equivalence` must name exactly CS, CC, IS, IC, axis_intercept, and axis_slope (missing CC, IS, IC, axis_intercept, axis_slope).

