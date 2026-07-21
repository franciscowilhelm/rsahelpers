test_that("coverage uses shared bins and retains empty cells", {
  data <- data.frame(
    self = c(1, 1.5, 2, 3, 4, 5, NA),
    other = c(1, 2, 2.5, 2.5, 4.5, 5, 3)
  )
  result <- predictor_coverage(
    data,
    self,
    other,
    breaks = 4,
    discrepancy = c(0.5, 1),
    min_count = 2
  )

  expect_s3_class(result, "predictor_coverage")
  expect_equal(nrow(result$cells), 16L)
  expect_equal(sum(result$cells$count), 6L)
  expect_equal(result$summary$n, 6L)
  expect_equal(result$summary$n_missing, 1L)
  expect_equal(result$summary$x_greater_than_y, 1 / 6)
  expect_equal(result$summary$y_greater_than_x, 3 / 6)
  expect_equal(result$summary$ties, 2 / 6)
  expect_equal(nrow(result$discrepancy), 4L)
  expect_equal(result$breaks, seq(1, 5, length.out = 5))
})

test_that("explicit shared breaks include both endpoints", {
  data <- data.frame(x = c(0, 1, 2), y = c(2, 1, 0))
  result <- predictor_coverage(data, "x", "y", breaks = c(0, 1, 2))

  expect_equal(sum(result$cells$count), 3L)
  expect_equal(
    result$cells$count[
      result$cells$x_bin == 1 &
        result$cells$y_bin == 2
    ],
    1L
  )
  expect_equal(
    result$cells$count[
      result$cells$x_bin == 2 &
        result$cells$y_bin == 1
    ],
    1L
  )
})

test_that("predictor columns accept bare names and strings", {
  data <- data.frame(x = 1:3, y = 3:1)

  bare <- predictor_coverage(data, x, y, breaks = 2)
  strings <- predictor_coverage(data, "x", "y", breaks = 2)

  expect_equal(bare, strings)
})

test_that("coverage reports invalid inputs clearly", {
  expect_snapshot(
    error = TRUE,
    predictor_coverage(data.frame(x = letters[1:3], y = 1:3), x, y)
  )
  expect_snapshot(
    error = TRUE,
    predictor_coverage(data.frame(x = 1:3, y = 1:3), x, y, breaks = c(0, 2))
  )
  expect_snapshot(
    error = TRUE,
    predictor_coverage(data.frame(x = 1:3, y = 1:3), x, y, discrepancy = 0)
  )
})

test_that("coverage plots single and faceted results", {
  skip_if_not_installed("ggplot2")
  data <- data.frame(x = 1:6, y = c(1, 2, 2, 4, 5, 5))
  coverage <- predictor_coverage(data, x, y, breaks = 3)

  single <- plot_predictor_coverage(coverage)
  faceted <- plot_predictor_coverage(list(a = coverage, b = coverage))

  expect_s3_class(single, "ggplot")
  expect_s3_class(faceted, "ggplot")
  expect_s3_class(faceted$facet, "FacetWrap")
})
