#' rsahelpers: Helpers for Response Surface Analysis
#'
#' `rsahelpers` provides tools for fitting and inspecting response surfaces,
#' including Edwards-Parry spline regression models and Mplus coefficient
#' extraction and plotting.
#'
#' @keywords internal
#' @importFrom stats coef logLik nobs vcov
"_PACKAGE"

# Column names referenced via non-standard evaluation in ggplot2 aes().
utils::globalVariables(c("x", "y", "z"))
