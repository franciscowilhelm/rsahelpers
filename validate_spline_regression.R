# Developer validation entrypoint.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required to run package validation.", call. = FALSE)
}

devtools::test()
