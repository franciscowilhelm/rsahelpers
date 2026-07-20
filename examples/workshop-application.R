# Example application for the Edwards-Parry workshop dataset.
#
# Run from the package root:
#   Rscript examples/workshop-application.R

if (!requireNamespace("haven", quietly = TRUE)) {
  stop("Package 'haven' is required to read spline.sma/spline.dta.", call. = FALSE)
}

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(rsahelpers)
}

workshop <- haven::read_dta(
  system.file("extdata", "spline.dta",
    package = "rsahelpers",
    mustWork = TRUE
  )
)

fit_attribute <- function(data, prefix, label) {
  wanted <- paste0(prefix, "WC")
  actual <- paste0(prefix, "HC")

  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(label, "\n", sep = "")
  cat(strrep("=", 72), "\n\n", sep = "")

  piecewise <- fit_piecewise_congruence(
    data,
    wanted,
    actual,
    "JOBSAT",
    n_seams = 2,
    center = FALSE
  )
  cat("OLS comparison models\n")
  print(tidy_lm_summary(piecewise))

  one_seam <- fit_spline_congruence(
    data,
    wanted,
    actual,
    "JOBSAT",
    n_seams = 1,
    center = FALSE
  )
  cat("\nOne-seam nonlinear spline model\n")
  print(one_seam)

  cat("\nOne-seam surface features\n")
  print(surface_features(one_seam))

  cat("\nOne-seam joint tests\n")
  print(spline_tests(one_seam)$joint)

  two_seam <- fit_spline_congruence(
    data,
    wanted,
    actual,
    "JOBSAT",
    n_seams = 2,
    center = FALSE
  )
  cat("\nTwo-seam nonlinear spline model\n")
  print(two_seam)

  cat("\nModel comparison\n")
  comparison <- data.frame(
    model = c("one_seam", "two_seam"),
    rss = c(one_seam$rss, two_seam$rss),
    r.squared = c(one_seam$r.squared, two_seam$r.squared),
    df.residual = c(one_seam$df.residual, two_seam$df.residual)
  )
  print(comparison)

  output_dir <- file.path("examples", "output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  png_file <- file.path(output_dir, paste0(tolower(label), "-one-seam-surface.png"))
  grDevices::png(png_file, width = 1200, height = 900, res = 140)
  plot_spline_surface(
    one_seam,
    main = paste(label, "one-seam spline surface"),
    xlab = "Wanted amount (centered)",
    ylab = "Actual amount (centered)",
    zlab = "Job satisfaction"
  )
  grDevices::dev.off()
  cat("\nSaved 3D surface plot:", png_file, "\n", sep = "")

  invisible(list(
    piecewise = piecewise,
    one_seam = one_seam,
    two_seam = two_seam
  ))
}

authority <- fit_attribute(workshop, "ATH", "Authority")
variety <- fit_attribute(workshop, "VAR", "Variety")
