#' Build, write, run, and read RSA models for Mplus
#'
#' These functions provide a composable workflow for response surface models
#' estimated in Mplus. [create_rsa_mplus_model()] uses tidySEM for ordinary
#' measurement and structural paths, then adds the `XWITH` interactions and
#' response-surface constraints required for latent polynomial models.
#'
#' @param data A data frame containing all variables exported to Mplus.
#' @param measurement A named list mapping latent variables to indicators, or a
#'   character vector of lavaan/tidySEM measurement statements.
#' @param structural Optional structural model. Supply a named list mapping
#'   outcomes to predictors, or a character vector of lavaan/tidySEM paths.
#'   The core X and Y paths to the RSA outcome are added when absent.
#' @param roles Named character vector with elements `x`, `y`, and `outcome`
#'   identifying the latent variables used in the RSA model.
#' @param model_type Either `"latent"` for multi-indicator factors or
#'   `"si_lms"` for reliability-corrected single-indicator factors.
#' @param reliability For `model_type = "si_lms"`, a named numeric vector of
#'   reliability estimates for the three latent variables in `roles`.
#' @param constraints Which built-in constraints to include. `"surface"`
#'   produces `cs`, `cc`, `is`, `ic`, and `a5`; `"principal"` additionally
#'   produces the stationary point and principal axes.
#' @param blocks Named list of Mplus input blocks. Values replace defaults for
#'   the corresponding block. `MODEL` and `MODELCONSTRAINT` are reserved; use
#'   `model_extra` and `constraint_extra` to extend them.
#' @param model_extra Additional raw statements appended to the Mplus `MODEL`
#'   block.
#' @param constraint_extra Additional raw statements appended to the Mplus
#'   `MODEL CONSTRAINT` block.
#' @param usevariables Character vector of data columns to export. Defaults to
#'   all columns in `data`.
#' @param quiet If `TRUE`, suppress status messages from MplusAutomation.
#'
#' @return An object of class `rsa_mplus_workflow`. The object is updated and
#'   returned by the write, run, and read functions.
#' @export
#'
#' @examples
#' \dontrun{
#' latent_model <- create_rsa_mplus_model(
#'   data = dat,
#'   measurement = list(
#'     X = c("x1", "x2", "x3"),
#'     Y = c("y1", "y2", "y3"),
#'     Z = c("z1", "z2", "z3")
#'   ),
#'   structural = list(Z = c("X", "Y"))
#' )
#' }
create_rsa_mplus_model <- function(
  data,
  measurement,
  structural = NULL,
  roles = c(x = "X", y = "Y", outcome = "Z"),
  model_type = c("latent", "si_lms"),
  reliability = NULL,
  constraints = c("surface", "principal"),
  blocks = list(),
  model_extra = NULL,
  constraint_extra = NULL,
  usevariables = NULL,
  quiet = TRUE
) {
  require_rsa_mplus_namespace("tidySEM")
  require_rsa_mplus_namespace("MplusAutomation")

  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  if (anyDuplicated(names(data))) {
    rlang::abort("`data` must have unique column names.")
  }

  roles <- validate_rsa_mplus_roles(roles)
  model_type <- match.arg(model_type)
  constraints <- match.arg(
    constraints,
    choices = c("surface", "principal"),
    several.ok = TRUE
  )
  constraints <- unique(c("surface", constraints))

  if (is.null(usevariables)) {
    usevariables <- names(data)
  }
  validate_rsa_mplus_variables(data, usevariables)

  tidy_model <- tidySEM::tidy_sem(data[, usevariables, drop = FALSE])
  tidy_model <- add_rsa_mplus_measurement(tidy_model, measurement)
  tidy_model <- add_rsa_mplus_structural(tidy_model, structural, roles)

  syntax_table <- tidySEM::syntax(tidy_model)
  syntax_table <- label_rsa_mplus_core(syntax_table, roles)

  factor_names <- unique(syntax_table$lhs[syntax_table$op == "=~"])
  missing_roles <- setdiff(unname(roles), factor_names)
  if (length(missing_roles) > 0L) {
    rlang::abort(sprintf(
      "The measurement model does not define RSA factor%s %s.",
      if (length(missing_roles) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_roles), collapse = ", ")
    ))
  }

  reserved_factors <- intersect(c("XS", "XY", "YS"), factor_names)
  if (length(reserved_factors) > 0L) {
    rlang::abort(sprintf(
      "Measurement factor name%s %s reserved for RSA interactions.",
      if (length(reserved_factors) == 1L) "" else "s",
      paste(sprintf("`%s`", reserved_factors), collapse = ", ")
    ))
  }

  syntax_table <- apply_rsa_mplus_measurement_type(
    syntax_table = syntax_table,
    data = data,
    roles = roles,
    model_type = model_type,
    reliability = reliability
  )
  tidySEM::syntax(tidy_model) <- syntax_table

  conventional_model <- tidySEM::as_mplus(tidy_model)
  custom_model <- c(
    sprintf("XY | %s XWITH %s;", roles[["x"]], roles[["y"]]),
    sprintf("XS | %s XWITH %s;", roles[["x"]], roles[["x"]]),
    sprintf("YS | %s XWITH %s;", roles[["y"]], roles[["y"]]),
    sprintf("%s ON XS (b3);", roles[["outcome"]]),
    sprintf("%s ON XY (b4);", roles[["outcome"]]),
    sprintf("%s ON YS (b5);", roles[["outcome"]])
  )
  model_syntax <- paste(
    c(conventional_model, custom_model, flatten_rsa_mplus_text(model_extra)),
    collapse = "\n"
  )

  constraint_syntax <- rsa_mplus_constraints(constraints)
  constraint_syntax <- paste(
    c(constraint_syntax, flatten_rsa_mplus_text(constraint_extra)),
    collapse = "\n"
  )

  block_args <- validate_rsa_mplus_blocks(blocks)
  defaults <- list(
    TITLE = "RSA model generated by rsahelpers;",
    ANALYSIS = "TYPE = RANDOM;\nALGORITHM = INTEGRATION;",
    OUTPUT = "STDYX CINTERVAL TECH1;"
  )
  defaults[names(block_args)] <- block_args
  block_args <- defaults

  if (!grepl("TYPE\\s*=\\s*RANDOM", block_args$ANALYSIS, ignore.case = TRUE)) {
    rlang::abort(
      "The Mplus `ANALYSIS` block must contain `TYPE = RANDOM` for `XWITH`."
    )
  }

  mplus_args <- c(
    block_args,
    list(
      MODEL = model_syntax,
      MODELCONSTRAINT = constraint_syntax,
      usevariables = usevariables,
      rdata = data[, usevariables, drop = FALSE],
      autov = FALSE,
      quiet = quiet
    )
  )
  mplus_object <- do.call(MplusAutomation::mplusObject, mplus_args)

  metadata <- list(
    outcome = unname(roles[["outcome"]]),
    pred_x = unname(roles[["x"]]),
    pred_y = unname(roles[["y"]]),
    pred_x2 = "XS",
    pred_xy = "XY",
    pred_y2 = "YS",
    b0 = 0,
    coefficient_labels = c(x = "b1", y = "b2", x2 = "b3", xy = "b4", y2 = "b5")
  )

  new_rsa_mplus_workflow(
    spec = list(
      roles = roles,
      model_type = model_type,
      reliability = reliability,
      constraints = constraints,
      usevariables = usevariables,
      metadata = metadata
    ),
    tidysem = tidy_model,
    mplus = mplus_object,
    status = "created"
  )
}

#' @param x An `rsa_mplus_workflow` object.
#' @param modelout Path to the Mplus `.inp` file to create.
#' @param dataout Path to the Mplus data file. By default, it uses the same
#'   stem as `modelout` with a `.dat` extension.
#' @param overwrite If `TRUE`, replace existing input and data files.
#' @rdname create_rsa_mplus_model
#' @export
write_rsa_mplus_model <- function(
  x,
  modelout,
  dataout = NULL,
  overwrite = FALSE,
  quiet = TRUE
) {
  require_rsa_mplus_namespace("MplusAutomation")
  validate_rsa_mplus_workflow(x)

  modelout <- normalizePath(modelout, mustWork = FALSE)
  if (!grepl("\\.inp$", modelout, ignore.case = TRUE)) {
    rlang::abort("`modelout` must have an `.inp` extension.")
  }
  if (is.null(dataout)) {
    dataout <- sub("\\.[^.]*$", ".dat", modelout)
    if (identical(dataout, modelout)) {
      dataout <- paste0(modelout, ".dat")
    }
  }
  dataout <- normalizePath(dataout, mustWork = FALSE)

  parent_dirs <- unique(dirname(c(modelout, dataout)))
  missing_dirs <- parent_dirs[!dir.exists(parent_dirs)]
  if (length(missing_dirs) > 0L) {
    rlang::abort(sprintf(
      "Output director%s %s not exist: %s.",
      if (length(missing_dirs) == 1L) "y" else "ies",
      if (length(missing_dirs) == 1L) "does" else "do",
      paste(sprintf("`%s`", missing_dirs), collapse = ", ")
    ))
  }

  existing <- c(modelout, dataout)[file.exists(c(modelout, dataout))]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    rlang::abort(sprintf(
      "Refusing to overwrite existing file%s: %s. Use `overwrite = TRUE` to replace %s.",
      if (length(existing) == 1L) "" else "s",
      paste(sprintf("`%s`", existing), collapse = ", "),
      if (length(existing) == 1L) "it" else "them"
    ))
  }

  modeler_args <- list(
    object = x$mplus,
    dataout = dataout,
    modelout = modelout,
    run = 0L,
    check = FALSE,
    writeData = "always",
    hashfilename = FALSE,
    quiet = quiet
  )
  modeled <- if (isTRUE(overwrite)) {
    invisible(utils::capture.output(
      result <- suppressMessages(suppressWarnings(
        do.call(MplusAutomation::mplusModeler, modeler_args)
      ))
    ))
    result
  } else {
    do.call(MplusAutomation::mplusModeler, modeler_args)
  }

  x$mplus <- modeled
  x$files <- list(
    input = modelout,
    data = dataout,
    output = sub("\\.inp$", ".out", modelout, ignore.case = TRUE)
  )
  x$status <- "written"
  x
}

#' @param Mplus_command Optional path or command used to invoke Mplus.
#' @param replace_outfile Passed to [MplusAutomation::runModels()].
#' @rdname create_rsa_mplus_model
#' @export
run_rsa_mplus_model <- function(
  x,
  Mplus_command = NULL,
  replace_outfile = c("always", "never", "modifiedDate"),
  quiet = TRUE
) {
  require_rsa_mplus_namespace("MplusAutomation")
  replace_outfile <- match.arg(replace_outfile)

  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    input_path <- normalizePath(x, mustWork = TRUE)
    if (!grepl("\\.inp$", input_path, ignore.case = TRUE)) {
      rlang::abort("A character `x` must point to an Mplus `.inp` file.")
    }
    x <- new_rsa_mplus_workflow(
      files = list(
        input = input_path,
        data = NULL,
        output = sub("\\.inp$", ".out", input_path, ignore.case = TRUE)
      ),
      status = "written"
    )
  } else {
    validate_rsa_mplus_workflow(x)
  }

  if (is.null(x$files$input) || !file.exists(x$files$input)) {
    rlang::abort("The workflow must be written before it can be run.")
  }

  command <- resolve_rsa_mplus_command(Mplus_command)
  input_path <- x$files$input

  MplusAutomation::runModels(
    target = input_path,
    recursive = FALSE,
    filefilter = NULL,
    showOutput = FALSE,
    replaceOutfile = replace_outfile,
    logFile = NULL,
    Mplus_command = command,
    quiet = quiet
  )

  output_path <- sub("\\.inp$", ".out", input_path, ignore.case = TRUE)
  if (!file.exists(output_path)) {
    rlang::abort(sprintf(
      "Mplus did not create the expected output file `%s`.",
      output_path
    ))
  }

  x$files$output <- output_path
  x <- read_rsa_mplus_model(x, quiet = quiet)

  model_errors <- x$results$errors
  if (!is.null(model_errors) && length(model_errors) > 0L) {
    x$status <- "failed"
    rlang::abort(sprintf(
      "Mplus reported an error. Output was preserved at `%s`: %s",
      output_path,
      paste(model_errors, collapse = " ")
    ))
  }

  x$status <- "run"
  x
}

#' @param target An `rsa_mplus_workflow`, Mplus `.inp` or `.out` path, an
#'   `mplus.model`, or an `mplusObject` containing parsed results.
#' @rdname create_rsa_mplus_model
#' @export
read_rsa_mplus_model <- function(target, quiet = TRUE) {
  require_rsa_mplus_namespace("MplusAutomation")

  if (inherits(target, "rsa_mplus_workflow")) {
    x <- target
    output_path <- x$files$output
    if (is.null(output_path) && !is.null(x$files$input)) {
      output_path <- sub("\\.inp$", ".out", x$files$input, ignore.case = TRUE)
    }
    if (is.null(output_path)) {
      rlang::abort("The workflow does not contain an Mplus output path.")
    }
    output_path <- normalizePath(output_path, mustWork = TRUE)
    x$results <- MplusAutomation::readModels(output_path, quiet = quiet)
    x$files$output <- output_path
    x$status <- "read"
    return(x)
  }

  if (inherits(target, "mplus.model")) {
    return(new_rsa_mplus_workflow(results = target, status = "read"))
  }

  if (inherits(target, "mplusObject")) {
    if (is.null(target$results)) {
      rlang::abort("The `mplusObject` does not contain parsed results.")
    }
    return(new_rsa_mplus_workflow(
      mplus = target,
      results = target$results,
      status = "read"
    ))
  }

  if (!is.character(target) || length(target) != 1L || is.na(target)) {
    rlang::abort(
      "`target` must be a workflow, Mplus object, or single file path."
    )
  }

  target <- normalizePath(target, mustWork = TRUE)
  output_path <- if (grepl("\\.inp$", target, ignore.case = TRUE)) {
    sub("\\.inp$", ".out", target, ignore.case = TRUE)
  } else {
    target
  }
  output_path <- normalizePath(output_path, mustWork = TRUE)

  new_rsa_mplus_workflow(
    files = list(input = NULL, data = NULL, output = output_path),
    results = MplusAutomation::readModels(output_path, quiet = quiet),
    status = "read"
  )
}

#' @param run If `TRUE`, run and read the generated Mplus model after writing
#'   its files. The default only writes reproducible input and data files.
#' @rdname create_rsa_mplus_model
#' @export
rsa_mplus_workflow <- function(
  data,
  measurement,
  structural = NULL,
  roles = c(x = "X", y = "Y", outcome = "Z"),
  model_type = c("latent", "si_lms"),
  reliability = NULL,
  constraints = c("surface", "principal"),
  blocks = list(),
  model_extra = NULL,
  constraint_extra = NULL,
  usevariables = NULL,
  modelout,
  dataout = NULL,
  overwrite = FALSE,
  run = FALSE,
  Mplus_command = NULL,
  quiet = TRUE
) {
  out <- create_rsa_mplus_model(
    data = data,
    measurement = measurement,
    structural = structural,
    roles = roles,
    model_type = model_type,
    reliability = reliability,
    constraints = constraints,
    blocks = blocks,
    model_extra = model_extra,
    constraint_extra = constraint_extra,
    usevariables = usevariables,
    quiet = quiet
  )
  out <- write_rsa_mplus_model(
    out,
    modelout = modelout,
    dataout = dataout,
    overwrite = overwrite,
    quiet = quiet
  )

  if (isTRUE(run)) {
    out <- run_rsa_mplus_model(
      out,
      Mplus_command = Mplus_command,
      quiet = quiet
    )
  }

  out
}

#' @export
print.rsa_mplus_workflow <- function(x, ...) {
  cat("<rsa_mplus_workflow>\n")
  cat("Status:", x$status, "\n")
  if (!is.null(x$spec$model_type)) {
    cat("Model type:", x$spec$model_type, "\n")
  }
  if (!is.null(x$files$input)) {
    cat("Input:", x$files$input, "\n")
  }
  if (!is.null(x$files$output) && file.exists(x$files$output)) {
    cat("Output:", x$files$output, "\n")
  }
  invisible(x)
}

new_rsa_mplus_workflow <- function(
  spec = NULL,
  tidysem = NULL,
  mplus = NULL,
  files = list(input = NULL, data = NULL, output = NULL),
  results = NULL,
  status = NULL
) {
  structure(
    list(
      spec = spec,
      tidysem = tidysem,
      mplus = mplus,
      files = files,
      results = results,
      status = status
    ),
    class = "rsa_mplus_workflow"
  )
}

validate_rsa_mplus_workflow <- function(x) {
  if (!inherits(x, "rsa_mplus_workflow")) {
    rlang::abort("`x` must be an `rsa_mplus_workflow` object.")
  }
  invisible(x)
}

require_rsa_mplus_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    rlang::abort(sprintf(
      "Package `%s` must be installed for the RSA-Mplus workflow.",
      package
    ))
  }
}

validate_rsa_mplus_roles <- function(roles) {
  required <- c("x", "y", "outcome")
  if (
    !is.character(roles) ||
      is.null(names(roles)) ||
      !all(required %in% names(roles))
  ) {
    rlang::abort(
      "`roles` must be a named character vector with `x`, `y`, and `outcome`."
    )
  }
  roles <- roles[required]
  if (anyNA(roles) || any(!nzchar(roles)) || anyDuplicated(toupper(roles))) {
    rlang::abort(
      "The three values in `roles` must be non-missing and distinct."
    )
  }
  roles
}

validate_rsa_mplus_variables <- function(data, usevariables) {
  if (
    !is.character(usevariables) ||
      anyNA(usevariables) ||
      any(!nzchar(usevariables))
  ) {
    rlang::abort("`usevariables` must be a non-missing character vector.")
  }
  missing_variables <- setdiff(usevariables, names(data))
  if (length(missing_variables) > 0L) {
    rlang::abort(sprintf(
      "Variable%s %s not found in `data`.",
      if (length(missing_variables) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_variables), collapse = ", ")
    ))
  }
  invalid <- usevariables[!grepl("^[A-Za-z][A-Za-z0-9_]*$", usevariables)]
  if (length(invalid) > 0L) {
    rlang::abort(sprintf(
      "Mplus variable name%s %s invalid.",
      if (length(invalid) == 1L) "" else "s",
      paste(sprintf("`%s`", invalid), collapse = ", ")
    ))
  }
  shortened <- toupper(substr(usevariables, 1L, 8L))
  collisions <- unique(shortened[duplicated(shortened)])
  if (length(collisions) > 0L) {
    rlang::abort(sprintf(
      "Variable names must be unique in their first eight characters for Mplus output; collision%s: %s.",
      if (length(collisions) == 1L) "" else "s",
      paste(collisions, collapse = ", ")
    ))
  }
  invisible(usevariables)
}

add_rsa_mplus_measurement <- function(tidy_model, measurement) {
  if (is.list(measurement) && !is.data.frame(measurement)) {
    if (is.null(names(measurement)) || any(!nzchar(names(measurement)))) {
      rlang::abort("A list `measurement` must be named by latent variable.")
    }
    if (anyDuplicated(toupper(names(measurement)))) {
      rlang::abort("Latent variable names in `measurement` must be unique.")
    }
    invalid <- !vapply(
      measurement,
      function(x) {
        is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
      },
      logical(1)
    )
    if (any(invalid)) {
      rlang::abort(
        "Every `measurement` list element must contain indicator names."
      )
    }

    indicators <- unlist(measurement, use.names = FALSE)
    if (anyDuplicated(indicators)) {
      rlang::abort(
        "Each observed indicator may belong to only one measurement factor."
      )
    }
    dictionary <- tidySEM::dictionary(tidy_model)
    missing_indicators <- setdiff(indicators, dictionary$name)
    if (length(missing_indicators) > 0L) {
      rlang::abort(sprintf(
        "Measurement indicator%s %s not found in exported data.",
        if (length(missing_indicators) == 1L) "" else "s",
        paste(sprintf("`%s`", missing_indicators), collapse = ", ")
      ))
    }
    dictionary$scale <- NA_character_
    for (factor_name in names(measurement)) {
      dictionary$scale[
        dictionary$name %in% measurement[[factor_name]]
      ] <- factor_name
    }
    tidySEM::dictionary(tidy_model) <- dictionary
    return(tidySEM::measurement(tidy_model))
  }

  statements <- flatten_rsa_mplus_text(measurement)
  if (length(statements) == 0L) {
    rlang::abort(
      "`measurement` must contain at least one measurement statement."
    )
  }
  do.call(tidySEM::add_paths, c(list(model = tidy_model), as.list(statements)))
}

add_rsa_mplus_structural <- function(tidy_model, structural, roles) {
  if (is.null(structural)) {
    statements <- character()
  } else if (is.list(structural) && !is.data.frame(structural)) {
    if (is.null(names(structural)) || any(!nzchar(names(structural)))) {
      rlang::abort("A list `structural` must be named by outcome variable.")
    }
    statements <- vapply(
      names(structural),
      function(outcome) {
        predictors <- structural[[outcome]]
        if (
          !is.character(predictors) ||
            length(predictors) == 0L ||
            anyNA(predictors)
        ) {
          rlang::abort(
            "Every `structural` list element must contain predictor names."
          )
        }
        sprintf("%s ~ %s", outcome, paste(predictors, collapse = " + "))
      },
      character(1)
    )
  } else {
    statements <- flatten_rsa_mplus_text(structural)
  }

  if (length(statements) > 0L) {
    tidy_model <- do.call(
      tidySEM::add_paths,
      c(list(model = tidy_model), as.list(statements))
    )
  }

  syntax_table <- tidySEM::syntax(tidy_model)
  for (predictor in roles[c("x", "y")]) {
    present <- any(
      syntax_table$op == "~" &
        syntax_table$lhs == roles[["outcome"]] &
        syntax_table$rhs == predictor
    )
    if (!present) {
      tidy_model <- tidySEM::add_paths(
        tidy_model,
        sprintf("%s ~ %s", roles[["outcome"]], predictor)
      )
      syntax_table <- tidySEM::syntax(tidy_model)
    }
  }
  tidy_model
}

label_rsa_mplus_core <- function(syntax_table, roles) {
  core <- list(b1 = roles[["x"]], b2 = roles[["y"]])
  reserved <- c("b1", "b2", "b3", "b4", "b5")
  labels <- tolower(syntax_table$label)

  for (label in names(core)) {
    hits <- which(
      syntax_table$op == "~" &
        syntax_table$lhs == roles[["outcome"]] &
        syntax_table$rhs == core[[label]]
    )
    if (length(hits) != 1L) {
      rlang::abort(sprintf(
        "Expected exactly one `%s ~ %s` structural path.",
        roles[["outcome"]],
        core[[label]]
      ))
    }
    conflicts <- which(labels == label & seq_len(nrow(syntax_table)) != hits)
    if (length(conflicts) > 0L) {
      rlang::abort(sprintf(
        "Parameter label `%s` is reserved for the RSA model.",
        label
      ))
    }
    syntax_table$label[hits] <- label
  }

  unused_conflicts <- setdiff(reserved, names(core))
  if (any(tolower(syntax_table$label) %in% unused_conflicts)) {
    conflict <- syntax_table$label[
      tolower(syntax_table$label) %in% unused_conflicts
    ][[1]]
    rlang::abort(sprintf(
      "Parameter label `%s` is reserved for the RSA model.",
      conflict
    ))
  }
  syntax_table
}

apply_rsa_mplus_measurement_type <- function(
  syntax_table,
  data,
  roles,
  model_type,
  reliability
) {
  counts <- vapply(
    unname(roles),
    function(factor_name) {
      sum(syntax_table$op == "=~" & syntax_table$lhs == factor_name)
    },
    integer(1)
  )

  if (identical(model_type, "latent")) {
    if (any(counts < 2L)) {
      invalid <- unname(roles)[counts < 2L]
      rlang::abort(sprintf(
        "Latent RSA factor%s %s must have at least two indicators; use `model_type = \"si_lms\"` for single indicators.",
        if (length(invalid) == 1L) "" else "s",
        paste(sprintf("`%s`", invalid), collapse = ", ")
      ))
    }
    return(syntax_table)
  }

  if (any(counts != 1L)) {
    invalid <- unname(roles)[counts != 1L]
    rlang::abort(sprintf(
      "SI-LMS factor%s %s must each have exactly one indicator.",
      if (length(invalid) == 1L) "" else "s",
      paste(sprintf("`%s`", invalid), collapse = ", ")
    ))
  }
  if (!is.numeric(reliability) || is.null(names(reliability))) {
    rlang::abort(
      "`reliability` must be a named numeric vector for SI-LMS models."
    )
  }
  missing_reliability <- setdiff(unname(roles), names(reliability))
  if (length(missing_reliability) > 0L) {
    rlang::abort(sprintf(
      "Missing reliability estimate%s for %s.",
      if (length(missing_reliability) == 1L) "" else "s",
      paste(sprintf("`%s`", missing_reliability), collapse = ", ")
    ))
  }
  selected_reliability <- reliability[unname(roles)]
  if (
    anyNA(selected_reliability) ||
      any(selected_reliability <= 0 | selected_reliability > 1)
  ) {
    rlang::abort(
      "SI-LMS reliability estimates must be greater than 0 and no greater than 1."
    )
  }

  for (factor_name in unname(roles)) {
    loading_row <- which(
      syntax_table$op == "=~" & syntax_table$lhs == factor_name
    )
    indicator <- syntax_table$rhs[loading_row]
    values <- data[[indicator]]
    if (!is.numeric(values)) {
      rlang::abort(sprintf("SI-LMS score `%s` must be numeric.", indicator))
    }
    observed_variance <- stats::var(values, na.rm = TRUE)
    if (is.na(observed_variance)) {
      rlang::abort(sprintf(
        "SI-LMS score `%s` needs at least two observed values.",
        indicator
      ))
    }
    residual_variance <- (1 - reliability[[factor_name]]) * observed_variance
    residual_row <- which(
      syntax_table$op == "~~" &
        syntax_table$lhs == indicator &
        syntax_table$rhs == indicator
    )
    if (length(residual_row) != 1L) {
      rlang::abort(sprintf(
        "Could not identify the residual variance for `%s`.",
        indicator
      ))
    }
    syntax_table$free[residual_row] <- 0
    syntax_table$ustart[residual_row] <- residual_variance
  }
  syntax_table
}

rsa_mplus_constraints <- function(constraints) {
  surface <- c(
    "NEW(cs cc is ic a5);",
    "cs = b1 + b2;",
    "cc = b3 + b4 + b5;",
    "is = b1 - b2;",
    "ic = b3 - b4 + b5;",
    "a5 = b3 - b5;"
  )
  if (!"principal" %in% constraints) {
    return(surface)
  }
  c(
    "NEW(cs cc is ic a5 x0 y0 p10 p11 p20 p21);",
    surface[-1],
    "x0 = (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4*b4);",
    "y0 = (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4*b4);",
    "p11 = (b5 - b3 + SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;",
    "p21 = (b5 - b3 - SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4;",
    "p10 = y0 - p11*x0;",
    "p20 = y0 - p21*x0;"
  )
}

validate_rsa_mplus_blocks <- function(blocks) {
  if (!is.list(blocks) || (length(blocks) > 0L && is.null(names(blocks)))) {
    rlang::abort("`blocks` must be a named list of Mplus input blocks.")
  }
  if (length(blocks) == 0L) {
    return(blocks)
  }
  names(blocks) <- toupper(names(blocks))
  allowed <- c(
    "TITLE",
    "DATA",
    "VARIABLE",
    "DEFINE",
    "MONTECARLO",
    "MODELPOPULATION",
    "MODELMISSING",
    "ANALYSIS",
    "MODELINDIRECT",
    "MODELTEST",
    "MODELPRIORS",
    "OUTPUT",
    "SAVEDATA",
    "PLOT"
  )
  invalid <- setdiff(names(blocks), allowed)
  if (length(invalid) > 0L) {
    rlang::abort(sprintf(
      "Unsupported or reserved Mplus block%s: %s.",
      if (length(invalid) == 1L) "" else "s",
      paste(sprintf("`%s`", invalid), collapse = ", ")
    ))
  }
  blocks <- lapply(blocks, function(x) {
    paste(flatten_rsa_mplus_text(x), collapse = "\n")
  })
  blocks
}

flatten_rsa_mplus_text <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }
  if (!is.character(x) || anyNA(x)) {
    rlang::abort(
      "Mplus syntax inputs must contain non-missing character values."
    )
  }
  x[nzchar(trimws(x))]
}

resolve_rsa_mplus_command <- function(command = NULL) {
  if (!is.null(command)) {
    if (!is.character(command) || length(command) != 1L || is.na(command)) {
      rlang::abort(
        "`Mplus_command` must be a single command or executable path."
      )
    }
    return(command)
  }

  candidates <- Sys.which(c("mplus", "Mplus"))
  candidates <- unname(candidates[nzchar(candidates)])
  if (length(candidates) > 0L) {
    return(candidates[[1]])
  }

  detected <- tryCatch(
    MplusAutomation::detectMplus(),
    error = function(e) ""
  )
  if (!is.character(detected) || length(detected) != 1L || !nzchar(detected)) {
    rlang::abort("Could not locate Mplus; supply `Mplus_command` explicitly.")
  }
  detected
}
