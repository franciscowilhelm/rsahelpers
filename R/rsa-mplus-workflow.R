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
#'   produces the stationary point and principal axes for the models that
#'   have them (see Details).
#' @param blocks Named list of Mplus input blocks. Values replace defaults for
#'   the corresponding block. `MODEL` and `MODELCONSTRAINT` are reserved; use
#'   `model_extra` and `constraint_extra` to extend them. The `VARIABLE` block
#'   may contain additional directives, but `USEVARIABLES` is generated from
#'   `usevariables` and must not be supplied through `blocks`.
#' @param model_extra Additional raw statements appended to the Mplus `MODEL`
#'   block.
#' @param constraint_extra Additional raw statements appended to the Mplus
#'   `MODEL CONSTRAINT` block.
#' @param usevariables Character vector of data columns to export and include
#'   in the Mplus `USEVARIABLES` statement. Defaults to all columns in `data`.
#' @param identification How the X, Y, and outcome factors are scaled, for
#'   `model_type = "latent"`. `"marker"` (default) fixes the first loading of
#'   each factor to 1. `"effects"` uses effects coding (Little, Slegers, &
#'   Card, 2006): loadings are constrained to average 1, so each factor is on
#'   the metric of its average indicator. This changes the units of b1-b5
#'   slightly but not the fit of the model. `"invariance"` keeps the marker
#'   loadings and puts the X and Y factors on one metric through scalar
#'   measurement invariance: the k-th X and Y indicators share their loading
#'   and their intercept, the two marker intercepts are fixed at 0, and both
#'   factor means are freed. It requires as many X as Y indicators, given in
#'   matching order. Unlike `"effects"`, it restricts the model, and it writes
#'   no `MODEL CONSTRAINT` statements for the measurement model, so it also
#'   runs under `ESTIMATOR = BAYES`. See Details.
#' @param center Origin of the X/Y plane. `"none"` (default) fixes both factor
#'   means at 0, so each factor is centred at its own mean and the line X = Y
#'   means "both at their average" rather than equality. `"pooled"`, which
#'   requires `identification = "effects"` or `"invariance"`, subtracts one
#'   constant (the mean of the X and Y scale means) from every X and Y
#'   indicator in a `DEFINE` statement and frees the two factor means; under
#'   effects coding it also constrains the intercepts to sum to zero. X = Y
#'   is then equality on the shared response scale, the
#'   latent counterpart of `RSA::RSA(center = "pooled")`. This assumes the X
#'   and Y items share a response metric. The move is a translation by the
#'   latent mean difference Delta = E(X) - E(Y): b3-b5, a2, a4, and a5 are
#'   unchanged, while a3 shifts by `-Delta * a4` and a1 by `-Delta * a5`.
#'   User `DEFINE` statements are placed after the shift.
#' @param free_intercepts For `identification = "invariance"`, names of
#'   indicators whose intercept is not shared between X and Y (partial scalar
#'   invariance). Naming either member of a matched pair frees the pair. The
#'   marker pair cannot be freed.
#' @param rsa_model Response surface model to generate, named as in the model
#'   family of the `RSA` package. `"full"` (the default) estimates all five
#'   polynomial terms. The restricted models are `"SRRR"`, `"SRR"`, `"RR"`,
#'   `"SRSQD"`, `"SSQD"`, `"SQD"`, `"IA"`, `"additive"`, `"mean"`, `"diff"`,
#'   `"onlyx"`, `"onlyy"`, `"onlyx2"`, `"onlyy2"`, and `"null"`. Four further
#'   models are not part of the `RSA` package: `"addx2"`, `"addy2"`, `"noIA"`,
#'   and `"a5zero"`. See Details. [fit_rsa_mplus_family()] fits several of them from one
#'   specification.
#' @param start Optional start values for the polynomial terms, typically the
#'   estimates of the full model. Supply a named numeric vector with names
#'   `b1` to `b5` (or `x`, `y`, `x2`, `xy`, and `y2`), a result of
#'   [RSA_mplus()], or a fitted workflow. Only the terms that `rsa_model`
#'   estimates freely receive a start value.
#' @param quiet If `TRUE`, suppress status messages from MplusAutomation.
#'
#' @details
#' `rsa_model` writes the restrictions of the second-order models in
#' `RSA::RSA()` (Schönbrodt, *Testing fit patterns with polynomial regression
#' models*) on the regression weights `b1` to `b5`. Terms set to zero are fixed
#' with `@0`. The other restrictions are placed in `MODEL CONSTRAINT`, together
#' with the `NEW` parameters that describe the restricted surface:
#'
#' | `rsa_model` | Restrictions | Model-specific `NEW` parameters |
#' |---|---|---|
#' | `"full"` | none | `x0`, `y0`, `p10`, `p11`, `p20`, `p21` (with `"principal"`) |
#' | `"SRRR"` | `b3 = b4*b4/(4*b5)` | `c`, `s`, `bm`, `a4r` |
#' | `"SRR"` | `b3 = b5`, `b4 = -2*b5` | `c`, `bm` |
#' | `"RR"` | `b1 = b2`, `b3 = b5`, `b4 = -2*b5` | `bm` |
#' | `"SRSQD"` | `b1 = b2*b4/(2*b5)`, `b3 = b4*b4/(4*b5)` | `c`, `s`, `a4r` |
#' | `"SSQD"` | `b1 = -b2`, `b3 = b5`, `b4 = -2*b5` | `c` |
#' | `"SQD"` | `b1 = b2 = 0`, `b3 = b5`, `b4 = -2*b5` | none |
#' | `"IA"` | `b3 = b5 = 0` | `x0`, `y0` (with `"principal"`) |
#' | `"additive"` | `b3 = b4 = b5 = 0` | none |
#' | `"mean"`, `"diff"` | `b3 = b4 = b5 = 0`, and `b1 = b2` or `b1 = -b2` | none |
#' | `"onlyx2"`, `"onlyy2"` | `b2 = b4 = b5 = 0` or `b1 = b3 = b4 = 0` | none |
#' | `"onlyx"`, `"onlyy"` | all terms but `b1` or `b2` are 0 | none |
#' | `"null"` | all five terms are 0 | none |
#' | `"addx2"`, `"addy2"` | `b4 = b5 = 0` or `b3 = b4 = 0` | none |
#' | `"noIA"` | `b4 = 0` | none |
#' | `"a5zero"` | `b3 = b5` | none |
#'
#' The last four models are additions to the family of `RSA::RSA()`.
#' `"addx2"` and `"addy2"` add one squared term to the plane; unlike
#' `"onlyx2"` and `"onlyy2"`, they keep both linear terms. `"noIA"` keeps
#' both squared terms and drops the product. `"a5zero"` sets
#' the two squared terms equal, so that `a5 = 0`; it is the smallest model
#' that nests both `"SRR"` and `"IA"`, which makes it a reference model for
#' testing either of them (see `reference` in [compare_rsa_mplus_models()]).
#'
#' The ridge parameters correspond to those of `RSA::RSA()`: `c` is the
#' lateral shift C of the ridge in units of X, `s` the scaling S of X that
#' rotates the ridge, `bm` the effect of the mean level along the ridge
#' (`meaneffect`), and `a4r` the curvature across the rotated ridge
#' (`a4.rescaled`, a4'). Each is written in its simplest closed form in
#' `b1` to `b5`; for example, a4' is `4*b5`. The rotated models `"SRRR"` and
#' `"SRSQD"` express `b3` (and `b1`) as functions of the other weights. One
#' run therefore covers both ridges and troughs, where `RSA::RSA()` fits an
#' inequality-constrained model twice. These models are not defined at
#' `b5 = 0`: a flat-ridge model that does not fit the data can drift towards
#' `b5 = 0` with a large `c`.
#'
#' Surface parameters (`cs`, `cc`, `is`, `ic`, `a5`) that a model fixes at
#' zero are not generated, and the others omit the terms fixed with `@0`,
#' because Mplus cannot refer to such a term in `MODEL CONSTRAINT`. With
#' `"principal"`, the full model gets its stationary point and principal axes
#' and `"IA"` its stationary point; the principal axes of `"IA"` have slopes
#' of 1 and -1 by construction, and the other models have no stationary
#' point. Restricted models keep the `XWITH` products, with fixed weights
#' where needed, so all models are nested in the full model and their
#' log-likelihoods can be compared. The `"weak"` and `"strong"` fit patterns
#' of `RSA::RSA()` are not generated.
#'
#' With `identification = "invariance"`, the latent line X = Y is equality on
#' the response scale because matched items measure X and Y with the same
#' loadings and intercepts. This is an assumption about the items: it holds
#' only if the X and Y items are parallel in wording and response format, and
#' it can be tested with [fit_rsa_mplus_invariance()]. Each factor mean is
#' then the mean of its marker item. With `center = "none"`, the origin
#' X = Y = 0 is the value 0 of the response scale, which usually lies outside
#' the data; `center = "pooled"` moves it to the mean of the two scale means.
#' To centre elsewhere, for example at the scale midpoint, use
#' `center = "none"` and subtract that value from the X and Y indicators in
#' `blocks$DEFINE`. The outcome factor keeps its marker identification.
#'
#' @return An object of class `rsa_mplus_workflow`. The object is updated and
#'   returned by the write, run, and read functions. Under pooled centering,
#'   `spec$centering` holds the subtracted constant and the two scale means it
#'   was computed from. Under invariance identification, `spec$invariance`
#'   holds the matched indicators and the freed pairs. `spec$rsa_model` names
#'   the generated model and
#'   `spec$start` holds the start values that were used.
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
  identification = c("marker", "effects", "invariance"),
  center = c("none", "pooled"),
  free_intercepts = NULL,
  rsa_model = "full",
  start = NULL,
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
  identification <- match.arg(identification)
  center <- match.arg(center)
  validate_rsa_mplus_identification(
    model_type,
    identification,
    center,
    free_intercepts
  )
  constraints <- match.arg(
    constraints,
    choices = c("surface", "principal"),
    several.ok = TRUE
  )
  constraints <- unique(c("surface", constraints))
  if (identical(rsa_model, "weak")) {
    rlang::abort(c(
      "The `weak` fit pattern is not generated.",
      i = "Its restriction `b3*b5 > 0` is an inequality."
    ))
  }
  if (identical(rsa_model, "strong")) {
    rlang::abort(c(
      "The `strong` fit pattern is not generated.",
      i = "It has the restrictions of `rsa_model = \"SRSQD\"`."
    ))
  }
  rsa_model <- rlang::arg_match(rsa_model, names(rsa_mplus_model_restrictions))
  model_spec <- rsa_mplus_model_spec(rsa_model)
  start_values <- rsa_mplus_start_values(start, model_spec)

  if (is.null(usevariables)) {
    usevariables <- names(data)
  }
  validate_rsa_mplus_variables(data, usevariables)

  tidy_model <- tidySEM::tidy_sem(data[, usevariables, drop = FALSE])
  tidy_model <- add_rsa_mplus_measurement(tidy_model, measurement)
  tidy_model <- add_rsa_mplus_structural(tidy_model, structural, roles)

  syntax_table <- tidySEM::syntax(tidy_model)
  syntax_table <- label_rsa_mplus_core(syntax_table, roles)
  syntax_table <- fix_rsa_mplus_linear_terms(syntax_table, model_spec$zero)

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

  effects <- NULL
  invariance <- NULL
  centering <- NULL
  if (identical(identification, "effects")) {
    effects <- apply_rsa_mplus_effects_coding(syntax_table, roles, center)
    syntax_table <- effects$syntax_table
  }
  if (identical(identification, "invariance")) {
    invariance <- apply_rsa_mplus_invariance(
      syntax_table,
      roles,
      free_intercepts
    )
    syntax_table <- invariance$syntax_table
  }
  if (identical(center, "pooled")) {
    centering <- rsa_mplus_pooled_centering(
      data,
      c(effects$indicators, invariance$indicators)
    )
  }
  tidySEM::syntax(tidy_model) <- syntax_table

  new_labels <- rsa_mplus_new_parameters(constraints, model_spec)$labels
  clashes <- intersect(new_labels, tolower(syntax_table$label))
  if (length(clashes) > 0L) {
    rlang::abort(sprintf(
      "Parameter label `%s` is reserved for the RSA model.",
      clashes[[1]]
    ))
  }

  conventional_model <- tidySEM::as_mplus(tidy_model)
  if (!is.null(effects)) {
    conventional_model <- free_rsa_mplus_first_loadings(
      conventional_model,
      effects$indicators
    )
  }
  for (label in intersect(names(start_values), c("b1", "b2"))) {
    conventional_model <- set_rsa_mplus_start(
      conventional_model,
      outcome = roles[["outcome"]],
      predictor = roles[[if (label == "b1") "x" else "y"]],
      label = label,
      value = start_values[[label]]
    )
  }
  modifiers <- rsa_mplus_term_modifiers(model_spec, start_values)
  custom_model <- c(
    sprintf("XY | %s XWITH %s;", roles[["x"]], roles[["y"]]),
    sprintf("XS | %s XWITH %s;", roles[["x"]], roles[["x"]]),
    sprintf("YS | %s XWITH %s;", roles[["y"]], roles[["y"]]),
    sprintf("%s ON XS%s (b3);", roles[["outcome"]], modifiers[["b3"]]),
    sprintf("%s ON XY%s (b4);", roles[["outcome"]], modifiers[["b4"]]),
    sprintf("%s ON YS%s (b5);", roles[["outcome"]], modifiers[["b5"]])
  )
  model_syntax <- paste(
    c(conventional_model, custom_model, flatten_rsa_mplus_text(model_extra)),
    collapse = "\n"
  )

  constraint_syntax <- rsa_mplus_constraints(constraints, model_spec)
  constraint_syntax <- paste(
    c(
      constraint_syntax,
      effects$constraints,
      flatten_rsa_mplus_text(constraint_extra)
    ),
    collapse = "\n"
  )
  # The null model has nothing to constrain, so the block is left out.
  if (!nzchar(constraint_syntax)) {
    constraint_syntax <- NULL
  }

  block_args <- validate_rsa_mplus_blocks(blocks)
  defaults <- list(
    TITLE = "RSA model generated by rsahelpers;",
    ANALYSIS = "TYPE = RANDOM;\nALGORITHM = INTEGRATION;",
    OUTPUT = "STDYX CINTERVAL TECH1;"
  )
  defaults[names(block_args)] <- block_args
  block_args <- defaults

  # The pooled shift comes first so that any user DEFINE statements see the
  # centred indicators, exactly as the model does.
  if (!is.null(centering)) {
    block_args$DEFINE <- paste(
      c(centering$define, block_args$DEFINE),
      collapse = "\n"
    )
  }

  block_args$VARIABLE <- rsa_mplus_variable_block(
    usevariables = usevariables,
    extra = block_args$VARIABLE
  )

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
    coefficient_labels = rsa_mplus_coefficient_labels,
    rsa_model = rsa_model
  )

  new_rsa_mplus_workflow(
    spec = list(
      roles = roles,
      model_type = model_type,
      reliability = reliability,
      constraints = constraints,
      usevariables = usevariables,
      identification = identification,
      center = center,
      centering = centering[c("constant", "scale_means")],
      invariance = invariance[c("indicators", "free_intercepts")],
      rsa_model = rsa_model,
      start = start_values,
      metadata = metadata
    ),
    tidysem = tidy_model,
    mplus = mplus_object,
    status = "created"
  )
}

#' @param x An `rsa_mplus_workflow` object.
#' @param modelout Path to the Mplus `.inp` file to create. Mplus lower-cases
#'   file names, so on a case-sensitive file system (Linux, WSL) the file name
#'   must be lower case for Mplus to find it.
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
    # Mplus lower-cases file names. On a case-sensitive file system it then
    # looks for a lower-case input file and reports the failure next to it.
    lower_path <- file.path(
      dirname(output_path),
      tolower(basename(output_path))
    )
    rlang::abort(c(
      sprintf(
        "Mplus did not create the expected output file `%s`.",
        output_path
      ),
      if (!identical(lower_path, output_path) && file.exists(lower_path)) {
        c(
          i = sprintf("Mplus wrote `%s` instead.", lower_path),
          i = "Mplus lower-cases file names; on a case-sensitive file system, give `modelout` a lower-case file name."
        )
      }
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
#' @param output Optional Mplus `.out` path to read into an existing workflow.
#'   This is useful for attaching a cached or relocated output while retaining
#'   the workflow's role metadata. It can only be used when `target` is an
#'   `rsa_mplus_workflow`.
#' @rdname create_rsa_mplus_model
#' @export
read_rsa_mplus_model <- function(target, output = NULL, quiet = TRUE) {
  require_rsa_mplus_namespace("MplusAutomation")

  if (inherits(target, "rsa_mplus_workflow")) {
    x <- target
    output_path <- output
    if (is.null(output_path)) {
      output_path <- x$files$output
    }
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

  if (!is.null(output)) {
    rlang::abort(
      "`output` can only be supplied when `target` is an RSA-Mplus workflow."
    )
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
  identification = c("marker", "effects", "invariance"),
  center = c("none", "pooled"),
  free_intercepts = NULL,
  rsa_model = "full",
  start = NULL,
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
    identification = identification,
    center = center,
    free_intercepts = free_intercepts,
    rsa_model = rsa_model,
    start = start,
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
  if (!is.null(x$spec$rsa_model)) {
    cat("RSA model:", x$spec$rsa_model, "\n")
  }
  if (!is.null(x$spec$centering)) {
    cat(
      "Centering: pooled at",
      format(x$spec$centering$constant, digits = 6),
      "\n"
    )
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

validate_rsa_mplus_identification <- function(
  model_type,
  identification,
  center,
  free_intercepts = NULL
) {
  if (identical(identification, "effects") && !identical(model_type, "latent")) {
    rlang::abort(
      "Effects coding requires `model_type = \"latent\"` (multi-indicator factors)."
    )
  }
  if (
    identical(identification, "invariance") && !identical(model_type, "latent")
  ) {
    rlang::abort(
      "Invariance identification requires `model_type = \"latent\"` (multi-indicator factors)."
    )
  }
  # Under plain marker identification the X and Y factors are not on a common
  # metric, so a shared origin would not make X = Y equality.
  if (identical(center, "pooled") && identical(identification, "marker")) {
    rlang::abort(
      "Pooled centering requires `identification = \"effects\"` or `\"invariance\"`."
    )
  }
  if (!is.null(free_intercepts)) {
    if (!identical(identification, "invariance")) {
      rlang::abort(
        "`free_intercepts` requires `identification = \"invariance\"`."
      )
    }
    if (
      !is.character(free_intercepts) ||
        anyNA(free_intercepts) ||
        any(!nzchar(free_intercepts))
    ) {
      rlang::abort("`free_intercepts` must be a character vector of indicator names.")
    }
  }
  invisible(NULL)
}

# Scalar invariance of the X and Y factors (shared loadings and intercepts of
# matched indicators) with marker intercepts at 0 and free factor means. The
# two factors are then on one metric with one origin, without MODEL CONSTRAINT
# statements, which Mplus does not accept for these parameters under Bayes.
apply_rsa_mplus_invariance <- function(syntax_table, roles, free_intercepts) {
  labelled <- !is.na(syntax_table$label) & nzchar(syntax_table$label)
  rows <- list()
  indicators <- list()
  for (role in c("x", "y")) {
    factor_name <- roles[[role]]
    loading_rows <- which(
      syntax_table$op == "=~" & syntax_table$lhs == factor_name
    )
    default_loadings <- syntax_table$free[loading_rows[1]] == 0 &&
      isTRUE(syntax_table$ustart[loading_rows[1]] == 1) &&
      all(syntax_table$free[loading_rows[-1]] == 1) &&
      !any(labelled[loading_rows])
    if (!default_loadings) {
      rlang::abort(sprintf(
        "Invariance identification needs default (unlabelled, marker-identified) loadings for `%s`.",
        factor_name
      ))
    }
    factor_indicators <- syntax_table$rhs[loading_rows]
    intercept_rows <- which(syntax_table$op == "~1")[match(
      factor_indicators,
      syntax_table$lhs[syntax_table$op == "~1"]
    )]
    mean_row <- which(syntax_table$op == "~1" & syntax_table$lhs == factor_name)
    if (anyNA(intercept_rows) || length(mean_row) != 1L) {
      rlang::abort(sprintf(
        "Could not identify the mean structure of `%s`.",
        factor_name
      ))
    }
    if (any(labelled[c(intercept_rows, mean_row)])) {
      rlang::abort(sprintf(
        "Invariance identification needs an unlabelled mean structure for `%s`.",
        factor_name
      ))
    }
    rows[[role]] <- list(
      loadings = loading_rows,
      intercepts = intercept_rows,
      mean = mean_row
    )
    indicators[[role]] <- stats::setNames(factor_indicators, NULL)
  }

  n_pairs <- length(indicators$x)
  if (length(indicators$y) != n_pairs) {
    rlang::abort(sprintf(
      "Invariance identification needs as many `%s` as `%s` indicators, in matching order (%d and %d given).",
      roles[["x"]],
      roles[["y"]],
      n_pairs,
      length(indicators$y)
    ))
  }

  free_pairs <- integer()
  if (length(free_intercepts) > 0L) {
    unknown <- setdiff(free_intercepts, unlist(indicators))
    if (length(unknown) > 0L) {
      rlang::abort(sprintf(
        "`free_intercepts` must name X or Y indicators; not found: %s.",
        paste(sprintf("`%s`", unknown), collapse = ", ")
      ))
    }
    free_pairs <- sort(unique(c(
      match(free_intercepts, indicators$x),
      match(free_intercepts, indicators$y)
    )))
    if (1L %in% free_pairs) {
      rlang::abort(
        "The intercepts of the marker indicators identify the factor means and cannot be freed."
      )
    }
  }

  pairs <- seq_len(n_pairs)[-1]
  loading_labels <- paste0("lxy", pairs)
  shared <- setdiff(pairs, free_pairs)
  intercept_labels <- paste0("ixy", shared)
  clashes <- intersect(
    c(loading_labels, intercept_labels),
    tolower(syntax_table$label[labelled])
  )
  if (length(clashes) > 0L) {
    rlang::abort(sprintf(
      "Parameter label `%s` is reserved for invariance identification.",
      clashes[[1]]
    ))
  }

  for (role in c("x", "y")) {
    r <- rows[[role]]
    syntax_table$label[r$loadings[pairs]] <- loading_labels
    syntax_table$label[r$intercepts[shared]] <- intercept_labels
    syntax_table$free[r$intercepts[[1]]] <- 0
    syntax_table$ustart[r$intercepts[[1]]] <- 0
    syntax_table$free[r$mean] <- 1
    syntax_table$ustart[r$mean] <- NA
  }

  list(
    syntax_table = syntax_table,
    indicators = indicators,
    free_intercepts = unlist(
      lapply(indicators, function(items) items[free_pairs]),
      use.names = FALSE
    )
  )
}

# Effects coding (Little, Slegers, & Card, 2006) for the three RSA factors:
# loadings average 1, so each factor is on the metric of its average indicator.
# Under pooled centering the X and Y intercepts also sum to zero and the two
# factor means are freed, so those means become the centred scale means.
apply_rsa_mplus_effects_coding <- function(syntax_table, roles, center) {
  role_keys <- c(x = "x", y = "y", outcome = "z")
  labelled <- !is.na(syntax_table$label) & nzchar(syntax_table$label)
  existing_labels <- tolower(syntax_table$label[labelled])
  indicators <- list()
  constraints <- character()

  for (role in names(role_keys)) {
    factor_name <- roles[[role]]
    key <- role_keys[[role]]
    loading_rows <- which(
      syntax_table$op == "=~" & syntax_table$lhs == factor_name
    )
    default_loadings <- syntax_table$free[loading_rows[1]] == 0 &&
      isTRUE(syntax_table$ustart[loading_rows[1]] == 1) &&
      all(syntax_table$free[loading_rows[-1]] == 1) &&
      !any(labelled[loading_rows])
    if (!default_loadings) {
      rlang::abort(sprintf(
        "Effects coding needs default (unlabelled, marker-identified) loadings for `%s`.",
        factor_name
      ))
    }

    factor_indicators <- syntax_table$rhs[loading_rows]
    loading_labels <- paste0("l", key, seq_along(loading_rows))
    intercept_labels <- paste0("i", key, seq_along(loading_rows))
    new_labels <- loading_labels
    syntax_table$free[loading_rows] <- 1
    syntax_table$ustart[loading_rows] <- NA
    syntax_table$label[loading_rows] <- loading_labels
    constraints <- c(
      constraints,
      rsa_mplus_sum_constraint(loading_labels, total = length(loading_rows))
    )

    if (identical(center, "pooled") && role != "outcome") {
      intercept_rows <- match(
        factor_indicators,
        syntax_table$lhs[syntax_table$op == "~1"]
      )
      intercept_rows <- which(syntax_table$op == "~1")[intercept_rows]
      mean_row <- which(
        syntax_table$op == "~1" & syntax_table$lhs == factor_name
      )
      if (anyNA(intercept_rows) || length(mean_row) != 1L) {
        rlang::abort(sprintf(
          "Could not identify the mean structure of `%s`.",
          factor_name
        ))
      }
      if (any(labelled[c(intercept_rows, mean_row)])) {
        rlang::abort(sprintf(
          "Pooled centering needs an unlabelled mean structure for `%s`.",
          factor_name
        ))
      }
      syntax_table$label[intercept_rows] <- intercept_labels
      syntax_table$free[mean_row] <- 1
      syntax_table$ustart[mean_row] <- NA
      new_labels <- c(new_labels, intercept_labels)
      constraints <- c(
        constraints,
        rsa_mplus_sum_constraint(intercept_labels, total = 0)
      )
    }

    clashes <- intersect(tolower(new_labels), existing_labels)
    if (length(clashes) > 0L) {
      rlang::abort(sprintf(
        "Parameter label `%s` is reserved for effects coding.",
        clashes[[1]]
      ))
    }
    indicators[[role]] <- stats::setNames(factor_indicators, NULL)
  }

  list(
    syntax_table = syntax_table,
    constraints = constraints,
    indicators = indicators
  )
}

# The first label is solved from the others, e.g. `lx1 = 3 - lx2 - lx3;`.
rsa_mplus_sum_constraint <- function(labels, total) {
  rest <- paste(labels[-1], collapse = " - ")
  rhs <- if (total == 0) paste0("-", rest) else paste(total, "-", rest)
  statement <- sprintf("%s = %s;", labels[[1]], rhs)
  # Mplus reads at most 90 characters per line.
  paste(strwrap(statement, width = 76L, exdent = 2L), collapse = "\n")
}

# tidySEM drops the asterisk for a freed first loading, and Mplus would then
# silently restore the marker constraint, so the asterisk is added here.
free_rsa_mplus_first_loadings <- function(model_syntax, indicators) {
  for (role in names(indicators)) {
    first <- indicators[[role]][[1]]
    pattern <- sprintf("(\\bBY\\s+%s)(\\s*\\()", first)
    hits <- gregexpr(pattern, model_syntax, ignore.case = TRUE, perl = TRUE)
    if (sum(unlist(hits) > 0L) != 1L) {
      rlang::abort(sprintf(
        "Could not free the first loading of indicator `%s`.",
        first
      ))
    }
    model_syntax <- sub(
      pattern,
      "\\1*\\2",
      model_syntax,
      ignore.case = TRUE,
      perl = TRUE
    )
  }
  model_syntax
}

# Pooled centering subtracts one constant, the mean of the X and Y scale means,
# from every X and Y indicator. With effects coding the latent means are then
# the two centred scale means, so the line X = Y is literal equality on the
# shared response scale - the latent counterpart of `RSA(center = "pooled")`.
rsa_mplus_pooled_centering <- function(data, indicators) {
  scale_means <- vapply(
    indicators[c("x", "y")],
    function(items) {
      values <- data[items]
      if (!all(vapply(values, is.numeric, logical(1)))) {
        rlang::abort("Pooled centering requires numeric X and Y indicators.")
      }
      mean(colMeans(values, na.rm = TRUE))
    },
    numeric(1)
  )
  if (anyNA(scale_means)) {
    rlang::abort("Every X and Y indicator needs at least one observed value.")
  }
  constant <- mean(scale_means)
  operator <- if (constant < 0) "+" else "-"
  value <- formatC(abs(constant), digits = 10, format = "fg")
  items <- unlist(indicators[c("x", "y")], use.names = FALSE)

  list(
    constant = constant,
    scale_means = scale_means,
    define = sprintf("%s = %s %s %s;", items, items, operator, trimws(value))
  )
}

rsa_mplus_coefficient_labels <- c(
  x = "b1",
  y = "b2",
  x2 = "b3",
  xy = "b4",
  y2 = "b5"
)

# Surface parameters a1-a5 as weights on b1-b5.
rsa_mplus_surface_weights <- list(
  cs = c(b1 = 1, b2 = 1),
  cc = c(b3 = 1, b4 = 1, b5 = 1),
  is = c(b1 = 1, b2 = -1),
  ic = c(b3 = 1, b4 = -1, b5 = 1),
  a5 = c(b3 = 1, b5 = -1)
)

# Order of the principal-axis parameters in `NEW`. They are defined in a
# different order, because p10 and p20 use the slopes.
rsa_mplus_principal_order <- c("x0", "y0", "p10", "p11", "p20", "p21")

# Restrictions that each second-order model of `RSA::RSA()` places on b1-b5
# beyond the zero terms in `rsa_mplus_zero_terms`. `linear` gives a dependent
# weight as weights on the free ones; `nonlinear` gives it as an Mplus
# expression in the free ones. The rotated ridge models use `b3 = b4^2/(4*b5)`
# because the inequality-constrained form of `RSA::RSA()` (b4^2 = 4*b3*b5 with
# b3, b5 < 0) collapses to b3 = b4 = b5 = 0 in Mplus. `derived` holds the ridge
# parameters C, S, bM, and a4', `principal` the stationary point and axes.
rsa_mplus_model_restrictions <- list(
  full = list(
    principal = c(
      x0 = "(b2*b4 - 2*b1*b5) / (4*b3*b5 - b4*b4)",
      y0 = "(b1*b4 - 2*b2*b3) / (4*b3*b5 - b4*b4)",
      p11 = "(b5 - b3 + SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4",
      p21 = "(b5 - b3 - SQRT((b3 - b5)*(b3 - b5) + b4*b4)) / b4",
      p10 = "y0 - p11*x0",
      p20 = "y0 - p21*x0"
    )
  ),
  SRRR = list(
    nonlinear = c(b3 = "b4*b4/(4*b5)"),
    derived = c(
      c = "-(2*b1*b5 + b2*b4)/(4*b4*b5)",
      s = "-b4/(2*b5)",
      bm = "(b2*b4 - 2*b1*b5)/b4",
      a4r = "4*b5"
    )
  ),
  SRR = list(
    linear = list(b3 = c(b5 = 1), b4 = c(b5 = -2)),
    derived = c(c = "(b1 - b2)/(4*b3)", bm = "b1 + b2")
  ),
  RR = list(
    linear = list(b1 = c(b2 = 1), b3 = c(b5 = 1), b4 = c(b5 = -2)),
    derived = c(bm = "b1 + b2")
  ),
  SRSQD = list(
    nonlinear = c(b1 = "b2*b4/(2*b5)", b3 = "b4*b4/(4*b5)"),
    derived = c(c = "-b2/(2*b5)", s = "-b4/(2*b5)", a4r = "4*b5")
  ),
  SSQD = list(
    linear = list(b1 = c(b2 = -1), b3 = c(b5 = 1), b4 = c(b5 = -2)),
    derived = c(c = "b1/(2*b3)")
  ),
  SQD = list(linear = list(b3 = c(b5 = 1), b4 = c(b5 = -2))),
  IA = list(principal = c(x0 = "-b2/b4", y0 = "-b1/b4")),
  additive = list(),
  mean = list(linear = list(b1 = c(b2 = 1))),
  diff = list(linear = list(b1 = c(b2 = -1))),
  onlyx2 = list(),
  onlyy2 = list(),
  onlyx = list(),
  onlyy = list(),
  null = list(),
  # Not in the RSA package: the plane with one squared term, the full model
  # without the product term, and the full model with equal squared terms
  # (a5 = 0), the smallest model that nests both SRR and IA.
  addx2 = list(),
  addy2 = list(),
  noIA = list(),
  a5zero = list(linear = list(b3 = c(b5 = 1)))
)

rsa_mplus_model_spec <- function(rsa_model) {
  restrictions <- rsa_mplus_model_restrictions[[rsa_model]]
  zero <- unname(rsa_mplus_coefficient_labels[
    rsa_mplus_zero_terms[[rsa_model]]
  ])
  dependent <- c(names(restrictions$linear), names(restrictions$nonlinear))
  list(
    rsa_model = rsa_model,
    zero = zero,
    linear = restrictions$linear,
    nonlinear = restrictions$nonlinear,
    free = unname(setdiff(rsa_mplus_coefficient_labels, c(zero, dependent))),
    principal = restrictions$principal,
    derived = restrictions$derived
  )
}

# Renders weights such as c(b3 = 1, b4 = -2) as "b3 - 2*b4".
rsa_mplus_linear_expression <- function(weights) {
  magnitudes <- vapply(abs(weights), format, character(1))
  terms <- ifelse(
    magnitudes == "1",
    names(weights),
    paste0(magnitudes, "*", names(weights))
  )
  signs <- ifelse(weights < 0, " - ", " + ")
  signs[[1]] <- if (weights[[1]] < 0) "-" else ""
  paste0(signs, terms, collapse = "")
}

# Surface parameters that are not identically zero under the model. Their
# definitions leave out the terms fixed with `@0`, because Mplus cannot refer
# to the label of a fixed parameter in MODEL CONSTRAINT.
rsa_mplus_surface_definitions <- function(spec) {
  basis <- setdiff(
    rsa_mplus_coefficient_labels,
    c(spec$zero, names(spec$linear))
  )
  expand <- function(term) {
    expanded <- stats::setNames(numeric(length(basis)), basis)
    if (term %in% spec$zero) {
      return(expanded)
    }
    weights <- spec$linear[[term]]
    if (is.null(weights)) {
      weights <- stats::setNames(1, term)
    }
    expanded[names(weights)] <- weights
    expanded
  }

  definitions <- vapply(
    rsa_mplus_surface_weights,
    function(weights) {
      total <- Reduce(
        `+`,
        Map(function(term, weight) weight * expand(term), names(weights), weights)
      )
      if (all(abs(total) < 1e-12)) {
        return(NA_character_)
      }
      rsa_mplus_linear_expression(weights[!names(weights) %in% spec$zero])
    },
    character(1)
  )
  definitions[!is.na(definitions)]
}

rsa_mplus_new_parameters <- function(constraints, spec) {
  surface <- rsa_mplus_surface_definitions(spec)
  principal <- if ("principal" %in% constraints) spec$principal
  list(
    labels = c(
      names(surface),
      intersect(rsa_mplus_principal_order, names(principal)),
      names(spec$derived)
    ),
    definitions = c(surface, principal, spec$derived)
  )
}

rsa_mplus_constraints <- function(constraints, spec) {
  definitions <- function(x) {
    if (length(x) == 0L) character() else sprintf("%s = %s;", names(x), x)
  }
  restrictions <- c(
    vapply(spec$linear, rsa_mplus_linear_expression, character(1)),
    spec$nonlinear
  )
  new_parameters <- rsa_mplus_new_parameters(constraints, spec)
  c(
    definitions(restrictions),
    if (length(new_parameters$labels) > 0L) {
      sprintf("NEW(%s);", paste(new_parameters$labels, collapse = " "))
    },
    definitions(new_parameters$definitions)
  )
}

# b1 and b2 come from tidySEM, which writes `@0` for a fixed path.
fix_rsa_mplus_linear_terms <- function(syntax_table, zero) {
  rows <- which(
    syntax_table$op == "~" &
      tolower(syntax_table$label) %in% intersect(zero, c("b1", "b2"))
  )
  syntax_table$free[rows] <- 0
  syntax_table$ustart[rows] <- 0
  syntax_table
}

# Suffix of each quadratic regression: `@0` for a term fixed at zero and
# `*value` for a start value.
rsa_mplus_term_modifiers <- function(spec, start_values) {
  modifiers <- stats::setNames(
    rep("", length(rsa_mplus_coefficient_labels)),
    rsa_mplus_coefficient_labels
  )
  modifiers[names(start_values)] <- paste0(
    "*",
    vapply(start_values, rsa_mplus_number, character(1))
  )
  modifiers[spec$zero] <- "@0"
  modifiers
}

rsa_mplus_number <- function(x) {
  trimws(formatC(x, digits = 6, format = "fg"))
}

rsa_mplus_start_values <- function(start, spec, call = rlang::caller_env()) {
  if (is.null(start)) {
    return(stats::setNames(numeric(), character()))
  }
  if (inherits(start, "rsa_mplus_workflow")) {
    start <- RSA_mplus(start, plot = FALSE)
  }
  if (inherits(start, "rsa_mplus")) {
    start <- start$coefficients
  }
  if (!is.numeric(start) || is.null(names(start))) {
    rlang::abort(
      "`start` must be a named numeric vector, an `RSA_mplus()` result, or a fitted workflow.",
      call = call
    )
  }

  start <- start[names(start) != "b0"]
  aliases <- c(
    rsa_mplus_coefficient_labels,
    stats::setNames(rsa_mplus_coefficient_labels, rsa_mplus_coefficient_labels)
  )
  unknown <- setdiff(names(start), names(aliases))
  if (length(unknown) > 0L) {
    rlang::abort(
      sprintf(
        "Unknown name%s in `start`: %s. Use `b1` to `b5`.",
        if (length(unknown) == 1L) "" else "s",
        paste(sprintf("`%s`", unknown), collapse = ", ")
      ),
      call = call
    )
  }
  names(start) <- unname(aliases[names(start)])
  if (anyDuplicated(names(start)) || any(!is.finite(start))) {
    rlang::abort(
      "`start` must name each term once with a finite value.",
      call = call
    )
  }
  start[intersect(spec$free, names(start))]
}

# Adds a start value to a path written by tidySEM, e.g. `Z ON X (b1);`.
set_rsa_mplus_start <- function(model_syntax, outcome, predictor, label, value) {
  pattern <- sprintf(
    "(\\b%s\\s+ON\\s+%s)(\\s*\\(%s\\)\\s*;)",
    outcome,
    predictor,
    label
  )
  hits <- gregexpr(pattern, model_syntax, ignore.case = TRUE, perl = TRUE)
  if (sum(unlist(hits) > 0L) != 1L) {
    rlang::abort(sprintf("Could not set the start value of `%s`.", label))
  }
  sub(
    pattern,
    paste0("\\1*", rsa_mplus_number(value), "\\2"),
    model_syntax,
    ignore.case = TRUE,
    perl = TRUE
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

rsa_mplus_variable_block <- function(usevariables, extra = NULL) {
  if (
    !is.null(extra) &&
      grepl("\\bUSEVARIABLES\\b", extra, ignore.case = TRUE, perl = TRUE)
  ) {
    rlang::abort(
      paste(
        "Do not set `USEVARIABLES` through `blocks$VARIABLE`;",
        "use the `usevariables` argument instead."
      )
    )
  }

  variable_lines <- c(
    "USEVARIABLES =",
    paste0("  ", strwrap(paste(usevariables, collapse = " "), width = 76L)),
    ";"
  )
  paste(c(variable_lines, flatten_rsa_mplus_text(extra)), collapse = "\n")
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

  # The `.exe` names find a Windows installation of Mplus from WSL, where the
  # Windows `PATH` is searched but the extension is not implied.
  candidates <- Sys.which(c("mplus", "Mplus", "mplus.exe", "Mplus.exe"))
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
