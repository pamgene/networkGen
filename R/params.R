#' Capture the parameters used for a run
#'
#' Generic capture only: evaluates `...` and records each parameter's value
#' plus a human-readable label for it, so a caller can recover *what* was
#' used without re-deriving it from context. Has no knowledge of what any
#' particular parameter name means -- see [build_param_folder()] for the
#' hardcoded rule that turns a captured `uka`/`ppi_network`/`spec_cutoff`/...
#' shaped `params` into a folder name, and [prepare_run_params()] for the
#' two composed together plus the filesystem side effects (folder creation,
#' writing `params.csv`).
#'
#' For a scalar (length-1 atomic) parameter, the label is the evaluated
#' *value* -- so `capture_params(spec_cutoff = sc)` inside a loop over `sc`
#' records `0.7`, not the loop variable's name `sc`. For anything else (a
#' data frame, a vector), the label is the call-site expression text (e.g.
#' the variable name `uka`) instead, since dumping a data frame's full
#' content wouldn't be readable.
#'
#' @param ... Named (or unnamed, using the argument's own variable name)
#'   parameters to capture.
#' @param .labels Named list/character vector overriding specific
#'   parameters' labels, instead of deriving them from the expression
#'   written at this call site. For when the immediate caller is itself a
#'   relay -- e.g. `run_network_grid()` re-binding its own `raw_uka`
#'   argument to `uka` before calling this -- the expression visible *here*
#'   would just be the relay's own local variable name (`raw_uka`), not
#'   the original caller's (`cleaned_batch`). A relay should capture its own
#'   caller's expression with `rlang::as_label(rlang::enquo(raw_uka))`
#'   *before* re-binding, and pass it through here.
#'
#' @return The evaluated parameters as a named list, with a `"symbol_names"`
#'   attribute recording each parameter's label as text.
#' @export
capture_params <- function(..., .labels = list()) {
  dots <- rlang::enquos(...)
  param_names <- names(dots)
  if (is.null(param_names)) {
    param_names <- rep("", length(dots))
  }

  expr_labels <- purrr::map_chr(dots, rlang::as_label)
  param_names <- ifelse(param_names == "" | is.na(param_names), expr_labels, param_names)

  values <- purrr::map(dots, rlang::eval_tidy)

  is_scalar <- purrr::map_lgl(values, ~ is.atomic(.x) && length(.x) == 1 && !is.null(.x))
  symbol_names <- expr_labels
  symbol_names[is_scalar] <- purrr::map_chr(values[is_scalar], as.character)

  override_idx <- match(names(.labels), param_names)
  has_override <- !is.na(override_idx)
  symbol_names[override_idx[has_override]] <- unlist(.labels)[has_override]

  params <- rlang::set_names(values, param_names)
  attr(params, "symbol_names") <- rlang::set_names(symbol_names, param_names)
  params
}

#' Build a parameter-encoded folder name for a networkGen run
#'
#' The hardcoded naming rule: reads a fixed set of recognized parameter
#' names out of `params` (as returned by [capture_params()]) and encodes
#' them into one folder-name string, one folder per unique parameter
#' combination -- e.g. 4 conditions x 2 `spec_cutoff` values is 8 folders,
#' not one folder covering a sweep of cutoffs. Pure function: no filesystem
#' side effects, no dependence on `respath` (joining onto a base path and
#' creating the directory is [prepare_run_params()]'s job).
#'
#' Only understands `networkGen`-native parameters. `networkScore` composes
#' its own folder-naming on top of this rather than this function knowing
#' about scoring-only parameters like `nperms_network` -- see
#' `networkScore::build_score_param_folder()`.
#'
#' @param params A parameters list from [capture_params()] (must carry the
#'   `"symbol_names"` attribute it sets). Recognized names: `uka`, `sens`
#'   (optional -- nests the result under a `uka`-named subfolder when
#'   present, for paired analyses), `ppi_network`, `spec_cutoff` (single
#'   value or vector; omitted from the name entirely when `0`, since `0`
#'   means "no filtering on this cutoff"), `perc_cutoff` (same omit-if-0
#'   rule), `sens_perc_cutoff` (optional -- the sensitivity-side percentile
#'   cutoff, gridded independently of `perc_cutoff`; same omit-if-0 rule),
#'   `rank_uka_abs`, `b`, `w`, `cs` (optional -- recorded in `params.csv`
#'   via [save_params()], but never shown in the folder name), `art_nodes`
#'   (optional).
#'
#' @return The folder name (a single path segment, or `uka_name/components`
#'   when `sens` is present), as a string.
#' @export
build_param_folder <- function(params) {
  symbol_names <- attr(params, "symbol_names")

  uka_name <- symbol_names[["uka"]]
  sens_name <- if ("sens" %in% names(symbol_names)) symbol_names[["sens"]] else NULL
  ppi_network_name <- symbol_names[["ppi_network"]]

  spec_cutoff_shown <- !is.null(params$spec_cutoff) && !(length(params$spec_cutoff) == 1 && params$spec_cutoff == 0)
  perc_cutoff_shown <- !is.null(params$perc_cutoff) && !identical(params$perc_cutoff, 0)
  sens_perc_cutoff_shown <- !is.null(params$sens_perc_cutoff) && !identical(params$sens_perc_cutoff, 0)

  components <- c(
    if (!is.null(sens_name)) sens_name,
    if (spec_cutoff_shown) paste0("spec", paste0(params$spec_cutoff, collapse = "-")),
    if (perc_cutoff_shown) paste0("perc", params$perc_cutoff),
    if (sens_perc_cutoff_shown) paste0("sensperc", params$sens_perc_cutoff),
    ppi_network_name,
    paste0("ukaabs", as.integer(params$rank_uka_abs)),
    paste0("b", params$b),
    if (!is.null(params$w)) paste0("w", params$w),
    if (!is.null(params$art_nodes)) "art"
  )

  folder_name <- paste(components, collapse = "_")
  if (!is.null(sens_name)) file.path(uka_name, folder_name) else folder_name
}

#' Capture a run's parameters, build its output folder, and write `params.csv`
#'
#' Composes [capture_params()] and [build_param_folder()], then the
#' filesystem side effects: creates the resulting folder under `respath`
#' and writes `params.csv` there via [save_params()]. This is what real
#' callers use -- e.g. before [generate_kinase_network()]/
#' [generate_paired_network()] with `write = TRUE`, passing this function's
#' resulting `respath` on as their `res.path`.
#'
#' The following build call should read its arguments back off this
#' function's return value (e.g. `generate_kinase_network(spec_cutoff =
#' run_params$spec_cutoff, ...)`) rather than retyping the same literals --
#' otherwise the two calls can drift apart, and `params.csv` would then
#' record different values than what was actually built.
#'
#' @param ... Passed to [capture_params()]. Must include `respath` (the base
#'   output directory this run's parameter-encoded folder is created under)
#'   and whatever [build_param_folder()] needs (see its docs).
#' @param .labels Passed to [capture_params()] -- see its docs. Needed when
#'   *this* call is itself made by a relay (e.g. `run_network_grid()`)
#'   rather than directly by the end caller.
#'
#' @return The evaluated parameters as a named list (see [capture_params()]),
#'   with `respath` replaced by the created parameter-encoded folder path.
#' @export
prepare_run_params <- function(..., .labels = list()) {
  params <- capture_params(..., .labels = .labels)
  folder_name <- build_param_folder(params)
  param_folder <- file.path(params$respath, folder_name)
  print(paste0("param folder: ", param_folder))
  if (!dir.exists(param_folder)) dir.create(param_folder, recursive = TRUE)
  params$respath <- param_folder
  save_params(params, respath = param_folder)
  params
}

#' Write a run's parameters to `params.csv`
#'
#' @param params A parameters list from [capture_params()] (must carry the
#'   `"symbol_names"` attribute it sets).
#' @param respath Folder to write `params.csv` into.
#'
#' @return `params`' parameter data frame, invisibly (also written to disk).
#' @export
save_params <- function(params, respath) {
  symbol_names <- attr(params, "symbol_names")

  df <- tibble::tibble(
    parameter = names(symbol_names),
    value = unname(symbol_names)
  ) %>%
    dplyr::mutate(value = gsub('\\"', '"', .data$value)) %>%
    dplyr::mutate(value = gsub("c\\(|\\)", "", .data$value)) %>%
    dplyr::mutate(value = gsub('^"|"$', "", .data$value)) %>%
    dplyr::mutate(value = stringr::str_remove_all(.data$value, "\\\\")) %>%
    dplyr::mutate(value = stringr::str_remove_all(.data$value, '"')) %>%
    dplyr::filter(!.data$parameter %in% c("uka_fam", "respath"))

  readr::write_csv(df, paste0(respath, "/params.csv"))
  print(df)
  invisible(df)
}
