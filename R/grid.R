#' Normalize `ppi_network` to a named list of one or more reference networks
#'
#' `ppi_network` is either a single data frame (the common case -- wrapped
#' into a length-1 named list using `single_label`) or an explicitly, fully
#' named list of them (for gridding across more than one reference network).
#' Exported so `networkScore`'s scoring entry points (which grid `ppi_network`
#' the same way) can reuse it rather than duplicating this normalization.
#'
#' @param ppi_network A data frame, or a fully named list of data frames.
#' @param single_label Label to use as the list name when `ppi_network` is a
#'   single data frame -- the caller's own correctly-captured expression
#'   (see [run_network_grid()]/[build_network_grid()] for why this can't
#'   just be derived here with `rlang::enquo()`).
#'
#' @return A fully named list of data frames.
#' @export
normalize_ppi_network_list <- function(ppi_network, single_label) {
  if (is.data.frame(ppi_network)) {
    return(stats::setNames(list(ppi_network), single_label))
  }
  if (!is.list(ppi_network) || is.null(names(ppi_network)) || any(names(ppi_network) == "")) {
    stop(
      "ppi_network must be a single data frame, or a fully named list of them ",
      "(e.g. list(v12_5 = string_12.5, kins502 = ppi_networkv12_502_kins)) -- ",
      "each name is used for output-folder naming.",
      call. = FALSE
    )
  }
  ppi_network
}

#' Expand a parameter grid into a flat list of per-network build inputs
#'
#' The one place `dataset x condition x spec_cutoff x perc_cutoff x b x w x
#' rank_uka_abs x ppi_network` gets expanded -- both [run_network_grid()]
#' (plain generation) and `networkScore`'s scoring entry points build on
#' this, so there is exactly one grid-construction mechanism in the whole
#' suite (see `docs/adr/0006-unified-grid-interface.md`). Pure: no
#' filesystem access, no network builds -- just data in, a flat list of
#' resolved cells out. Every one of `spec_cutoff`, `perc_cutoff`, `b`, `w`,
#' `rank_uka_abs`, `ppi_network` can be a vector (or, for `ppi_network`, a
#' named list) -- every combination is gridded, not just the ones a
#' particular caller happens to demonstrate.
#'
#' @param raw_uka Raw kinase-activity data. If it has a `dataset_col` column
#'   (e.g. from [load_uka_dataset_files()]), it's split by that column
#'   before cleaning -- each dataset cleaned and gridded independently. If
#'   not, treated as a single implicit dataset.
#' @param clean_fn Function `(raw_uka_subset) -> cleaned` producing a data
#'   frame with `fscore`/`uniprotname`/`LogFC` and `comparison_col` columns
#'   (e.g. [prep_uka()], or `networkScore::prep_uka_paired()` -- whichever
#'   convention the caller's comparison-identifying column uses). Default
#'   `identity`: `raw_uka` is already reshaped.
#' @param comparison_col Name of the column identifying each comparison in
#'   `clean_fn`'s output (e.g. `"Sgroup_contrast"`, `"Sample"`, `"cell_line"`
#'   -- differs by which reshaping convention `clean_fn` uses).
#' @param spec_cutoff,perc_cutoff,b,w,rank_uka_abs Vectors -- every
#'   combination is gridded, not just paired elementwise. `spec_cutoff`/
#'   `perc_cutoff` are passed to [uka_top()] (which also determines
#'   `uka_filt`); `rank_uka_abs` too. `b` and `w` don't affect filtering,
#'   only the later PCSF build (see [generate_kinase_network()] for what
#'   they do), but are gridded here alongside the others for one uniform
#'   mechanism rather than a second, separate expansion step. `b`/`w`
#'   default to `2` (see [PCSF_rand_pg()]).
#' @param ppi_network A data frame with columns `head`, `tail`, `cost`, or a
#'   fully named list of them (e.g. `list(v12_5 = string_12.5,
#'   kins502 = ppi_networkv12_502_kins)`) to grid across more than one
#'   reference network. Defaults to the bundled
#'   [string_12.5].
#' @param dataset_col Name of the dataset-identifying column, if present in
#'   `raw_uka`. Default `"dataset"`.
#'
#' @return A list, one element per `(dataset, condition, spec_cutoff,
#'   perc_cutoff, b, w, rank_uka_abs, ppi_network)` combination, each a list
#'   with `dataset`, `condition`, `spec_cutoff`, `perc_cutoff`, `b`, `w`,
#'   `rank_uka_abs`, `ppi_network_name`, `ppi_network` (the actual data
#'   frame for that name), `uka_filt` (top-hit terminal-node data frame,
#'   from [uka_top()]), `uka_cell_all` (that condition's full, unfiltered
#'   rows -- for permutation building).
#' @export
build_network_grid <- function(raw_uka, clean_fn = identity, comparison_col,
                                spec_cutoff, perc_cutoff, b = 2, w = 2, rank_uka_abs = TRUE,
                                ppi_network = string_12.5, dataset_col = "dataset") {
  ppi_label <- rlang::as_label(rlang::enquo(ppi_network))
  ppi_list <- normalize_ppi_network_list(ppi_network, ppi_label)

  datasets <- if (dataset_col %in% colnames(raw_uka)) {
    split(raw_uka, raw_uka[[dataset_col]])
  } else {
    list(raw_uka)
  }
  dataset_names <- if (dataset_col %in% colnames(raw_uka)) names(datasets) else NA_character_

  cells <- purrr::map(seq_along(datasets), function(i) {
    ds_name <- dataset_names[[i]]
    cleaned <- clean_fn(datasets[[i]])
    conditions <- unique(cleaned[[comparison_col]])
    combos <- tidyr::expand_grid(
      condition = conditions, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff,
      b = b, w = w, rank_uka_abs = rank_uka_abs, ppi_network_name = names(ppi_list)
    )

    purrr::pmap(combos, function(condition, spec_cutoff, perc_cutoff, b, w, rank_uka_abs, ppi_network_name) {
      uka_cell_all <- cleaned[cleaned[[comparison_col]] == condition, ]
      uka_filt <- uka_top(uka_cell_all, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, rank_uka_abs = rank_uka_abs)
      list(
        dataset = ds_name, condition = condition, spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff,
        b = b, w = w, rank_uka_abs = rank_uka_abs, ppi_network_name = ppi_network_name,
        ppi_network = ppi_list[[ppi_network_name]], uka_filt = uka_filt, uka_cell_all = uka_cell_all
      )
    })
  })

  purrr::flatten(cells)
}

#' Create one output folder per distinct grid combination
#'
#' Shared by [run_network_grid()] and `networkScore`'s scoring entry points
#' -- the one place "one [prepare_run_params()] folder per distinct
#' `(spec_cutoff, perc_cutoff, b, w, rank_uka_abs, ppi_network)` combination,
#' given a grid" gets implemented, so it isn't duplicated between plain
#' generation and scoring.
#'
#' @param grid Output of [build_network_grid()].
#' @param prepare_fn [prepare_run_params()], or a package's own extension
#'   of it (e.g. `networkScore::prepare_score_run_params()`, which adds
#'   `nperms_network`/`relative_to` to the naming rule).
#' @param base_labels Named list of label overrides constant across every
#'   combination (typically `uka`/`sens` -- see [capture_params()]'s
#'   `.labels`). A per-combination `ppi_network` label (from the grid
#'   cell's own `ppi_network_name`) is added automatically; don't include
#'   one here.
#' @param ... Passed to `prepare_fn` for every combination (`respath`,
#'   `uka`, and whatever else `prepare_fn` needs) -- `spec_cutoff`/
#'   `perc_cutoff`/`b`/`rank_uka_abs`/`ppi_network` come from each grid
#'   cell itself, not from here.
#'
#' @return A list: `combo_key` (character vector, same length/order as
#'   `grid`, tagging each cell with which combination it belongs to) and
#'   `folders` (named character vector, the folder path for each distinct
#'   `combo_key` value).
#' @export
prepare_grid_folders <- function(grid, prepare_fn = prepare_run_params, base_labels = list(), ...) {
  # sens_perc_cutoff is an optional extra grid dimension, only present on
  # cells run_network_grid() has layered on top of build_network_grid()'s
  # output for the paired path (see its own docs) -- folded into the combo
  # key and forwarded to prepare_fn() whenever it's there, without this
  # shared function needing to know what it means.
  combo_key_of <- function(cell) {
    parts <- c(cell$spec_cutoff, cell$perc_cutoff, cell$b, cell$w, cell$rank_uka_abs, cell$ppi_network_name)
    if (!is.null(cell$sens_perc_cutoff)) parts <- c(parts, cell$sens_perc_cutoff)
    paste(parts, collapse = "||")
  }
  combo_key <- purrr::map_chr(grid, combo_key_of)
  unique_cells <- grid[!duplicated(combo_key)]

  folders <- purrr::map_chr(unique_cells, function(cell) {
    cell_labels <- utils::modifyList(base_labels, list(ppi_network = cell$ppi_network_name))
    extra_args <- if (!is.null(cell$sens_perc_cutoff)) list(sens_perc_cutoff = cell$sens_perc_cutoff) else list()
    do.call(prepare_fn, c(
      list(
        spec_cutoff = cell$spec_cutoff, perc_cutoff = cell$perc_cutoff,
        b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs, ppi_network = cell$ppi_network,
        .labels = cell_labels
      ),
      extra_args,
      list(...)
    ))$respath
  })

  list(
    combo_key = combo_key,
    folders = stats::setNames(folders, purrr::map_chr(unique_cells, combo_key_of))
  )
}

#' Build every network in a parameter grid, in one flattened batch
#'
#' The unified entry point for plain (non-scoring) generation across a
#' `dataset x condition x spec_cutoff x perc_cutoff x b x w x rank_uka_abs x
#' ppi_network` grid -- the caller specifies the grid, not loops. Replaces
#' `generate_networks_for_conditions()` (removed -- see
#' `docs/adr/0006-unified-grid-interface.md`): the single-file,
#' single-combination case is just this function called with length-1
#' vectors, so keeping a separate "simple" function only duplicated logic
#' for no real gain.
#'
#' One [prepare_run_params()] folder is created per distinct combination
#' when `write = TRUE` -- every dataset/condition sharing that combination
#' writes into the same folder, distinguished by file name, not one folder
#' per cell. All cells across the whole grid are still submitted as a
#' single [generate_networks_batch()] call (full flattening, see
#' `docs/adr/0003-full-flattening-parallelization.md`).
#'
#' @inheritParams build_network_grid
#' @param sens Optional raw sensitivity data frame (columns `uniprotname`,
#'   `LogFC` -- the shape [sens_top()] expects); if given, builds paired
#'   networks for every cell. If `NULL` (default), builds kinase-only
#'   networks. One sensitivity profile is shared across every condition
#'   (unlike `uka`, sensitivity here isn't per-condition) -- filtered
#'   per `sens_perc_cutoff` value, not per condition.
#' @param sens_perc_cutoff Required when `sens` is given. A vector -- every
#'   value is gridded (crossed against every other dimension) via
#'   [sens_top()], the same way `perc_cutoff` grids `uka`. There is no
#'   sensitivity-side equivalent of `spec_cutoff`: specificity is a
#'   kinase-activity-only concept.
#' @param sens_balance Passed to [sens_top()] as `balance` for every
#'   `sens_perc_cutoff` value -- fixed, not gridded. Default `TRUE`.
#' @param respath Base output directory for `write = TRUE`.
#' @param write If `TRUE` (default), write each build's `nodes`/`edges`/
#'   `wc_df` CSVs via [prepare_run_params()]-tracked folders.
#' @param max_tasks Refuse to proceed (`stop()`, without building anything)
#'   if the grid expands to more than this many networks -- a safety guard
#'   against an unintentionally huge overnight run, not a time estimate.
#'   Raise explicitly if the size is actually intended. Default 500.
#' @param ... Additional arguments passed to every build, uniform across
#'   the whole grid (not gridded) -- e.g. `n`, `w`, `r`, `mu`, `seed`.
#'
#' @return Same shape as [generate_networks_batch()]: a list, one element
#'   per grid cell, each with `result` and `meta` (`dataset`, `condition`,
#'   `spec_cutoff`, `perc_cutoff`, `b`, `w`, `rank_uka_abs`,
#'   `ppi_network_name`, and, for the paired path, `sens_perc_cutoff`).
#' @export
run_network_grid <- function(raw_uka, clean_fn = identity, comparison_col,
                              spec_cutoff, perc_cutoff, b = 2, w = 2,
                              ppi_network = string_12.5,
                              rank_uka_abs = TRUE, dataset_col = "dataset",
                              sens = NULL, sens_perc_cutoff = NULL, sens_balance = TRUE,
                              respath, write = TRUE, max_tasks = 500, ...) {
  # Captured immediately, before anything else (including the is.null()
  # check below) forces these arguments -- forcing a promise before
  # enquo() silently degrades the captured label to a generic value
  # placeholder (e.g. "<df[,4]>") instead of the caller's actual
  # expression. capture_params(), called several frames down inside
  # prepare_run_params(), would otherwise only see this function's own
  # local variable names anyway (uka/sens), not what our caller wrote.
  uka_label <- rlang::as_label(rlang::enquo(raw_uka))
  ppi_network_label <- rlang::as_label(rlang::enquo(ppi_network))
  sens_label <- rlang::as_label(rlang::enquo(sens))

  paired <- !is.null(sens)
  if (paired && is.null(sens_perc_cutoff)) {
    stop("sens_perc_cutoff is required when sens is given.", call. = FALSE)
  }

  # Normalized here (using our own correctly-captured label above), not
  # left to build_network_grid()'s own normalization -- if we forwarded
  # ppi_network through unnormalized, build_network_grid()'s enquo() would
  # only see our local variable name ("ppi_network"), the same relay
  # problem the label capture above exists to avoid.
  ppi_list <- normalize_ppi_network_list(ppi_network, ppi_network_label)

  base_labels <- list(uka = uka_label)
  if (paired) base_labels$sens <- sens_label

  uka_grid <- build_network_grid(
    raw_uka, clean_fn = clean_fn, comparison_col = comparison_col,
    spec_cutoff = spec_cutoff, perc_cutoff = perc_cutoff, b = b, w = w,
    rank_uka_abs = rank_uka_abs, ppi_network = ppi_list, dataset_col = dataset_col
  )

  # sens_perc_cutoff grids independently of every uka-side dimension --
  # every uka cell is built once per sens_perc_cutoff value, each time
  # against that value's own sens_top() filtering. Layered on top of
  # build_network_grid()'s output here, not inside it: sens has no
  # per-condition structure in this plain-generation path (one profile
  # shared across every uka condition), unlike networkScore's golden-score
  # path, where sensitivity is matched per cell_line -- a different enough
  # shape that build_network_grid() shouldn't need to know about either one.
  grid <- if (paired) {
    sens_filts <- stats::setNames(
      purrr::map(sens_perc_cutoff, ~ sens_top(sens, perc_cutoff = .x, balance = sens_balance)),
      as.character(sens_perc_cutoff)
    )
    purrr::flatten(purrr::map(uka_grid, function(cell) {
      purrr::map(sens_perc_cutoff, function(spc) {
        c(cell, list(sens_perc_cutoff = spc, sens_filt = sens_filts[[as.character(spc)]]))
      })
    }))
  } else {
    uka_grid
  }

  if (length(grid) > max_tasks) {
    stop(
      "run_network_grid() would build ", length(grid), " networks (dataset x condition x ",
      "spec_cutoff x perc_cutoff x b x w x rank_uka_abs x ppi_network",
      if (paired) " x sens_perc_cutoff",
      "), over max_tasks = ", max_tasks,
      ". Review the grid before proceeding -- narrow the grid, or pass a higher max_tasks ",
      "explicitly if this size is actually intended.",
      call. = FALSE
    )
  }

  generate_fn <- if (paired) generate_paired_network else generate_kinase_network

  # When the grid spans a single reference network (the common case), hand
  # it to generate_networks_batch() once as its shared `ppi_network` rather
  # than copying the multi-MB data frame into every task's `args` -- each
  # per-task copy is serialized afresh to every `future` worker (tens of
  # GiB for a grid with many permutation-style tasks). Only a genuine
  # multi-network grid keeps `ppi_network` per-task, and those are kept
  # deliberately small.
  single_ppi_network <- if (length(ppi_list) == 1) ppi_list[[1]] else NULL

  combos <- if (write) {
    prepare_args <- list(respath = respath, uka = raw_uka)
    if (paired) prepare_args$sens <- sens
    do.call(
      prepare_grid_folders,
      c(list(grid, prepare_fn = prepare_run_params, base_labels = base_labels), prepare_args)
    )
  } else {
    NULL
  }

  tasks <- purrr::map(seq_along(grid), function(i) {
    cell <- grid[[i]]
    args <- list(
      uka = cell$uka_filt, condition = cell$condition, spec_cutoff = cell$spec_cutoff,
      b = cell$b, w = cell$w, write = write, ...
    )
    if (is.null(single_ppi_network)) args$ppi_network <- cell$ppi_network
    if (write) args$res.path <- combos$folders[[combos$combo_key[i]]]
    if (paired) args$sens <- cell$sens_filt
    meta <- list(
      dataset = cell$dataset, condition = cell$condition, spec_cutoff = cell$spec_cutoff,
      perc_cutoff = cell$perc_cutoff, b = cell$b, w = cell$w, rank_uka_abs = cell$rank_uka_abs,
      ppi_network_name = cell$ppi_network_name
    )
    if (paired) meta$sens_perc_cutoff <- cell$sens_perc_cutoff
    list(args = args, meta = meta)
  })

  generate_networks_batch(tasks, generate_fn = generate_fn, ppi_network = single_ppi_network)
}
