#' Build many independent networks in parallel
#'
#' The single parallelization primitive for this whole network suite (see
#' `docs/adr/0003-full-flattening-parallelization.md` -- that decision is
#' substantively about `networkScore`'s design, filed here since suite-wide
#' ADRs are anchored in this package's repo).
#' Has no concept of what the build requests represent -- whether they're
#' many real conditions (plain batch generation) or a permutation task list
#' (`networkScore`'s observed + permuted builds, flattened across every
#' condition being scored) is entirely up to the caller. Each task carries
#' its own `meta`, returned untouched alongside its result, so the caller can
#' regroup results afterward without this function needing to understand the
#' tags.
#'
#' Runs against whatever `future::plan()` is currently active in the calling
#' session -- set one (e.g. `future::plan(future::multisession, workers = 4)`)
#' before calling this for actual parallelism; the default `sequential` plan
#' runs the batch one build at a time. Set up the plan once per session and
#' reuse it across calls, rather than creating/tearing down a plan per batch
#' -- worker startup (including loading `PCSF`'s compiled code, see below)
#' has real cost.
#'
#' `PCSF`'s compiled `.Call("_PCSF_call_sr", ...)` requires the `PCSF`
#' namespace to be loaded in each worker process -- a real requirement of
#' compiled/Rcpp-style packages under socket-cluster parallelism, not an
#' instance of global-environment coupling. This is declared explicitly via
#' `furrr_options(packages = "PCSF")` below rather than relying on `future`'s
#' automatic package re-attachment, since this exact failure mode
#' (`"_PCSF_call_sr" not available for .Call()`) has been hit before in this
#' project's history.
#'
#' @param tasks A list of tasks, each a list with:
#'   - `args`: named list of arguments to pass to `generate_fn` for this
#'     build (e.g. `uka`, `sens`, `spec_cutoff`, `b`, `condition`).
#'   - `meta`: arbitrary data to carry through untouched, for the caller's
#'     own bookkeeping (e.g. `list(condition = "A_vs_B", role = "permutation", perm_index = 3)`).
#' @param generate_fn The single-network build function to call for every
#'   task -- [generate_paired_network()] or [generate_kinase_network()].
#' @param ppi_network Data frame with columns `head`, `tail`, `cost` -- the
#'   PPI network used for any task that doesn't specify its own via
#'   `task$args$ppi_network` (the common case: one network shared across
#'   the whole batch). A task carrying its own `ppi_network` -- e.g. when
#'   gridding across more than one reference network -- overrides this.
#' @param extra_args Named list of additional arguments merged into every
#'   task's `args` (task-specific values in `args` take precedence over
#'   same-named entries here). Typically constant across the whole batch,
#'   e.g. `list(write = FALSE)`.
#' @param progress If `TRUE`, show a `furrr` progress bar.
#'
#' @return A list, same length and order as `tasks`, each element a list of
#'   `result` (the `"networkGen_result"` from `generate_fn`, or `NULL` if
#'   that build failed) and `meta` (unchanged from the input task).
#' @export
generate_networks_batch <- function(tasks, generate_fn, ppi_network = NULL, extra_args = list(), progress = FALSE) {
  furrr::future_map(
    tasks,
    function(task) {
      build_args <- utils::modifyList(extra_args, task$args)
      if (is.null(build_args$ppi_network)) build_args$ppi_network <- ppi_network
      list(result = do.call(generate_fn, build_args), meta = task$meta)
    },
    .options = furrr::furrr_options(seed = TRUE, packages = "PCSF"),
    .progress = progress
  )
}

