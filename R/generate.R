#' Build a network-result object
#'
#' Internal constructor for the shared network-result contract every
#' `networkGen` build returns: `network`, `nodes`, `edges`, `missing_nodes`,
#' `wc_df`, `maintitle`, `params`. Normalizes the `missing`/`missing_nodes`
#' field-name inconsistency present in [kinograte_pg_pcsf()]'s return shape
#' (differs between its `cluster = TRUE`/`cluster = FALSE` branches).
#'
#' @param raw The list returned by [kinograte_pg_pcsf()], or `NULL` if the
#'   build failed.
#' @param params A named list of the parameters used for this build.
#'
#' @return A list of class `"networkGen_result"`, or `NULL` if `raw` is `NULL`.
#' @keywords internal
new_network_result <- function(raw, params) {
  if (is.null(raw)) {
    return(NULL)
  }

  missing_nodes <- if (!is.null(raw$missing_nodes)) raw$missing_nodes else raw$missing

  structure(
    list(
      network = raw$network,
      nodes = raw$nodes,
      edges = raw$edges,
      missing_nodes = missing_nodes,
      wc_df = raw$wc_df,
      maintitle = raw$maintitle,
      params = params
    ),
    class = "networkGen_result"
  )
}

#' Generate a PCSF network from paired kinase-activity + sensitivity data
#'
#' The pure generation step for the "golden" (paired) analysis: combines
#' pre-filtered kinase and sensitivity terminal-node data frames (e.g. from
#' [uka_top()] / [sens_top()]) and builds one PCSF network. Computes no
#' topology statistics (that's `networkScore`'s job) and does no pathway
#' enrichment or plotting (that's `networkEnrich`/`networkPlot`'s job).
#'
#' @param uka Pre-filtered kinase terminal-node data frame (columns `name`,
#'   `prize`, `type`, `LogFC`), e.g. from [uka_top()].
#' @param sens Pre-filtered sensitivity terminal-node data frame, same shape,
#'   e.g. from [sens_top()].
#' @param ppi_network Data frame with columns `head`, `tail`, `cost` -- the
#'   PPI network to search over. Defaults to the bundled [ppi_networkv12]
#'   reference network (lazy-loaded from the package's own data, not a
#'   global variable -- pass a different data frame explicitly, e.g.
#'   [ppi_networkv12_filt] or [ppi_networkv12_502_kins], to use another one).
#' @param spec_cutoff Specificity cutoff, recorded in `maintitle`/output file
#'   names.
#' @param b,w The two main PCSF cost knobs -- `b` multiplies node prizes
#'   (higher `b` = bigger network), `w` is the cost of connecting a terminal
#'   straight to the root (higher `w` = isolated terminals dropped). See
#'   [PCSF_rand_pg()] for the full picture. Both default 2.
#' @param condition Condition/comparison label, recorded in `maintitle`/output
#'   file names.
#' @param res.path Output folder, used only when `write = TRUE`.
#' @param art_nodes,art_lfc Optional artificial terminal nodes to add before
#'   building: `art_nodes` is a character vector of names, `art_lfc` their
#'   `LogFC` values (prize fixed at 1, type `"Artificial"`).
#' @param write If `TRUE`, write `nodes_*.csv`/`missing_nodes_*.csv` to
#'   `res.path`. Default `FALSE`.
#' @param ... Passed through to [PCSF_rand_pg()] (`n`, `r`, `mu`, `seed`).
#'
#' @return A `"networkGen_result"` object (see [new_network_result()]), or
#'   `NULL` if PCSF could not find a subnetwork for these inputs.
#' @export
generate_paired_network <- function(uka, sens, ppi_network = ppi_networkv12, spec_cutoff, b = 2, w = 2,
                                     condition = NULL, res.path = NULL,
                                     art_nodes = NULL, art_lfc = NULL,
                                     write = FALSE, ...) {
  if (!is.null(art_nodes)) {
    art_df <- data.frame(name = art_nodes, prize = 1, type = "Artificial", LogFC = art_lfc)
    uka <- dplyr::bind_rows(uka, art_df)
  }

  combined_df <- dplyr::bind_rows(list(Sensitivity = sens, Kinase = uka), .id = "type") %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarize(
      type = paste0(.data$type, collapse = "-"),
      prize = mean(.data$prize),
      LogFC_all = paste0(round(.data$LogFC, 2), collapse = ", ")
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(uka[, c("name", "LogFC")], by = "name") %>%
    dplyr::left_join(sens[, c("name", "LogFC")] %>% dplyr::rename(LogFC_s = "LogFC"), by = "name") %>%
    dplyr::mutate(LogFC = round(dplyr::coalesce(.data$LogFC, .data$LogFC_s), 2)) %>%
    dplyr::select(-"LogFC_s")

  sumdf <- combined_df %>% dplyr::group_by(.data$type) %>% dplyr::summarize(n_nodes = dplyr::n())
  print(sumdf)

  raw <- kinograte_pg_pcsf(
    df = combined_df, ppi_network = ppi_network, spec_cutoff = spec_cutoff,
    res.path = res.path, condition = condition, cluster = TRUE, b = b, w = w, write = write, ...
  )

  new_network_result(raw, params = list(
    spec_cutoff = spec_cutoff, b = b, w = w, condition = condition
  ))
}

#' Generate a PCSF network from kinase-activity data only
#'
#' The pure generation step for the kinase-only analysis: builds one PCSF
#' network from a pre-filtered kinase terminal-node data frame (e.g. from
#' [uka_top()]). Computes no topology statistics (that's `networkScore`'s
#' job) and does no pathway enrichment or plotting (that's
#' `networkEnrich`/`networkPlot`'s job).
#'
#' @param uka Pre-filtered kinase terminal-node data frame (columns `name`,
#'   `prize`, `type`, `LogFC`), e.g. from [uka_top()].
#' @inheritParams generate_paired_network
#'
#' @return A `"networkGen_result"` object (see [new_network_result()]), or
#'   `NULL` if PCSF could not find a subnetwork for these inputs.
#' @export
generate_kinase_network <- function(uka, ppi_network = ppi_networkv12, spec_cutoff, b = 2, w = 2,
                                     condition = NULL, res.path = NULL,
                                     art_nodes = NULL, art_lfc = NULL,
                                     write = FALSE, ...) {
  if (!is.null(art_nodes)) {
    art_df <- data.frame(name = art_nodes, prize = 1, type = "Artificial", LogFC = art_lfc)
    uka <- dplyr::bind_rows(uka, art_df)
  }

  raw <- tryCatch(
    {
      kinograte_pg_pcsf(
        df = uka, ppi_network = ppi_network, spec_cutoff = spec_cutoff,
        res.path = res.path, condition = condition, cluster = TRUE, b = b, w = w, write = write, ...
      )
    },
    error = function(e) {
      cat("\nCan't generate network for ", condition, " with spec cutoff ", spec_cutoff, ": ", e$message, "\n")
      NULL
    }
  )

  new_network_result(raw, params = list(
    spec_cutoff = spec_cutoff, b = b, w = w, condition = condition
  ))
}
