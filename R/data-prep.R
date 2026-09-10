#' Filter kinase-activity (UKA) data to its top hits
#'
#' Filters `uka` to rows passing `spec_cutoff`, percentile-ranks the
#' survivors by `LogFC`, and keeps those at or above `perc_cutoff`.
#'
#' @param uka Data frame with columns `fscore`, `uniprotname`, `LogFC`.
#' @param spec_cutoff Minimum `fscore` (specificity score) to keep a row.
#' @param rank_uka_abs If `TRUE` (default), rank by `abs(LogFC)` (fast path,
#'   [percentile_score_fast()]); if `FALSE`, rank by signed `LogFC` via
#'   [percentile_score_noabs()].
#' @param perc_cutoff Minimum percentile rank (0-1) to keep a row.
#' @param cs Unused here, kept for call-site compatibility with the historic
#'   `uka_top()` signature (specificity-score column selection happens
#'   upstream, in [clean_uka_to_kinograte()]).
#'
#' @return Data frame with columns `name`, `prize`, `type` ("Kinase"), `LogFC`.
#' @export
uka_top <- function(uka, spec_cutoff, rank_uka_abs = TRUE, perc_cutoff, cs = FALSE) {
  uka <- uka %>% dplyr::filter(.data$fscore >= spec_cutoff)

  if (rank_uka_abs) {
    uka_rank <- percentile_score_fast(uka, uniprotname, LogFC)
  } else {
    uka_rank <- percentile_score_noabs(uka, symbol = uniprotname, metric = LogFC, rank_lowest_highest = FALSE)
  }

  uka_rank %>%
    dplyr::filter(.data$percentile_score >= perc_cutoff) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(type = "Kinase") %>%
    dplyr::select("name", prize = "percentile_score", "type", "LogFC")
}

#' Filter sensitivity data to its top hits
#'
#' Percentile-ranks `sens` by (signed) `LogFC`, most-sensitive-first, and
#' keeps rows at or above `perc_cutoff`.
#'
#' @param sens Data frame with columns `uniprotname`, `LogFC`.
#' @param perc_cutoff Minimum percentile rank (0-1) to keep a row.
#' @param balance If `TRUE` (default), lowers the effective cutoff by 0.2 to
#'   balance the typically smaller pool of sensitivity hits against kinase
#'   hits.
#'
#' @return Data frame with columns `name`, `prize`, `type` ("Sensitivity"), `LogFC`.
#' @export
sens_top <- function(sens, perc_cutoff, balance = TRUE) {
  if (balance) {
    perc_cutoff <- perc_cutoff - 0.2
  }

  sens_rank <- percentile_score_noabs(sens, uniprotname, LogFC, rank_lowest_highest = TRUE)

  sens_rank %>%
    dplyr::filter(.data$percentile_score >= perc_cutoff) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(type = "Sensitivity") %>%
    dplyr::select("name", prize = "percentile_score", "type", "LogFC")
}

#' Relative overlap between two sets of hit names
#'
#' @param uka Character vector of names (e.g. top kinase hits).
#' @param sens Character vector of names (e.g. top sensitivity hits).
#'
#' @return The overlap size divided by `max(length(sens), 1)`.
#' @export
overlap_uka_sens <- function(uka, sens) {
  overlap <- length(intersect(uka, sens))
  overlap / max(length(sens), 1)
}

#' Clean a raw Tercen-exported UKA table into kinograte's expected shape
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param cs `TRUE` to use the per-comparison "Specificity Score"/"Kinase
#'   Statistic" columns (csUKA), `FALSE` to use the "Mean"/"Median" aggregate
#'   variants. `NULL` (default) auto-detects via [detect_csuka()].
#'
#' @return Data frame with columns `Sgroup_contrast`, `uniprotname`, `LogFC`, `fscore`.
#' @export
clean_uka_to_kinograte <- function(uka, cs = NULL) {
  if (is.null(cs)) cs <- detect_csuka(uka)
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"

  uka %>%
    clean_tercen_columns() %>%
    dplyr::select("Sgroup_contrast", "Kinase Name", dplyr::all_of(stat_col), dplyr::all_of(finalscore_col)) %>%
    dplyr::rename(
      "uniprotname" = "Kinase Name", "LogFC" = dplyr::all_of(stat_col),
      "fscore" = dplyr::all_of(finalscore_col)
    ) %>%
    dplyr::distinct()
}

#' Clean a raw UKA table into kinograte's expected shape, with a fixed control
#'
#' Handles both "X vs control" and "control vs X" contrast naming, flipping
#' the statistic's sign when control is on the left so `LogFC` always reads
#' "treatment relative to control."
#'
#' @param uka Raw UKA data frame, Tercen-style dotted column names.
#' @param control Name of the control condition as it appears in `contrast`.
#' @param spec_cutoff Minimum specificity-score to keep a row. Default 0.
#' @param del_cell Unused, kept for call-site compatibility.
#' @param cs `TRUE`/`FALSE` to force per-comparison (csUKA) vs. mean/median
#'   columns; `NULL` (default) auto-detects via [detect_csuka()].
#'
#' @return Data frame with columns `cell_line`, `uniprotname`, `LogFC`, `fscore`.
#' @export
clean_uka_to_kinograte1 <- function(uka, control, spec_cutoff = 0, del_cell = NULL, cs = NULL) {
  if (is.null(cs)) cs <- detect_csuka(uka)
  finalscore_col <- if (cs) "Specificity Score" else "Mean Specificity Score"
  stat_col <- if (cs) "Kinase Statistic" else "Median Kinase Statistic"

  uka %>%
    clean_tercen_columns() %>%
    dplyr::filter(.data[[finalscore_col]] > spec_cutoff) %>%
    dplyr::filter(
      stringr::str_detect(.data$contrast, paste0(" vs ", control)) |
        stringr::str_detect(.data$contrast, paste0(control, " vs "))
    ) %>%
    dplyr::mutate(contrast = gsub(.data$contrast, pattern = "-", replacement = "")) %>%
    dplyr::mutate(
      control_left = stringr::str_detect(.data$contrast, paste0("^", control, " vs ")),
      cell_line = ifelse(.data$control_left,
        stringr::str_replace(.data$contrast, paste0(control, " vs "), ""),
        stringr::str_replace(.data$contrast, paste0(" vs ", control), "")
      ),
      MKS = ifelse(.data$control_left, -.data[[stat_col]], .data[[stat_col]])
    ) %>%
    dplyr::select("cell_line", "Kinase Name", "MKS", dplyr::all_of(finalscore_col)) %>%
    dplyr::rename("uniprotname" = "Kinase Name", "LogFC" = "MKS", "fscore" = dplyr::all_of(finalscore_col)) %>%
    dplyr::mutate(cell_line = gsub(.data$cell_line, pattern = "-", replacement = "")) %>%
    dplyr::distinct()
}

#' Clean sensitivity data into kinograte's expected shape
#'
#' @param sens Raw sensitivity data frame with columns `TARGET_1`, `cell_line`, `LN_IC50`.
#' @param control If given, `LogFC` is computed as `LN_IC50 - LN_IC50` of
#'   `control` for the same target (fold change vs. control). If `NULL` and
#'   `zscore = FALSE`, `sens` is returned as-is (already-computed values, e.g.
#'   median sensitivity). If `NULL` and `zscore = TRUE`, a per-target z-score
#'   across cell lines is computed instead.
#' @param zscore See `control`.
#' @param del_cell Optional character vector of `cell_line` values to drop.
#'
#' @return Data frame with columns `cell_line`, `uniprotname`, `LogFC`.
#' @export
clean_sens_to_kinograte <- function(sens, control, zscore = FALSE, del_cell = NULL) {
  sens_filt <- sens %>% dplyr::rename("uniprotname" = "TARGET_1")

  if (!is.null(control) && !zscore) {
    control_df <- sens_filt %>% dplyr::filter(.data$cell_line == control)
    sens_clean <- sens_filt %>%
      dplyr::left_join(control_df, by = c("uniprotname"), suffix = c("", ".control")) %>%
      dplyr::mutate(LogFC = .data$LN_IC50 - .data$LN_IC50.control) %>%
      dplyr::filter(.data$cell_line != control, !is.na(.data$LogFC)) %>%
      dplyr::select("cell_line", "uniprotname", "LogFC")
  } else if (is.null(control) && !zscore) {
    sens_clean <- sens_filt
  } else if (is.null(control) && zscore) {
    sens_clean <- sens_filt %>%
      tidyr::pivot_wider(id_cols = "cell_line", names_from = "uniprotname", values_from = "LN_IC50") %>%
      tibble::column_to_rownames("cell_line") %>%
      as.matrix() %>%
      scale() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_line") %>%
      tidyr::pivot_longer(cols = -"cell_line", names_to = "uniprotname", values_to = "LogFC")
  }

  if (!is.null(del_cell)) {
    sens_clean <- sens_clean %>% dplyr::filter(!.data$cell_line %in% del_cell)
  }

  sens_clean
}

#' Strip Tercen's dotted column-name prefixes
#'
#' Tercen exports columns as e.g. `some.path.ColumnName`; this keeps only
#' the last dot-separated segment.
#'
#' @param df A data frame with Tercen-style dotted column names.
#'
#' @return `df` with column names shortened to their last dot-segment.
#' @export
clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- stringr::str_split(cols, pattern = "\\.")
  colnames(df) <- sapply(split, utils::tail, 1)
  df
}

#' Detect whether a raw UKA export is a csUKA (per-comparison) export
#'
#' A regular UKA export aggregates its metric columns across comparisons and
#' names them accordingly -- `Mean Specificity Score`, `Median Kinase
#' Statistic`, `Mean Significance Score`, `Median Final score`. A csUKA
#' (per-comparison) export carries one value per comparison and names the
#' same columns without the `Mean`/`Median` prefix -- `Specificity Score`,
#' `Kinase Statistic`, etc. The `clean_uka_to_kinograte*()` functions read
#' different columns for each; this picks which based on what's present, so
#' callers don't have to pass a `cs` flag by hand.
#'
#' @param uka A raw UKA data frame (Tercen-style dotted column names are
#'   fine -- only the last dot-segment is inspected, as [clean_tercen_columns()]
#'   would produce).
#'
#' @return `TRUE` if the columns look like a csUKA export (a bare
#'   `Specificity Score` present and no `Mean Specificity Score`), `FALSE`
#'   otherwise -- including when neither is present, so a malformed input
#'   fails later at the column-selection step with a clearer message rather
#'   than here.
#' @export
detect_csuka <- function(uka) {
  cols <- colnames(uka)
  if (is.null(cols)) {
    return(FALSE)
  }
  bare <- vapply(strsplit(cols, ".", fixed = TRUE), function(x) x[length(x)], character(1))
  "Specificity Score" %in% bare && !("Mean Specificity Score" %in% bare)
}
