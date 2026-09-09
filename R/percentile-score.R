#' Percentile-rank a data frame by a metric
#'
#' Ranks rows by the absolute value of `metric`, producing a `percentile_score`
#' column in \[0, 1\] (1 = strongest). Copied in from the (now unused,
#' dormant since 2022) `kinograte` package's `percentile_rank()` -- see
#' `vignette("networkGen")` / project ADRs for why this isn't a live
#' dependency. Renamed from `score` to `percentile_score` to avoid colliding
#' with networkScore's unrelated "significance score" concept.
#'
#' @param df Data frame with the columns named by `symbol` and `metric`.
#' @param symbol Column (unquoted) with the row identifier (e.g. protein name).
#' @param metric Column (unquoted) to rank by.
#' @param desc If `TRUE`, rank descending by `-abs(metric)` instead of `abs(metric)`.
#'
#' @return `df`, deduplicated by `symbol` (renamed to `name`), with a
#'   `percentile_score` column, sorted by that column descending.
#' @export
percentile_score <- function(df, symbol, metric, desc = FALSE) {
  df <- df %>%
    dplyr::filter(!is.na({{ metric }})) %>%
    dplyr::rename(name = {{ symbol }}) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)

  if (desc) {
    df <- dplyr::mutate(df, percentile_score = dplyr::percent_rank(dplyr::desc(base::abs({{ metric }}))))
  } else {
    df <- dplyr::mutate(df, percentile_score = dplyr::percent_rank(base::abs({{ metric }})))
  }

  dplyr::arrange(df, dplyr::desc(.data$percentile_score))
}

#' Percentile-rank a data frame by a metric, without taking the absolute value
#'
#' Like [percentile_score()], but ranks directly by `metric` (not its
#' absolute value) -- used where sign matters, e.g. sensitivity data where
#' "most sensitive" means most negative.
#'
#' @param df Data frame with the columns named by `symbol` and `metric`.
#' @param symbol Column (unquoted) with the row identifier.
#' @param metric Column (unquoted) to rank by.
#' @param rank_lowest_highest If `TRUE`, the lowest values of `metric` get the
#'   highest `percentile_score` (e.g. most negative LogFC = most sensitive).
#'
#' @return `df`, deduplicated by `symbol` (renamed to `name`), with a
#'   `percentile_score` column, sorted by that column descending.
#' @export
percentile_score_noabs <- function(df, symbol, metric, rank_lowest_highest = FALSE) {
  df <- df %>%
    dplyr::filter(!is.na({{ metric }})) %>%
    dplyr::rename(name = {{ symbol }}) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)

  if (rank_lowest_highest) {
    df <- dplyr::mutate(df, percentile_score = dplyr::percent_rank(dplyr::desc({{ metric }})))
  } else {
    df <- dplyr::mutate(df, percentile_score = dplyr::percent_rank({{ metric }}))
  }

  dplyr::arrange(df, dplyr::desc(.data$percentile_score))
}

#' Percentile-rank a data frame by a metric (data.table-backed, fast path)
#'
#' Same semantics as [percentile_score()] (ranks by `abs(metric)`), but
#' implemented with `data.table` for speed on large inputs (e.g. permutation
#' loops in networkScore, called once per permutation).
#'
#' @param df Data frame with the columns named by `symbol` and `metric`.
#' @param symbol Column (unquoted) with the row identifier.
#' @param metric Column (unquoted) to rank by.
#' @param desc If `TRUE`, rank descending by `-abs(metric)` instead of `abs(metric)`.
#'
#' @return A `data.table`, deduplicated by `symbol` (renamed to `name`), with
#'   a `percentile_score` column, sorted by that column descending.
#' @export
percentile_score_fast <- function(df, symbol, metric, desc = FALSE) {
  dt <- data.table::as.data.table(df)

  sym_col <- deparse(substitute(symbol))
  met_col <- deparse(substitute(metric))

  dt <- dt[!is.na(get(met_col))]
  data.table::setnames(dt, sym_col, "name", skip_absent = TRUE)
  dt <- unique(dt, by = "name")

  abs_metric <- abs(dt[[met_col]])

  if (desc) {
    dt[, percentile_score := data.table::frank(-abs_metric, ties.method = "average") / .N]
  } else {
    dt[, percentile_score := data.table::frank(abs_metric, ties.method = "average") / .N]
  }

  data.table::setorder(dt, -percentile_score)
  dt[]
}
