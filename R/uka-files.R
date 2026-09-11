#' Load and tag raw UKA files by dataset, merging PTK/STK pairs
#'
#' The one place file-path/filename knowledge lives in this package --
#' everything downstream ([build_network_grid()], the cleaning functions)
#' stays filesystem-blind and works on plain data frames.
#'
#' Real UKA exports are named `UKA_<PTK|STK>_<number>_<dataset name>.csv`
#' (the dataset name itself may contain `_`; it always comes after the
#' number and before `.csv`). A PTK file and an STK file are the same
#' dataset when they share both a folder and a `<number>_<dataset name>`
#' suffix -- these get row-bound together, tagged `kinase_type` for
#' provenance. The folder is part of the dataset identity: two different
#' experiments can reuse the same `<number>_<dataset name>` (e.g. both named
#' `01_TvsC`), and are not the same dataset just because the suffix matches.
#'
#' A path that doesn't match this naming convention falls back to using its
#' own file name (without extension) as a standalone dataset -- no PTK/STK
#' merging attempted, since there's nothing to merge.
#'
#' Each file's columns are cleaned ([clean_tercen_columns()]) *before* the
#' files are row-bound. Real UKA exports carry a Tercen app-instance number
#' in their column prefixes (`csUKA_app0.UKA0.Sgroup_contrast` in one file,
#' `csUKA_app1.UKA00.Sgroup_contrast` in another) that differs from file to
#' file. Binding the *raw*, still-prefixed frames first -- and cleaning only
#' afterward -- used to produce one real column per distinct prefix variant
#' plus an all-NA leftover column for every other file's variant; cleaning
#' then collapsed all of them onto the same short name (e.g.
#' `Sgroup_contrast`), and a downstream `select()` could silently grab one
#' of the all-NA duplicates instead of the real column for a given dataset.
#' Cleaning per file first means every file already shares the one
#' canonical name by the time they're bound, so no duplicate columns are
#' ever created.
#'
#' @param paths Character vector of raw UKA CSV file paths.
#'
#' @return One combined data frame (Tercen prefixes already stripped -- see
#'   [clean_tercen_columns()]), with `dataset` and `kinase_type` (`"PTK"`,
#'   `"STK"`, or `NA` for non-matching files) columns added. Pass to
#'   [build_network_grid()].
#' @export
load_uka_dataset_files <- function(paths) {
  pattern <- "^UKA_(PTK|STK)_[0-9]+_(.+)\\.csv$"

  tagged <- purrr::map(paths, function(path) {
    fname <- basename(path)
    m <- regmatches(fname, regexec(pattern, fname))[[1]]

    if (length(m) == 3) {
      suffix <- sub("^UKA_(PTK|STK)_", "", tools::file_path_sans_ext(fname))
      dataset <- paste0(basename(dirname(path)), "/", suffix)
      kinase_type <- m[2]
    } else {
      dataset <- tools::file_path_sans_ext(fname)
      kinase_type <- NA_character_
    }

    df <- readr::read_csv(path, show_col_types = FALSE)
    df <- clean_tercen_columns(df)
    df$dataset <- dataset
    df$kinase_type <- kinase_type
    df
  })

  dplyr::bind_rows(tagged)
}
