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
#' @param paths Character vector of raw UKA CSV file paths.
#'
#' @return One combined data frame (raw Tercen-style columns, as read),
#'   with `dataset` and `kinase_type` (`"PTK"`, `"STK"`, or `NA` for
#'   non-matching files) columns added. Pass to [build_network_grid()].
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
    df$dataset <- dataset
    df$kinase_type <- kinase_type
    df
  })

  dplyr::bind_rows(tagged)
}
