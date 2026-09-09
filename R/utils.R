#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`

#' @importFrom rlang .data
#' @importFrom data.table :=
NULL

# Bare column-name symbols passed positionally into tidy-eval helpers
# (percentile_score()/percentile_score_noabs()/percentile_score_fast()'s
# `symbol`/`metric` args), plus the magrittr `.` pronoun and data.table's
# `.N` pronoun -- not real global variables, but R CMD check's static
# analysis can't tell the difference.
#
# ppi_networkv12 is the bundled reference dataset (see data.R), used as the
# default `ppi_network` argument in generate_paired_network()/
# generate_kinase_network() -- lazy-loaded package data, not a real global,
# but again indistinguishable to static analysis.
utils::globalVariables(c("uniprotname", "LogFC", ".", ".N", "ppi_networkv12"))

#' Build the `_spec<cutoff>` filename suffix, omitted when there's no filter
#'
#' `spec_cutoff = 0` means "no specificity filtering" -- encoding it in every
#' output filename as `_spec0` would just be noise, so it's dropped entirely
#' in that case.
#'
#' @param spec_cutoff The specificity cutoff used for this run.
#'
#' Vectorized: safe to use on a single value or a vector of them.
#'
#' @return `""` where `spec_cutoff` is `0`; `"_spec<spec_cutoff>"` elsewhere.
#' @keywords internal
spec_suffix <- function(spec_cutoff) {
  ifelse(spec_cutoff == 0, "", paste0("_spec", spec_cutoff))
}
