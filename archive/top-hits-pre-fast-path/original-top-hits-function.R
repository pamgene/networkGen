# Archived pre-fast-path top-hits filter.
#
# Removed from R/top-hits.R -- kept here, not deleted outright, as a
# historical record of the per-omic-layer filtering step used before
# uka_top()/sens_top() (R/data-prep.R) took over as the live path. See
# docs/adr/0007-remove-legacy-top-hits.md for why this was removed from the
# package.
#
# Verbatim copy of top_hits_pg() as it existed in networkGen before removal.
# Not runnable standalone by itself -- needs dplyr's pipe (%>%) and
# dplyr::filter()/mutate()/select() in scope, e.g. library(dplyr), since
# this file is sourced outside the package namespace.
#
# extract_data_percentile() (the multi-omic-layer combiner that called this)
# was NOT archived alongside it -- confirmed to be early exploratory code
# with no real caller ever, safe to delete outright rather than keep as a
# historical reference.

top_hits_pg <- function(df, perc_cutoff, omic_type, spec_cutoff) {
  if (omic_type != "Kinase") {
    df %>%
      dplyr::filter(.data$percentile_score >= perc_cutoff) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rank = seq_len(nrow(.)), prize = .data$percentile_score, type = omic_type) %>%
      dplyr::select("name", "prize", "type", "LogFC")
  } else {
    df %>%
      dplyr::filter(.data$fscore >= spec_cutoff) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rank = seq_len(nrow(.)), prize = .data$percentile_score, type = omic_type) %>%
      dplyr::select("name", "prize", "type", "LogFC")
  }
}
