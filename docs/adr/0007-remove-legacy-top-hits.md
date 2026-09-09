# Remove the legacy top-hits filter and multi-omic combiner

`networkGen` carried two parallel "filter data to its top hits" systems:
`uka_top()`/`sens_top()` (`R/data-prep.R`), the live path used by every real
caller (`run_network_grid()`, every `networkScore` entry point), and
`top_hits_pg()`/`extract_data_percentile()` (`R/top-hits.R`), which nothing
in the package actually called -- confirmed by grep, not inference. Despite
having no real caller, `top_hits_pg()`/`extract_data_percentile()` were
cross-linked from `kinograte_pg_pcsf()`'s own `@param df` docs as if they
were a normal, live way to build its input, misleading a reader into
thinking both systems were real alternatives.

The two orphaned functions aren't equivalent to each other, and got
different treatment:

- **`top_hits_pg()`** filters one omic layer to its top hits by percentile
  (non-Kinase) or specificity (Kinase) -- the direct predecessor of
  `uka_top()`/`sens_top()`, confirmed by tracing what each calls
  internally: `uka_top()` ranks via `percentile_score_fast()`
  (`data.table`-backed), while `top_hits_pg()`'s non-Kinase branch (via
  `extract_data_percentile()`) still used the older, plain `percentile_score()`
  (`dplyr`/`percent_rank()`-based). Same relationship as the already-archived
  igraph-based PCSF path superseded by the dataframe-based fast path (see
  `docs/adr/0005-remove-original-pcsf-path.md`) -- so archived the same way:
  kept as a plain, non-package file under `archive/top-hits-pre-fast-path/`
  (excluded from the package build via `.Rbuildignore`'s existing `^archive$`
  rule, not discovered by `testthat`), rather than deleted outright.
- **`extract_data_percentile()`** (the multi-omic-layer combiner that called
  `top_hits_pg()` internally, merging several omic layers -- RNA, protein,
  kinase, sensitivity -- into one terminal-node set) was confirmed to be
  early exploratory code from before the current pipeline existed, with no
  real caller ever and no equivalent-but-faster successor to preserve as a
  historical record. Deleted outright, not archived.

Also removed: `tests/testthat/test-top-hits.R` (tested only these two
functions, now gone), and the misleading `@param df` cross-reference in
`kinograte_pg_pcsf()`'s docs, repointed to `uka_top()`/`sens_top()`.

Not touched: `percentile_score()` (the plain, non-`data.table` percentile
ranker `top_hits_pg()` used internally) stays -- it's an independently
exported, independently tested general-purpose utility
(`test-percentile-score.R`), not something that existed only to serve
`top_hits_pg()`, even though it now has no internal caller either.
