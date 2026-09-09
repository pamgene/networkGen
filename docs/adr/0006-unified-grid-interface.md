# Unified grid interface (`run_network_grid()`)

Two functions previously handled different, inconsistent subsets of the
real `dataset x condition x spec_cutoff x perc_cutoff` grid:
`generate_networks_for_conditions()` looped over conditions internally but
required the caller to loop over `spec_cutoff`/`perc_cutoff` externally;
`networkScore`'s `make_golden_score_*()` looped over conditions and
`perc_cutoff` internally but not `spec_cutoff`. No consistent mental model
across the two, and whichever dimension a given function didn't handle
internally, the caller had to write that loop themselves. Neither handled
multiple raw UKA files (many single-comparison exports, or one/few
multi-comparison files) as a first-class case either.

Replaced with three pieces, each with one job:

- **`build_network_grid()`** -- pure grid expansion: `dataset x condition x
  spec_cutoff x perc_cutoff` in, a flat list of resolved cells out, no
  filesystem access, no network builds. The one place this expansion
  happens -- `networkScore`'s scoring entry points build on this same
  function rather than reimplementing their own loop, so plain generation
  and scoring can't drift into different grid semantics.
- **`run_network_grid()`** -- the plain-generation entry point: expands the
  grid via `build_network_grid()`, refuses to proceed past `max_tasks`
  (default 500) without an explicit override (a safety guard against an
  unintentionally huge overnight run, not a time estimate), creates one
  `prepare_run_params()` folder per distinct combination, and submits the
  whole grid as a single `generate_networks_batch()` call (full flattening,
  see `docs/adr/0003-full-flattening-parallelization.md`, unchanged by this
  decision).
- **`load_uka_dataset_files()`** -- parses the real UKA export naming
  convention (`UKA_<PTK|STK>_<number>_<dataset name>.csv`), merges each
  dataset's real PTK/STK pair, and tags every row with its source dataset
  (folder-qualified, since two different experiments can reuse the same
  `<number>_<dataset name>` suffix) -- the one place that convention gets
  parsed, feeding `run_network_grid()`'s dataset-column handling.

Decided: remove `generate_networks_for_conditions()` rather than keep it
alongside `run_network_grid()` -- the single-file, single-combination case
is just `run_network_grid()` called with length-1 vectors, so a separate
"simple" function would only have duplicated logic for no real gain.

Why now: the caller-side inconsistency (which dimension had to be
hand-looped depended on which function you called) was a real, recurring
source of bugs waiting to happen, not a style preference -- fixing it
required one shared grid-expansion mechanism both plain generation and
scoring could build on, not a per-function patch.

(This grid was later generalized further -- `b`, `rank_uka_abs`, and
`ppi_network` became griddable dimensions alongside `spec_cutoff`/
`perc_cutoff`, and a `prepare_grid_folders()` helper was extracted from
`run_network_grid()` for `networkScore` to reuse -- see commit history
rather than this ADR for that later work; the three-piece split and the
reasoning for it above are unchanged.)
