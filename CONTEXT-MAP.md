# Context Map — Network Suite

Four R packages that build, score, enrich, and plot PCSF (Prize-Collecting
Steiner Forest) networks from kinase-activity data, replacing the loose
`R/` scripts in `Network_generation`. Each package is its own git repo
under `DevOpti/` (`DevOpti/` itself is just a folder holding many unrelated
projects, not a project or repo of its own).

This file, and `docs/adr/`, live in `networkGen`'s repo (the foundation
package) rather than at the `DevOpti` root or in a separate docs-only repo —
they describe the whole suite, not just `networkGen`, so a decision recorded
here or in `docs/adr/` may be substantively about a different package (e.g.
`networkScore`) even though it's filed here.

## Contexts

- [networkGen](./CONTEXT.md) — builds a PCSF network from kinase-activity
  data (optionally paired with sensitivity data). The foundation; no
  scoring, enrichment, or plotting logic. *(built)*
- [networkScore](../networkScore/CONTEXT.md) — computes golden score
  (paired kinase+sensitivity) and kinase-only score via permutation
  testing. *(built)*
- networkEnrich — Reactome pathway enrichment of a generated network, plus
  visualizing that enrichment output (pathway heatmaps, hierarchy trees).
  *(designed, not yet built)*
- networkPlot — the interactive network HTML visualization only.
  *(designed, not yet built)*

## Relationships

- **networkScore → networkGen**: calls `networkGen::generate_networks_batch()`
  to build every network it needs (one observed + many permuted networks per
  condition, flattened into a single batch call). `networkGen` has no concept
  of scoring or permutation; `networkScore` owns that entirely and is purely
  a caller of `networkGen`'s batch-build capability.
- **networkEnrich → networkGen**: consumes a `networkGen` result's
  `network`/`nodes` (the shared network-result contract) plus explicit
  pathway reference tables. No dependency in the other direction.
- **networkPlot → networkGen**: consumes a `networkGen` result's
  `nodes`/`edges`. Optionally consumes a nodes-with-pathway data frame (the
  shape `networkEnrich` produces) but has no package dependency on
  `networkEnrich` — any data frame with the right columns works, so it stays
  usable for un-enriched networks too.
- **Network_generation** (the original analysis repo, not part of this
  suite) is the orchestration layer: it will depend on `networkGen` (and
  later `networkEnrich`/`networkPlot`) instead of keeping its own copies of
  this logic. See that repo's own docs for its migration plan.

## Shared vocabulary

Terms used identically across all four contexts — not redefined per package:

**Network-result object**:
The structured, persistable object every `networkGen` build returns:
`network` (igraph), `nodes`, `edges`, `missing_nodes`, `wc_df` (cluster
assignment), `maintitle`, `params`. The contract all three downstream
packages consume.
_Avoid_: "kinograte_res" (old ad hoc variable name)

**Condition**:
One comparison/sample within a kinase-activity (UKA) dataset — e.g. one
cell line's treated-vs-control comparison. A dataset may contain one
condition (one file, one comparison) or many (one file, many comparisons);
the packages in this suite treat every condition as an independent build
input regardless of which file it came from.
_Avoid_: cell, sample, comparison (used loosely elsewhere in older code)

**PPI network**:
The protein-protein interaction network (edges with `head`/`tail`/`cost`
columns) that a PCSF build searches over. Supplied by the caller — never a
default/global value.
_Avoid_: ppi_networkv12 (a specific loaded object's name, not a concept)
