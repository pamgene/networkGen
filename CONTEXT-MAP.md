# Context Map — Network Suite

Three R packages that build, score, and enrich-and-plot PCSF (Prize-Collecting
Steiner Forest) networks from kinase-activity data, replacing the loose
`R/` scripts in `Network_generation`. Each package is its own git repo
under `DevOpti/` (`DevOpti/` itself is just a folder holding many unrelated
projects, not a project or repo of its own).

(Originally planned as four — `networkEnrich` and `networkPlot` were merged
into one `networkPlot`; see `docs/adr/0008-merge-enrich-into-plot.md`.)

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
- [networkPlot](../networkPlot/README.md) — pathway enrichment of a
  generated network (Enrichr → Reactome/KEGG/WikiPathways, with
  Reactome-hierarchy collapsing and cross-comparison reconciliation)
  **and** every rendering of it: the interactive visNetwork HTML and the
  kinase×pathway heatmaps. Works on un-enriched networks too (enrichment is
  a separate, skippable call). *(built — see
  `docs/plans/networkplot-package.md` for the design)*

## Relationships

- **networkScore → networkGen**: calls `networkGen::generate_networks_batch()`
  to build every network it needs (one observed + many permuted networks per
  condition, flattened into a single batch call). `networkGen` has no concept
  of scoring or permutation; `networkScore` owns that entirely and is purely
  a caller of `networkGen`'s batch-build capability.
- **networkPlot → networkGen**: consumes a `networkGen` result's
  `network`/`nodes`/`edges` (the shared network-result contract), plus
  explicit pathway reference tables for the enrichment step. No dependency
  in the other direction, and none on `networkScore`.
- **Network_generation** (the original analysis repo, not part of this
  suite) is the orchestration layer: it will depend on `networkGen` and
  `networkPlot` instead of keeping its own copies of this logic. The
  loose `R/` scripts `networkPlot` replaced (`network_enrichment_and_vis.R`,
  the `plot_kinase_pathway_heatmaps` family in `plotting_functions.R`, and
  the dead per-cluster/topGO enrichment paths) are still to be removed from
  that repo. See that repo's own docs for its migration plan.

## Shared vocabulary

Terms used identically across all three contexts — not redefined per package:

**Network-result object**:
The structured, persistable object every `networkGen` build returns:
`network` (igraph), `nodes`, `edges`, `missing_nodes`, `wc_df` (cluster
assignment), `maintitle`, `params`. The contract both downstream packages
(`networkScore`, `networkPlot`) consume.
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
