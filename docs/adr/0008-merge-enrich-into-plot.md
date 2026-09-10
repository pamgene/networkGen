# Merge `networkEnrich` into `networkPlot` — three packages, not four

Supersedes the package count in
[0001](./0001-four-package-split.md): the suite is now **three** packages
(`networkGen`, `networkScore`, `networkPlot`), not four.

ADR 0001 split the `Network_generation` logic into `networkGen` (foundation),
`networkScore`, `networkEnrich` (Reactome/Enrichr pathway enrichment of a
built network, plus pathway heatmaps and hierarchy trees), and `networkPlot`
(the interactive visNetwork HTML). `networkGen` and `networkScore` are
built; `networkEnrich` and `networkPlot` were designed-not-built.

While planning `networkPlot` it became clear the enrichment output has no
consumer other than the plot: the pathway table exists to add a `pathway`
column to `nodes` (which colours and filters the HTML) and to feed the
kinase×pathway heatmaps, which are themselves a visualization of that same
network. There is one linear workflow — build → enrich → plot — with no
point where a caller wants enrichment without going on to plot, or plots
from enrichment produced by a different tool.

Decided: fold `networkEnrich` into `networkPlot`. `networkPlot` owns both
the enrichment (Enrichr call, Reactome-hierarchy collapse, cross-comparison
pathway reconciliation, WikiPathways ontology filtering) and every rendering
of it (the interactive HTML, the pathway heatmaps). It still consumes only a
`networkGen` result plus explicit pathway reference tables, and it stays
usable for un-enriched networks — `plot_network()` works with or without a
`pathway` column, and enrichment is a separate call the caller may skip.

Why: a package seam across a single linear workflow with no independent
consumer of the intermediate is pure overhead — one more repo to version and
keep contract-compatible, for a boundary nothing crosses. The
generation/scoring seams in 0001 are real (`Network_generation` wants
generation without scoring's parallel/permutation stack; other consumers
want one piece); the enrich/plot seam is not.

Unchanged from 0001: `networkGen` is still the foundation, `networkScore`
still depends only on `networkGen`, and `networkPlot` still depends only on
`networkGen` (for the network-result contract), not on `networkScore`.
