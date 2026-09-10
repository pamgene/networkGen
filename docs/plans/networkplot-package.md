# Plan — `networkPlot` package (HTML visualization + pathway enrichment)

## 1. Goal

Stand up `networkPlot`: one package that takes a `networkGen` result and
produces both

- the **interactive HTML network** (visNetwork), and
- the **pathway enrichment** of that network (Enrichr → Reactome/KEGG/
  WikiPathways), including the enrichment-driven pathway heatmaps and the
  per-node pathway annotation that colours/filters the HTML.

Source today: `Network_generation/R/network_enrichment_and_vis.R` (both the
`visualize_network_pg()` plotter and the `do_network_enrichment()` family)
plus the `plot_kinase_pathway_heatmaps()` family in
`Network_generation/R/plotting_functions.R`, wired together by
`make_network_and_stats*()` in `kinograte_PG.R`.

Three things to settle while planning:

1. **Prune the enrichment code** carried from the `kinograte`/`PCSF`
   lineage. Specific decisions in §3.1.
2. **Lift plot styling** (node shapes per type, colours, highlight rules)
   out of the function body into a config object. §3.2.
3. **Design for a future feature**: edge boldness = link strength. §3.3.

## 2. Overview

### One package, two halves

Enrichment and visualization ship together as `networkPlot`, because in
practice the enrichment output (a per-node `pathway` column, pathway↔gene
tables) only exists to drive the plot and its heatmaps — splitting them just
adds a package boundary across one workflow. Recorded in
`docs/adr/0008-merge-enrich-into-plot.md`; `CONTEXT-MAP.md` and the
superseded-note on `docs/adr/0001` are updated to match.

### Public surface

| Function | Does | Pure? |
|---|---|---|
| `enrich_network(network, nodes, databases = "Reactome_Pathways_2024", min_hits = 2, ontology_filter = TRUE)` | one Enrichr call on the network's kinase nodes → filtered, Reactome-hierarchy-collapsed pathway table; joins a `pathway` column onto `nodes` | returns data frames, no writes |
| `reconcile_pathways(pathway_tables)` | re-pick one pathway per gene-set across several comparisons so heatmaps line up | pure |
| `plot_network(nodes, edges, title, clusters = NULL, colour_by = "pathway", theme = networkplot_theme())` | build the visNetwork htmlwidget | pure |
| `save_network_html(widget, path)` | the `saveWidget` → `visSave` fallback dance, one place | write |
| `plot_pathway_heatmaps(pathway_tables, nodes_tables, out_dir, min_kinases = 3)` | per-comparison + combined kinase×pathway `ComplexHeatmap` PNGs | write |
| `networkplot_theme(...)` | config object holding every styling knob, defaults = today's look | — |

`make_network_and_stats*()` do **not** move in whole — they are monolith
orchestration (build + enrich + plot + disk I/O). The monolith keeps a thin
orchestrator that calls `networkGen::generate_*` → `networkPlot::enrich_network`
→ `networkPlot::plot_network` → `save_network_html`.

### Link-strength feature — `w` is not it

`w` is one scalar per build (cost of parking a terminal on the artificial
root); it cannot vary per edge. Two real per-edge quantities:

- `edges$weight` from `PCSF_rand_pg()` = how many of the `n` randomized runs
  kept the edge (a robustness count). Already reaches the plotter untouched
  — a theme-driven map to visNetwork's `width` needs **no `networkGen`
  change**.
- the **PPI edge `cost`** (biological interaction confidence) — the
  preferred signal, but it is dropped by the PCSF solver and must be
  reattached at the build. **Small additive `networkGen` change.** See §3.3.

## 3. Details

### 3.1 Enrichment code — decisions

| Item | Decision | Notes |
|---|---|---|
| **topGO / `mode = 1`** path in `call_enr_pg()` | **DELETE** | Offline GO enrichment via `topGO` + `org.Hs.eg.db`. Never needed; `mode` is never passed by any caller. Removing it drops a heavy Bioconductor dependency. |
| **`reactome_pw_hierarchy_vis()`** | **DELETE** | `ggraph` tree of the Reactome hierarchy. Never called (one commented call site); has a live `browser()` mid-body. |
| **`augment_by_threshold_steps()`** (`kinograte_PG.R`) | **DELETE, don't migrate** | Grows a PCSF result outward along low-cost edges. Stale copy — the maintained version lives in a different repo; this one references an undefined `core_nodes0` and would error. |
| **WikiPathways: `filter_wps_by_ontology()` + its SPARQL query** | **KEEP** | Filters WikiPathways hits to signalling/metabolic/regulatory ontology tags, and drops disease/drug/syndrome pathways, via the live `sparql.wikipathways.org` endpoint. Currently only wired into the *per-cluster* path (§ below) — **must be rewired** so the live enrichment path can use it whenever `WikiPathways_2024_Human` is among `databases`. |
| **Per-cluster enrichment** (`network_enrichment_pg` → `enrichment_analysis_pg` → `call_enr_pg` mode 0) | **DELETE** — explained below | Unreachable today; not being kept as an option. Salvage `filter_wps_by_ontology()` out of it first. |
| **The ~50-line commented block** at the tail of `do_network_enrichment()` | **DELETE** — explained below | Superseded post-processing recipe. |
| **`add_reactome_hierarchy()` trailing line** | **KEEP — it is deliberate** | Dropping top-level / shallow (`hierarchy_n <= 1`) Reactome pathways is intended: they are too generic to be useful in the result. Make it an explicit `return()` with a comment so it is not re-flagged as a debug leftover. |

#### What "per-cluster enrichment" is

Two ways `do_network_enrichment()` can enrich a network:

- **Simple path** (`per_cluster = FALSE`, the default — and the only one any
  caller actually uses): take the network's kinase nodes as *one* gene
  list, make *one* Enrichr call, keep pathways with more than `min_n_hits`
  overlapping genes.

- **Per-cluster path** (`per_cluster = TRUE`):
  1. `igraph::cluster_edge_betweenness(network)` splits the network into
     communities — it repeatedly removes the edge with the highest
     betweenness (the edges that act as bridges between dense regions),
     and the pieces that fall apart are the clusters.
  2. Each community's gene list is sent to Enrichr **separately**, so every
     community gets its own pathway table, tagged with a `Cluster` number.
  3. `enrichment_analysis_pg()` also stamps `V(subnet)$group <-
     clusters$membership` and builds a per-cluster HTML table (for a
     tooltip when you hover a node), and derives a `"PCSFe"` class — this
     is lifted almost verbatim from the `PCSF` package's own
     `enrichment_analysis()`.
  4. Back in `do_network_enrichment()`, the per-cluster tables are pooled:
     group by pathway `Term` across clusters, **sum** the overlapping-gene
     counts, keep pathways with ≥ `min_n_hits` total.

  Rationale it was built on: different network modules = different biology,
  so enrich them separately and you can annotate each module. In practice
  the whole-network call was enough, the per-module annotation was dropped,
  and per-cluster is slower (one Enrichr round-trip per community) — so it
  went unused. The one useful thing only it currently calls is
  `filter_wps_by_ontology()`.

  **Decision: delete it.** Remove `network_enrichment_pg()`,
  `enrichment_analysis_pg()`, `call_enr_pg()` (the edge-betweenness +
  per-cluster-Enrichr machinery) and the `per_cluster` argument of
  `do_network_enrichment()` entirely. Move `filter_wps_by_ontology()` out
  first (it is standalone — takes a pathway data frame, returns it filtered)
  and call it from the live enrichment path whenever WikiPathways is among
  the requested databases.

#### What the ~50-line commented block was

It sits at the end of `do_network_enrichment()`, all commented out. It was
an **older way to de-duplicate the per-cluster pathway results**:

1. keep only the pathways that survived filtering
   (`enrichr_res_pw_clusters`),
2. per `(Cluster, Genes)`, collapse the pathway names into one `" | "`-joined
   string (`enrichr_res_pw_short`),
3. `remove_subset_rows()` — drop any row whose gene set is a *subset* of
   another row's gene set, i.e. keep only the maximal gene sets,
4. write `enrich_results_short_<cond>_spec<cutoff>.csv`,
5. pivot to one row per pathway with a count of how many gene-sets/clusters
   it spans, write `enrich_results_per_pathway_<cond>_spec<cutoff>.csv`.

This "collapse pathways whose genes are a subset of another pathway's"
approach was **replaced by** `add_reactome_hierarchy()`, which instead
collapses redundant pathways by their shared parent in the Reactome tree.
So the block is dead, superseded code — safe to delete.

#### The `add_reactome_hierarchy()` trailing filter — deliberate, keep it

The function's last expression is `pw_df_n_grouped2 %>% filter(hierarchy_n
> 1)`, so that filtered frame is what it returns and what
`do_network_enrichment()` writes to `pathways_<cond>_spec<cutoff>_all.csv`.
`hierarchy_n` = `lengths(strsplit(hierarchy_ids, ";"))` = how many levels
deep the pathway sits in the Reactome tree, so `hierarchy_n <= 1` = top-level
/ shallow pathways.

**This is intended** — top-level Reactome pathways (e.g. "Signal
Transduction", "Metabolism") are too generic to be useful in the result and
are dropped on purpose. When migrating, keep the filter but make it an
explicit `return(pw_df_n_grouped2 %>% filter(hierarchy_n > 1))` with a
one-line comment saying *why*, so it isn't mistaken for a stray
console-inspection line again.

### 3.2 Plot styling → `networkplot_theme()`

Everything `visualize_network_pg()` hardcodes, and its config field:

| Hardcoded now | Value | `networkplot_theme()` field |
|---|---|---|
| type → node shape | 14-entry map (Kinase=`dot`, Sensitivity=`square`, Artificial=`star`, `Hidden`=`text`, …), fallback `ellipse` | `shapes` (named vector, merged over the default) + `shape_default` |
| LogFC diverging palette | `c("#0000CC", "#BDBDBD", "#C62828")`, 201 steps | `logfc_palette` |
| LogFC colour limits | clamp to `±max(|p25 of negatives|, |p75 of positives|)` | `logfc_limit_fn` (fn of the LogFC vector; default = that percentile rule) |
| high-degree highlight | border `#4C9900` vs `#454545`, width `2.5` vs `1`, threshold = `highlight_degree` arg | `highlight_border`, `border_default`, `highlight_border_width`, `border_width_default` |
| edge colour | `#000000`, every edge | `edge_colour` |
| edge width | none (all equal) | `edge_width` — see §3.3 |
| node font size | `10` | `node_font_size` |
| legend font size | `12` | `legend_font_size` |
| layout / seed | `"layout_with_fr"`, `123` | `layout`, `seed` |

Defaults reproduce today's output exactly (un-themed call = no visible
change). Overrides are by name and *merge* into the default map, e.g.
`plot_network(..., theme = networkplot_theme(shapes = c(Kinase = "diamond")))`
changes only Kinase.

**R list constructor, not a config file, for v1** (see Rejected). A
`read_networkplot_theme(path)` YAML/JSON reader producing the same list can
be layered on later without touching `plot_network()`.

### 3.3 Edge boldness = link strength — how `edges$weight` is calculated

`weight` comes out of `PCSF_rand_pg()` (`networkGen/R/pcsf-core.R`). The
algorithm runs PCSF `n` times and counts edge recurrence:

1. **Per run `i` (of `n`):**
   - every PPI edge cost is perturbed up by random noise:
     `cost_i = base_cost * (1 + U(0, r))`, `r = 0.1` default — noise only
     ever *raises* cost, by ≤ 10%.
   - the dummy root edges are all given cost `w`.
   - the compiled Steiner solver (`call_sr`) is run once on those costs,
     returning the edge list of that run's optimal subnetwork.
   - dummy-root edges are stripped; the run's real edges are stored.
2. **Accumulate:** an all-nodes × all-nodes integer matrix `adj_matrix` is
   zeroed. For every run, for every edge `(a, b)` in that run, the pair is
   put in a canonical order (so `(a,b)` and `(b,a)` are the same slot) and
   `adj_matrix[a, b] += 1`.
3. **Emit:** `weight = adj_matrix[a, b]` for every pair with a non-zero
   count.

So **`weight` = the number of the `n` randomized runs in which that edge
appeared in the optimal subnetwork** (integer, 1..`n`). An edge that
survived every run has `weight = n`; one that showed up once has
`weight = 1`. It is a **stability / robustness count** — not a biological
interaction strength, and not the PPI edge's `cost`. (The same run-count
idea drives the output node `prize`: `table()` of node appearances across
runs.)

`weight` is already: returned in `$edges`, written to
`edges_<condition>_spec<cutoff>.csv`, read back when a cached network is
reloaded, and used as the edge weight for
`igraph::edge.betweenness.community()` clustering.
`visualize_network_pg()` is the only place that drops it — it hardcodes
`edges_vis$color <- "#000000"` and sets no `width`.

**The feature (entirely inside `networkPlot`):** `theme$edge_width` is
either `NULL` (constant width — today) or
`list(column = "weight", scale = function(x) scales::rescale(x, to = c(1, 6)))`
→ `edges_vis$width <- scale(edges[[column]])`. visNetwork renders `width`
as edge boldness natively. Monotonic with `weight`, so no inversion. Add a
one-line legend note ("line thickness = runs supporting the edge") when it
is active. Fallback to constant width, with a warning, if the requested
column is absent (older CSVs).

#### Preferred metric — biological interaction cost — DONE in `networkGen`

Boldness should show the **PPI edge `cost`** (the reference network's
interaction confidence; lower = stronger), not the run-frequency `weight`.

`cost` is a column of the input PPI, but the PCSF solver returns only
`(from, to)` name pairs and drops it. It is now reattached at the build, in
`PCSF_rand_pg()` (`networkGen/R/pcsf-core.R`):

1. After the consensus `result_edges` (`from`, `to`, `weight`) is built, it
   is keyed **direction-agnostically** — `paste(pmin(from, to), pmax(from,
   to))`, so `(A,B)` and `(B,A)` collapse to one key (the solver's output
   orientation comes from `which(adj_matrix > 0, arr.ind = TRUE)` and is
   arbitrary).
2. That key is matched against the same key built from the **simplified**
   interactome (`edges_simplified` / `edge_weights_simplified`) — the exact
   edge list and costs the solver searched over, so whatever
   `igraph::simplify()` did with any parallel edges is inherited for free;
   no separate parallel-edge rule.
3. `cost` then rides automatically into the returned `$edges`, into
   `edges_<condition>_spec<cutoff>.csv` (additive column — old readers that
   `select()` are unaffected; `networkScore` only touches `weight`), and the
   reloaded-from-disk path. `CONTEXT.md` gains an "Edge cost / edge weight"
   vocab entry.
4. `networkPlot`'s `theme$edge_width` then supports
   `list(column = "cost", scale = ..., invert = TRUE)` — `invert` because
   low cost = bold; default scale e.g.
   `function(x) scales::rescale(-x, to = c(1, 6))`.

Caveats:

- **Only PPI edges have a cost.** Every edge in a PCSF result *is* a PPI
  edge (dummy-root edges are stripped; artificial terminal nodes still
  connect through real PPI edges), so in practice there should be no `NA` —
  but the plot code should fall back to the default width for any `NA` and
  warn, in case of an unmatched pair.
- **Older `edges_*.csv` files** (written before this change) won't have a
  `cost` column; `plot_network()` warns and falls back to constant width
  when the requested column is absent.
- Both metrics coexist — `weight` stays in the data frame alongside `cost`,
  and `theme$edge_width$column` picks which drives boldness.

### 3.4 Pluggable pathway database

Today only **Reactome** is wired end to end. Enrichr itself is already
database-agnostic — `enrich_network(databases = ...)` is passed straight
through to the Enrichr call, and Enrichr ships many built-in libraries
(KEGG, WikiPathways, GO, MSigDB Hallmark, …). What is Reactome-specific is
the **post-processing**: `add_reactome_hierarchy()` collapses redundant
terms by shared Reactome parent and drops top-level terms, using the
`ReactomePathways.txt` / `ReactomePathwaysRelation.txt` tables.

Design so the reference database can be swapped:

- `enrich_network(network, nodes, databases = "Reactome_Pathways_2024",
  postprocess = reactome_postprocess(pathways_txt, relations_txt),
  min_hits = 2, ...)`.
- `postprocess` is a function `(pathway_df) -> pathway_df` — one step that
  takes the raw Enrichr hits and returns the filtered/collapsed table.
  `reactome_postprocess()` is the built-in that does hierarchy collapsing +
  shallow-term drop + WikiPathways-ontology filtering; passing
  `postprocess = NULL` (or a lighter generic filter) is what a KEGG /
  Hallmark run uses until a database-specific step exists for it.
- Reference tables are arguments to `reactome_postprocess()`, never
  package-bundled or global.
- **Next step (not v1):** try other Enrichr libraries by passing
  `databases = c("KEGG_2021_Human", ...)` and building matching
  `*_postprocess()` steps; the interface above is what makes that additive.
- Open question for the design pass: keep the current hand-rolled Enrichr
  HTTP client (`call_enr_simple()`), or move to the `enrichR` CRAN package
  (maintained wrapper, same service). See §5.

### 3.5 Input / output data

**Input** — a `networkGen` result:

- `nodes`: `Protein`, `type`, `prize`, `LogFC` (+ `LogFC_all` for the hover
  label), `degree`; `pathway` added by `enrich_network()`.
- `edges`: `from`, `to`, `weight` (run-frequency count), `cost` (PPI
  interaction cost, lower = stronger).
- `wc_df` (optional): `id`, `cluster`.

**Reference tables** for enrichment (passed in, not bundled): Reactome
`ReactomePathways.txt`, `ReactomePathwaysRelation.txt`.

**Output**: `pathways_<cond>_spec<cutoff>.csv` (+ `_all.csv`),
`nodes_<cond>_spec<cutoff>_with_pathways.csv`,
`<cond>_spec<cutoff>.html`, and the kinase×pathway heatmap PNGs.

## 4. Code reference

Monolith source → destination:

| Monolith (`Network_generation/R/…`) | → | Notes |
|---|---|---|
| `network_enrichment_and_vis.R::visualize_network_pg()` | `networkPlot::plot_network()` | styling → `networkplot_theme()` |
| `kinograte_PG.R` HTML-save tryCatch block | `networkPlot::save_network_html()` | |
| `network_enrichment_and_vis.R::do_network_enrichment()` (simple path) | `networkPlot::enrich_network()` | drop the `per_cluster` arg |
| `network_enrichment_and_vis.R::call_enr_simple()` | `networkPlot` (internal) | |
| `network_enrichment_and_vis.R::filter_wps_by_ontology()` + SPARQL | `networkPlot` (internal) | **rewire** into `enrich_network()` when WikiPathways is requested |
| `network_enrichment_and_vis.R::add_reactome_hierarchy()`, `id_of()`, `get_ancestors()` | `networkPlot` (internal, behind `reactome_postprocess()`) | keep the deliberate `hierarchy_n > 1` drop as an explicit `return()` (§3.1) |
| `network_enrichment_and_vis.R::reconcile_pathway_selection()` | `networkPlot::reconcile_pathways()` | |
| `plotting_functions.R::plot_kinase_pathway_heatmaps()` + `create_combined_heatmap()` / `create_pairwise_combined_heatmaps()` / `create_unique_heatmap_comparison()` (lines ~1009–1655) | `networkPlot::plot_pathway_heatmaps()` + internals | `create_union_combined_heatmap()` is already dead (call site commented) — drop |
| `network_enrichment_and_vis.R`: `network_enrichment_pg`, `enrichment_analysis_pg`, `call_enr_pg`, `reactome_pw_hierarchy_vis`, tail comment block | **delete from monolith** | §3.1 |
| `kinograte_PG.R::augment_by_threshold_steps()` | **delete from monolith** | §3.1 |
| `plotting_functions.R` lines ~10–1008 (MOFA `plot_var_explained` / `plot_factors_pg` / `.set_*`; PCA `plot_pca*`; pathway tilemaps `make_pathway_tilemaps` / `plot_pathway_tilemap` / `plot_heatmap`; `plot_rna_pep_hist`) | **neither** | unrelated leftovers from another project; leave in monolith |

## 5. Open questions — resolve before implementation

1. **`enrich_network()` API shape.** Take loose `nodes`/`edges` data frames
   only, or also a bare `networkGen` result object via an S3 method? (Same
   question for `plot_network()`.)
2. **Enrichr client.** Keep the hand-rolled HTTP client (`call_enr_simple()`,
   POST `addList` + GET `export`), or switch to the `enrichR` CRAN package
   (maintained wrapper on the same service, gives `listEnrichrDbs()` etc.)?
   Also: the current URLs are the legacy `amp.pharm.mssm.edu` host — move to
   `maayanlab.cloud/Enrichr` regardless.
3. **Which nodes are enriched.** The live path sends only the *kinase*-typed
   nodes to Enrichr (Kinase, Kinase-Peptide, Protein-Kinase, …). The user's
   intent is "enrichment on **all** nodes" — confirm that means every node
   in the network regardless of type (including Steiner/Hidden connectors),
   and change the gene-list selection accordingly.
4. **Caching / resume.** `enrich_network()` currently skips the Enrichr call
   and re-reads `pathways_*.csv` if the file already exists. Keep that
   implicit resume, or make it an explicit `refresh = FALSE` arg?
5. **`min_hits` comparison.** Simple path keeps pathways with
   `n_hits > min_n_hits` (strictly greater); per-cluster used `>=`. Pick one
   — `>=` reads more naturally for a "minimum".
6. **Heatmap scope.** `plot_pathway_heatmaps()` today reads a whole result
   *folder* of `pathways_*.csv` / `nodes_*.csv`. In the package, should it
   take an explicit list of `(pathway_df, nodes_df, label)` triples instead,
   so it has no filesystem-layout knowledge?
7. **Legend construction.** Rebuild the `visLegend addNodes` frame from
   `theme` (shapes + colour stops + the edge-width note) — confirm the
   pathway multi-select (`visOptions(selectedBy = "pathway")`) still works
   after the theme refactor.
8. **Bundled vs. supplied reference tables.** `ReactomePathways.txt` /
   `ReactomePathwaysRelation.txt` — always caller-supplied (like `ppi_network`
   in `networkGen`), or shipped as package data with a refresh helper?

## 6. Vignette

One vignette, `networkPlot.Rmd`, covering the two outputs — **not** a deep
explanation of the enrichment algorithm the way `networkGen`'s vignette
explains PCSF. Sections:

1. **Inputs** — a `networkGen` result (`nodes`/`edges`), and the Reactome
   reference tables for enrichment.
2. **Enrichment** — one `enrich_network()` call. State plainly: it runs
   Enrichr on **all nodes of the network** and returns a pathway table plus
   a `pathway` column joined onto `nodes`. **Currently only Reactome is
   supported**, but the reference database is a swap point: `databases` goes
   straight to Enrichr (which has many built-in libraries — KEGG,
   WikiPathways, GO, …) and the Reactome-specific hierarchy collapsing is a
   pluggable `postprocess` step. Note as a next step: experiment with the
   other Enrichr libraries.
3. **Reading the pathway output** — what the pathway table columns mean
   (pathway, overlapping genes, hit count, adjusted p, combined score), and
   the kinase×pathway heatmap.
4. **The interactive HTML** — `plot_network()` → `save_network_html()`.
   Walk the tunable HTML features: node shape per node type, the
   LogFC diverging colour scale, high-degree node highlighting, edge colour
   and **edge boldness = interaction cost / run-frequency**, layout, the
   pathway multi-select and node-search box. Say **where the theme lives**:
   every one of those is a field on the list `networkplot_theme()` returns,
   passed as `theme =`; defaults reproduce the standard look, override by
   name.
5. A tiny end-to-end example on a bundled toy network.

## 7. Verification

- **Plot parity**: render one small real result folder with the monolith's
  `visualize_network_pg()` and with `plot_network(theme =
  networkplot_theme())`; diff `widget$x$nodes` / `widget$x$edges` — must be
  identical.
- **Enrichment parity**: `enrich_network()` on a fixture network vs. the
  monolith's `do_network_enrichment(per_cluster = FALSE)` → identical
  `pathways_*.csv` (allowing for the all-nodes vs kinase-nodes gene-list
  change in open question 3, once decided — record the expected diff).
- **Theme override**: `networkplot_theme(shapes = c(Kinase = "diamond"))`
  changes only Kinase rows' `shape`.
- **Edge width**: `edge_width = list(column = "weight")` → a `weight = n`
  edge renders wider than a `weight = 1` edge; `edge_width = NULL`
  reproduces the constant-width baseline exactly.
- **Biological cost on edges**: after the `PCSF_rand_pg()` join, every row
  of a fixture network's `edges` has a non-`NA` `cost` that matches the
  input PPI's cost for that pair (both directions); with `edge_width =
  list(column = "cost", invert = TRUE)` the lowest-cost edge renders
  boldest. A pre-change `edges_*.csv` (no `cost` column) → warning + constant
  width, no error.
- **WikiPathways rewire**: `enrich_network(databases = c("Reactome_...",
  "WikiPathways_2024_Human"))` calls `filter_wps_by_ontology()` and the
  result has disease/drug pathways removed; Reactome-only does not call it.
- **Dead-code removal**: monolith still runs `plot_networks.Rmd`'s "Run
  kinase network" chunk end-to-end with no reference to a deleted symbol.

## 8. Rejected approaches

- **`networkEnrich` and `networkPlot` as separate packages** (the original
  ADR 0001 split). Enrichment output exists only to drive the plot and its
  heatmaps; a package seam across one linear workflow buys nothing. Merged
  into `networkPlot`.
- **`w` as the link-strength signal.** `w` is a single build-wide scalar
  (terminal-to-root parking cost), not an edge attribute. The per-run edge
  frequency (`edges$weight`) is the real per-edge robustness measure and is
  already available to the plotter.
- **Styling config as a YAML/JSON file (v1).** Adds parser + schema +
  missing-file handling for no v1 gain; callers are R code that can pass an
  R list. A file reader can layer on later.
- **Keeping the per-cluster / topGO enrichment paths "just in case".** They
  duplicate the live simple path, carry a heavy Bioconductor dependency
  (topGO), and have demonstrably never run (a `browser()` sits in one).
  Recoverable from git history if ever needed. `filter_wps_by_ontology()`
  is the one piece worth salvaging.
- **Moving `make_network_and_stats*()` in whole.** They mix build, enrich,
  plot, and disk I/O; `networkPlot` takes the enrich + plot halves, the
  monolith keeps the orchestration as thin package calls.
