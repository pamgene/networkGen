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
| `enrich_network(nodes, databases = "Reactome_Pathways_2024", min_hits = 2, postprocess = reactome_postprocess(reactome_refs()), out_dir = NULL, refresh = FALSE)` | one `enrichR` call on every network node except `Hidden` (Steiner connectors) → filtered, post-processed pathway table + a `pathway` column joined onto `nodes`; if `out_dir` given, writes `pathways_*_all.csv` / `pathways_*.csv` and, unless `refresh`, resumes from an existing `pathways_*.csv` | data frames; writes only when `out_dir` given |
| `reconcile_pathways(all_tables)` | given several comparisons' `_all` candidate tables, re-pick one pathway per gene-set favouring pathways recurring across the most comparisons; returns each comparison's reconciled short table | pure |
| `plot_network(nodes, edges, title, clusters = NULL, colour_by = "pathway", theme = networkplot_theme())` | build the visNetwork htmlwidget | pure |
| `save_network_html(widget, path)` | the `saveWidget` → `visSave` fallback dance, one place | write |
| `plot_pathway_heatmaps(comparisons, out_dir, min_kinases = 3)` | core: `comparisons` = named list of `list(pathways = df, nodes = df)`; writes per-comparison + common-to-all + unique/pairwise-overlap `ComplexHeatmap` PNGs | write |
| `plot_pathway_heatmaps_dir(res_dir, spec_cutoff, ...)` | thin wrapper: discovers `pathways_*.csv` / `nodes_*.csv` in a result folder and calls `plot_pathway_heatmaps()` | write |
| `networkplot_theme(...)` | config object holding every styling knob, defaults = today's look | — |
| `reactome_refs(pathways_txt = <bundled>, relations_txt = <bundled>)` | load the Reactome reference tables (defaults to the bundled snapshot) | read |
| `download_reactome_refs(dest)` | fetch the current Reactome release to `dest`; run manually to refresh, never automatic | write + network |

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
  caller actually uses): take the network's nodes as *one* gene list, make
  *one* Enrichr call, keep pathways with enough overlapping genes. (In
  `networkPlot` this becomes: every node except `Hidden`; `>= min_hits`
  overlap — see §3.4.)

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

### 3.4 Enrichment — resolved shape

**Enrichr client.** Use the **`enrichR` CRAN package** (`enrichr()`,
`listEnrichrDbs()`), not the hand-rolled `call_enr_simple()` HTTP client.
Point it at `maayanlab.cloud/Enrichr` (the legacy `amp.pharm.mssm.edu` host
the old code used is deprecated).

**Which nodes.** Enrich **every node in the network except `type ==
"Hidden"`** (the Steiner connectors PCSF added, which carry no input
identity). This is a deliberate change from the current code's
kinase-types-only filter (`Kinase`, `Kinase-Peptide`, `Protein-Kinase`,
…) — kept explicit in the docs and a verification diff, not silent.

**Filter.** Keep a pathway when its overlap is `>= min_hits` (the current
simple path uses strictly `>`; `>=` is what a "minimum" should mean).

**Caching.** `refresh = FALSE` (default): if `out_dir` already holds this
comparison's `pathways_*.csv`, read it and skip the Enrichr call.
`refresh = TRUE` forces a fresh call. No implicit file-exists magic beyond
that one documented arg.

**Pluggable reference database.** Only **Reactome** is wired end to end
today. `enrichR` is already database-agnostic — `databases` is passed
straight through, and Enrichr ships many libraries (KEGG, WikiPathways, GO,
MSigDB Hallmark, …). What is Reactome-specific is the **post-processing**:
hierarchy collapsing by shared parent + top-level-term drop
(`add_reactome_hierarchy()`), using the `ReactomePathways.txt` /
`ReactomePathwaysRelation.txt` tables. So:

- `postprocess` is a function `(pathway_df) -> pathway_df`.
  `reactome_postprocess(refs)` is the built-in (hierarchy collapse +
  shallow-term drop + WikiPathways-ontology filter when WikiPathways is
  among `databases`). `postprocess = NULL` skips it — what a KEGG /
  Hallmark run uses until a database-specific step exists.
- Reference tables come from `reactome_refs()`, which defaults to a
  **bundled snapshot** shipped as package data (works out of the box, fast).
  `enrich_network(postprocess = reactome_postprocess(reactome_refs(my_txt,
  my_rel)))` overrides. `download_reactome_refs(dest)` fetches the current
  release — run manually when Reactome cuts a new version (~quarterly),
  never automatically.
- **Next step (not v1):** exercise other Enrichr libraries via `databases =`
  and matching `*_postprocess()` steps; the interface above is what makes
  that additive.

### 3.5 Input / output data

**Input** — a `networkGen` result:

- `nodes`: `Protein`, `type`, `prize`, `LogFC` (+ `LogFC_all` for the hover
  label), `degree`; `pathway` added by `enrich_network()`.
- `edges`: `from`, `to`, `weight` (run-frequency count), `cost` (PPI
  interaction cost, lower = stronger).
- `wc_df` (optional): `id`, `cluster`.

**Reference tables** for enrichment: `reactome_refs()` — a bundled Reactome
snapshot (`ReactomePathways.txt`, `ReactomePathwaysRelation.txt`) shipped as
package data; overridable, refreshable via `download_reactome_refs()`.

**Two-pass pathway selection** (why `reconcile_pathways()` exists):

1. Per comparison, `enrich_network()` writes `pathways_<cond>_spec<cutoff>_all.csv`
   (every surviving candidate) **and** a *provisional*
   `pathways_<cond>_spec<cutoff>.csv` — one pathway per gene-set, picked by
   the hierarchy tiebreak alone, **no cross-comparison awareness**.
2. Once every comparison for a `spec_cutoff` is enriched,
   `reconcile_pathways()` reads all the `_all.csv` files, counts how many
   comparisons each candidate pathway appears in, and **overwrites** every
   comparison's short `pathways_<cond>_spec<cutoff>.csv` with a pick that
   favours pathways recurring across the most comparisons (tie → hierarchy
   depth → combined score → name). This is what makes the same biology line
   up under the same pathway label across comparisons.
3. `plot_pathway_heatmaps()` then just reads the reconciled short tables —
   it does **no** selection of its own.

**Output**: `pathways_<cond>_spec<cutoff>.csv` (+ `_all.csv`),
`nodes_<cond>_spec<cutoff>_with_pathways.csv`,
`<cond>_spec<cutoff>.html`, and the kinase×pathway heatmap PNGs (per
comparison, common-to-all, and unique/pairwise-overlap).

## 4. Code reference

Monolith source → destination:

| Monolith (`Network_generation/R/…`) | → | Notes |
|---|---|---|
| `network_enrichment_and_vis.R::visualize_network_pg()` | `networkPlot::plot_network()` | styling → `networkplot_theme()` |
| `kinograte_PG.R` HTML-save tryCatch block | `networkPlot::save_network_html()` | |
| `network_enrichment_and_vis.R::do_network_enrichment()` (simple path) | `networkPlot::enrich_network()` | drop `per_cluster`; enrich all-but-`Hidden` nodes; `>= min_hits`; `refresh` arg |
| `network_enrichment_and_vis.R::call_enr_simple()` | **replaced** by the `enrichR` package | |
| `network_enrichment_and_vis.R::filter_wps_by_ontology()` + SPARQL | `networkPlot` (internal, inside `reactome_postprocess()`) | called when `WikiPathways_2024_Human` is among `databases` |
| `network_enrichment_and_vis.R::add_reactome_hierarchy()`, `id_of()`, `get_ancestors()` | `networkPlot` (internal, inside `reactome_postprocess()`) | keep the deliberate `hierarchy_n > 1` drop as an explicit `return()` (§3.1) |
| `network_enrichment_and_vis.R::reconcile_pathway_selection()` | `networkPlot::reconcile_pathways()` | folder-reading split into a `_dir()` wrapper; core takes the `_all` tables |
| `plotting_functions.R::plot_kinase_pathway_heatmaps()` + `create_combined_heatmap()` / `create_pairwise_combined_heatmaps()` / `create_unique_heatmap_comparison()` (lines ~1009–1655) | `networkPlot::plot_pathway_heatmaps()` (core, takes dfs) + `plot_pathway_heatmaps_dir()` (discovery wrapper) + internals | `create_union_combined_heatmap()` is already dead (call site commented) — drop |
| `network_enrichment_and_vis.R`: `network_enrichment_pg`, `enrichment_analysis_pg`, `call_enr_pg`, `reactome_pw_hierarchy_vis`, tail comment block | **delete from monolith** | §3.1 |
| `kinograte_PG.R::augment_by_threshold_steps()` | **delete from monolith** | §3.1 |
| `plotting_functions.R` lines ~10–1008 (MOFA `plot_var_explained` / `plot_factors_pg` / `.set_*`; PCA `plot_pca*`; pathway tilemaps `make_pathway_tilemaps` / `plot_pathway_tilemap` / `plot_heatmap`; `plot_rna_pep_hist`) | **neither** | unrelated leftovers from another project; leave in monolith |

## 5. Resolved design decisions

1. **API shape.** `plot_network()` / `enrich_network()` take loose data
   frames only (`nodes`, `edges`). No S3 method on `networkGen_result` —
   the caller unpacks `res$nodes` / `res$edges` in one line; nothing varies
   across that seam.
2. **Enrichr client.** The `enrichR` CRAN package, `maayanlab.cloud`
   endpoint. `call_enr_simple()` is dropped.
3. **Nodes enriched.** Every node except `type == "Hidden"` (Steiner
   connectors). Deliberate change from the current kinase-types-only
   filter; documented + covered by a verification diff.
4. **Caching.** Explicit `refresh = FALSE` arg; default resumes from an
   existing `pathways_*.csv` in `out_dir`.
5. **Filter.** Keep pathways with overlap `>= min_hits`.
6. **Heatmaps.** `plot_pathway_heatmaps()` core takes a named list of
   `list(pathways = df, nodes = df)`; `plot_pathway_heatmaps_dir()` is the
   folder-discovery wrapper. Same output set (per-comparison, common-to-all,
   unique/pairwise-overlap). The two-pass provisional→`reconcile_pathways()`
   selection (§3.5) is unchanged in behaviour.
7. **Legend.** Rebuilt from `theme` (shapes + colour stops + edge-width
   note); pathway multi-select (`visOptions(selectedBy = "pathway")`)
   preserved — a build-time check, not a design fork.
8. **Reference tables.** Bundled Reactome snapshot as package data
   (`reactome_refs()` default), `reactome_refs(path, path)` override,
   `download_reactome_refs(dest)` manual refresh.

## 6. Vignette

One vignette, `networkPlot.Rmd`, covering the two outputs — **not** a deep
explanation of the enrichment algorithm the way `networkGen`'s vignette
explains PCSF. Sections:

1. **Inputs** — a `networkGen` result (`nodes`/`edges`), and the Reactome
   reference tables for enrichment.
2. **Enrichment** — one `enrich_network()` call. State plainly: it runs
   Enrichr (via the `enrichR` package) on **every node of the network
   except the Steiner connectors** (`type == "Hidden"`) and returns a
   pathway table plus a `pathway` column joined onto `nodes`. **Currently
   only Reactome is supported**, but the reference database is a swap point:
   `databases` goes straight to Enrichr (which has many built-in libraries —
   KEGG, WikiPathways, GO, …) and the Reactome-specific hierarchy collapsing
   is a pluggable `postprocess` step. Note as a next step: experiment with
   the other Enrichr libraries.
3. **Reading the pathway output** — what the pathway table columns mean
   (pathway, overlapping genes, hit count, adjusted p, combined score); the
   two-pass selection (provisional per comparison, then reconciled across
   comparisons so labels line up); and the kinase×pathway heatmap.
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
  monolith's `do_network_enrichment(per_cluster = FALSE)` — the pathway
  table matches **except** where the wider gene list (all-but-`Hidden` vs
  kinase-types-only) legitimately adds pathways; assert the extra rows all
  trace to non-kinase, non-Hidden nodes.
- **Node set**: `enrich_network()`'s gene list = `nodes$Protein[nodes$type
  != "Hidden"]`; a fixture with a `Hidden` node confirms it is excluded and
  every other type is included.
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
- **An S3 method taking a `networkGen_result` directly.** One shape
  (in-memory result) vs another (reconstructed from CSVs) both reduce to
  `(nodes, edges)`; the wrapper would be a seam with nothing varying across
  it. Caller unpacks in one line.
- **Enriching every node including `Hidden` connectors.** Steiner nodes
  PCSF invented to bridge terminals carry no input identity; including them
  adds noise to the gene list. Every *other* type is a real hit worth
  enriching, so the rule is all-but-`Hidden`, not kinase-only and not
  literally-all.
- **Keeping the hand-rolled Enrichr HTTP client.** `enrichR` is a
  maintained CRAN wrapper on the same service with `listEnrichrDbs()` —
  strictly better for the "try other libraries" goal.
