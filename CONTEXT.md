# networkGen

Builds a PCSF (Prize-Collecting Steiner Forest) network from kinase-activity
data, optionally paired with sensitivity data. The foundation package of the
network suite — no scoring, enrichment, or plotting logic lives here. See
the suite-wide [Context Map](./CONTEXT-MAP.md) for shared vocabulary
(network-result object, condition, PPI network) and how this package relates
to the others.

## Language

**Terminal node**:
A node with a known prize — supplied directly from the input data (a
kinase, or a kinase+sensitivity hit), as opposed to a node PCSF adds to
connect terminals together.
_Avoid_: hit, seed node

**Steiner node**:
A node PCSF adds to the network to connect terminal nodes, not itself a
direct hit in the input data.
_Avoid_: hidden node, connector

**Prize**:
The numeric weight PCSF assigns a terminal node, driving how strongly the
algorithm wants to include it in the output network.
_Avoid_: weight (ambiguous with edge weight/cost)

**Edge cost / edge weight** (the two numeric columns on a result's `edges`):
`cost` is the reference PPI's interaction cost for that edge (lower = a
stronger/more confident interaction), carried straight over from the input
PPI network. `weight` is the number of the `n` randomized PCSF runs the
edge survived in (1..`n`) — a robustness count, not a biological quantity;
it's what the community-detection step weights by. Keep the two names
distinct; never call either one just "weight" unqualified.

**Percentile score**:
A node's percentile rank (0–1) on its input metric (e.g. absolute LogFC),
used via a cutoff to decide which nodes become terminals. Renamed from the
legacy `score` to avoid collision with `networkScore`'s significance score,
an unrelated concept that happened to share the word.
_Avoid_: score (reserved for networkScore's meaning), rank

**Fast path**:
The dataframe-based PCSF build (`kinograte_pg_pcsf()`/`PCSF_rand_pg()`),
skipping `PCSF::construct_interactome()`. A slower igraph-based
implementation (formerly `kinograte_pg()`/`PCSF_rand()`, the "original path")
existed alongside it and was removed after confirming bit-exact equivalence
-- see `archive/pcsf-equivalence/` and
`docs/adr/0005-remove-original-pcsf-path.md`.
_Avoid_: "new"/"old" PCSF (ambiguous over time)

**Batch build**:
A single call that builds many independent networks in parallel from a flat
list of build inputs. Has no concept of what the inputs represent (real
conditions, permuted inputs, etc.) — that's the caller's concern.
_Avoid_: parallel run, cluster job
