# Parallelize scoring by fully flattening the task list, not by condition alone

*Filed here rather than in `networkScore`'s own repo: suite-wide ADRs are
anchored in `networkGen` (see `../../CONTEXT-MAP.md`), even where, as here,
the decision is substantively about `networkScore`'s design.*

`networkScore` needs many networks built per run: one observed + ~50
permuted per condition, across potentially many conditions. A working
precedent exists (`Dev2024-17-GDS_Dora/R/helper_only_kinase_paralell.R`'s
`make_golden_score(parallel_conditions = TRUE)`): parallelize across
conditions via `foreach %dopar%`, with each condition's 50 permutations run
sequentially inside that worker.

Considered: adopt that condition-level pattern as-is (simpler, proven) vs.
fully flatten every build (all conditions' observed + permutation tasks)
into one list and parallel-map the whole thing in a single call.

Decided: full flattening. Condition-level parallelism only keeps every
worker busy when `n_conditions >= n_workers`; with fewer conditions than
cores it leaves cores idle despite there being plenty of independent
permutation work available. Full flattening gets the same throughput when
conditions are plentiful and additionally uses all cores when they aren't,
through the same `networkGen::generate_networks_batch()` primitive — no
nested worker pools.

Consequences: `networkScore` must build the flat task list itself (tagging
each item with condition, role, and permutation index) and regroup results
by condition afterward — bookkeeping the condition-level precedent got for
free by construction. Checkpointing (resume a partially-completed run)
becomes per-build-task instead of per-condition.

Also established while building this: PCSF's compiled
`.Call("_PCSF_call_sr", ...)` genuinely requires `library(PCSF)` to be
loaded in each worker process (a real, previously-hit failure on Windows
PSOCK clusters, documented in that same ancestor project's
`PCSF_PARALLEL_SOLUTION.md`) — not an instance of the global-environment
coupling this migration otherwise removes. `future`/`furrr`'s default
package re-attachment handles it, verified by an explicit test rather than
assumed.
