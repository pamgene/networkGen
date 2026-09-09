# Rename the per-node `score` column to `percentile_score`

`Network_generation`'s existing per-node `score` column (percentile rank
0–1, inherited from `kinograte::percentile_rank()`, used via
`spec_cutoff`/`perc_cutoff` to pick terminal nodes) and the new
`networkScore` package's significance score are unrelated concepts that
happen to share a word — one is a per-node input-ranking value, the other a
permutation-test p-value-like output for a whole network. With `networkScore`
now a real package name, the collision would be a live source of confusion
for anyone calling both.

Decided: rename the per-node value to `percentile_score` throughout
`networkGen`'s API, freeing "score" exclusively for `networkScore`'s
meaning.

Why now: `networkGen`'s public API is being defined from scratch in this
migration anyway (and `percentile_rank()` is being copied in fresh from the
now-dropped `kinograte` dependency — see ADR 0002), so the rename is free
here. Doing it later, after both packages have shipped and callers exist for
each meaning of "score," would be a breaking change with real migration
cost.
