# Split network generation, scoring, enrichment, and plotting into four packages

> **Partially superseded by [0008](./0008-merge-enrich-into-plot.md):**
> `networkEnrich` and `networkPlot` were merged into a single `networkPlot`
> package, so the suite is now three packages, not four. The foundation +
> layered-packages reasoning below still stands.

`Network_generation`'s network-building, scoring, enrichment, and plotting
logic is wanted in other projects and needs to interoperate cleanly, so it's
becoming installable packages instead of `source()`d scripts. We considered
one combined package (generation + scoring + enrichment + plotting) versus
several smaller ones layered on a foundation.

Decided: four packages — `networkGen` (foundation), `networkScore`,
`networkEnrich`, `networkPlot` (each depending on `networkGen`, not on each
other, except `networkPlot` optionally consuming `networkEnrich`'s output
shape without a package dependency) — built in two phases: `networkGen` +
`networkScore` first, `networkEnrich`/`networkPlot` designed against the
same contract now but built later.

Why: `Network_generation` itself only needs generation (not scoring's
permutation/parallel dependencies, nor necessarily enrichment/plotting as
separate installs). Other future consumers may want only one piece. A single
combined package would force every consumer to take all of it. The cost —
four repos to version and keep compatible instead of one — was judged worth
it for that flexibility.
