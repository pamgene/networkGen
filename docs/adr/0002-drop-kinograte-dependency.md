# Drop the external `kinograte` package dependency

`Network_generation`'s `kinograte_PG_core.R` did `library(kinograte)` and
called its `percentile_rank()`, depending on an external GitHub package
(`CogDisResLab/Kinograte`, reached historically via `kalganem/Kinograte`
before a repo-owner transfer). Investigation found it dormant since April
2022, single-author, and `percentile_rank()` itself is ~15 lines with no
meaningful logic — the package as a whole already covers roughly the same
four-way split this suite is adopting (generation/enrichment/plotting/
scoring), at much smaller scale, but nothing else in it is a candidate for
reuse given how far `Network_generation`'s own code has since diverged and
grown past it.

Decided: copy `percentile_rank()` directly into `networkGen` (renamed
`percentile_score()` — see ADR 0004) rather than taking `kinograte` as a
`networkGen` dependency. `library(kinograte)` is removed from every migrated
file.

Why: a live dependency on a single-author, four-year-dormant package isn't
worth it for 15 lines with no real logic. This also removes the last
external reason `networkGen` would need `PCSF` (`Imports: PCSF`) transitively
through `kinograte` rather than directly and deliberately.
