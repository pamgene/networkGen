# networkGen

Builds PCSF (prize-collecting Steiner forest) networks from kinase-activity
data, optionally paired with sensitivity data. Foundation package of a
three-package suite -- see [`networkScore`](../networkScore) and
[`networkPlot`](../networkPlot), and the suite-wide
[Context Map](./CONTEXT-MAP.md).

## Reference PPI network

The default `ppi_network` (`ppi_network_human_filtered_v12.5`, bundled as
package data) is built from **STRING v12.5** (*Homo sapiens*) by the
`DevOpti/STRING_download` pipeline: an edge is kept if it has direct,
non-orthology-transferred evidence from at least one of `experiments`,
`database` (curated pathway/complex databases -- KEGG, Reactome, MetaCyc,
EBI Complex Portal, GO Complexes), `textmining`, or `coexpression`
(STRING's genomic-context channels -- neighborhood, fusion, cooccurrence,
homology -- aren't used as an inclusion rule), then filtered to STRING's
overall confidence `>= 0.5`. `cost = max(0.01, 1 - confidence)`. See
`?ppi_network_human_filtered_v12.5` for details, and `DevOpti/STRING_download`
(`generate_string_ppi.rmd` / `R/helper.R`) for the build script.

## Install

```r
# install.packages("remotes")
remotes::install_github("IOR-Bioinformatics/PCSF")  # compiled dependency, needs a C++ toolchain
remotes::install_local("path/to/networkGen")
```

## Local development notes

`PCSF` is a compiled package -- if it's installed under an R version with no
matching Rtools, `install_github()` silently produces a package with no
compiled binary at all (no `libs/` dir), and every real PCSF call fails with
`"_PCSF_call_sr" not available for .Call()"`, with no error at install time
to flag it. If you hit this, check which R version actually has a working,
compiled PCSF install (`libs/` present) with a matching Rtools version, and
run dev commands (`R CMD build`/`check`, `testthat`,
`roxygen2::roxygenise()`, vignette rendering) via that R version's explicit
path rather than bare `Rscript`/`R` (which resolve to whichever R happens
to be first on `PATH`), e.g.:

```
"C:\Program Files\R\R-4.3.0\bin\Rscript.exe" -e '...'
"C:\Program Files\R\R-4.3.0\bin\R.exe" CMD check ...
```

Before trusting any "it works" result (tests pass, vignette renders,
`R CMD check` is clean), check the actual output for real PCSF activity, not
just that files/folders got created -- folder/`params.csv` creation happens
regardless of whether the underlying PCSF solve actually succeeded. A silent
`NULL` result and a "no subnetwork" folder look identical to a real success
at the file-listing level.

## Getting started

See `vignette("networkGen")` for a full walkthrough on toy data.

```r
library(networkGen)

uka_hits <- uka_top(uka, spec_cutoff = 1.0, perc_cutoff = 0.7)
sens_hits <- sens_top(sens, perc_cutoff = 0.7)

result <- generate_paired_network(
  uka = uka_hits, sens = sens_hits, ppi_network = ppi_network,
  spec_cutoff = 1.0, b = 1.5, condition = "my_condition"
)
```
