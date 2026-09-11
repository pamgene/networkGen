# networkGen

Builds PCSF (prize-collecting Steiner forest) networks from kinase-activity
data, optionally paired with sensitivity data. Foundation package of a
three-package suite -- see [`networkScore`](../networkScore) and
[`networkPlot`](../networkPlot), and the suite-wide
[Context Map](./CONTEXT-MAP.md).

## Reference PPI network

The reference PPI network -- `ppi_network_human_filtered_v12.5`, the
default `ppi_network`, bundled as package data -- is built from STRING
v12.5. An edge is kept if it has direct, non-transferred support from at
least one of the `experiments`, `database`, `textmining`, or
`coexpression` evidence channels. The STRING combined score is scaled
between 0-1, then filtered for combined score `>= 0.5`. Combined score is
converted to cost: `cost = 1 - (combined score/1000)`. See:
<https://github.com/pamgene/STRING_download>.

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
