#' Reference human protein-protein interaction network (STRING v12.5)
#'
#' The PPI network used as `ppi_network` throughout this project's analyses
#' -- a filtered, high-confidence subset of STRING v12.5 human interactions.
#' Built by the STRING-download pipeline (see `DevOpti/STRING_download`'s
#' `generate_string_ppi.rmd` / `R/helper.R::generate_string_ppi()`), not
#' regenerated here; bundled as-is so `networkGen` ships a working reference
#' network out of the box instead of requiring every caller to independently
#' rebuild one. This is the default `ppi_network` for
#' [generate_paired_network()] and [generate_kinase_network()].
#'
#' Built from STRING's "full" per-evidence-channel links file for
#' *Homo sapiens* (taxon 9606): an edge is kept if it has direct,
#' non-orthology-transferred support from at least one of `experiments`,
#' `database` (curated pathway/complex databases -- KEGG, Reactome,
#' MetaCyc, EBI Complex Portal, GO Complexes), `textmining`, or
#' `coexpression` -- STRING's genomic-context channels (neighborhood,
#' fusion, cooccurrence, homology) are not used as an inclusion criterion.
#' Kept edges are then filtered to STRING's overall `combined_score >= 500`
#' (0.5 on the 0-1 scale, i.e. roughly medium-to-high confidence). `cost` is
#' `max(0.01, 1 - combined_score)`, so lower cost = a stronger interaction.
#'
#' @format A tibble with 1,218,528 rows and 3 columns:
#' \describe{
#'   \item{head}{Gene symbol, interaction partner 1.}
#'   \item{tail}{Gene symbol, interaction partner 2.}
#'   \item{cost}{Edge cost for PCSF -- lower means a stronger/more confident interaction.}
#' }
#' @source STRING v12.5 (<https://string-db.org>), filtered as described
#'   above; see the STRING publication (<https://academic.oup.com/nar/article/53/D1/D730/7903368>)
#'   for what each evidence channel means.
"ppi_network_human_filtered_v12.5"

#' Reference PPI network, kinase-only subset (502 kinases)
#'
#' A restriction of an earlier STRING v12 build of the reference network
#' (see [ppi_network_human_filtered_v12.5]) to interactions involving a
#' curated set of 502 kinases, 107,508 rows. Same `head`/`tail`/`cost`
#' schema. Not the default -- pass explicitly as `ppi_network` when working
#' with kinases only (e.g. faster iteration during development). Not yet
#' regenerated from the current v12.5 default network.
#'
#' @format A tibble with 107,508 rows and 3 columns (`head`, `tail`, `cost`).
#' @source STRING v12, filtered the same way as
#'   [ppi_network_human_filtered_v12.5] (an earlier build).
"ppi_networkv12_502_kins"
