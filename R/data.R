#' Reference human protein-protein interaction network (STRING v12.5)
#'
#' The reference PPI network -- `ppi_network_human_filtered_v12.5` -- is
#' built from STRING v12.5. An edge is kept if it has direct, non-transferred
#' support from at least one of the `experiments`, `database`, `textmining`,
#' or `coexpression` evidence channels. The STRING combined score is scaled
#' between 0-1, then filtered for combined score `>= 0.5`. Combined score is
#' converted to cost: `cost = 1 - (combined score/1000)`. This is the
#' default `ppi_network` for [generate_paired_network()] and
#' [generate_kinase_network()].
#'
#' @format A tibble with 1,218,528 rows and 3 columns:
#' \describe{
#'   \item{head}{Gene symbol, interaction partner 1.}
#'   \item{tail}{Gene symbol, interaction partner 2.}
#'   \item{cost}{Edge cost for PCSF -- lower means a stronger/more confident interaction.}
#' }
#' @source <https://github.com/pamgene/STRING_download>.
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
