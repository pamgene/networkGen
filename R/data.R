#' Reference human protein-protein interaction network, direct evidence (STRING v12.5)
#'
#' The reference PPI network -- `string_12.5_direct` -- is built from STRING
#' v12.5 (human) with direct evidence only: the `experiments`, `database`
#' and `textmining` channels. Coexpression, evidence transferred from other
#' species and genomic context (neighborhood, fusion, co-occurrence) are not
#' used. The edge score is recomputed from those three channels with
#' STRING's own combination rule (STRING FAQ "How are the scores
#' computed?", `combine_subscores.v2.py`), kept if `>= 0.5`, and converted
#' to cost: `cost = 1 - score`. 18,440 proteins, 431,946 interactions,
#' each stored in both directions. This is the default `ppi_network` for
#' [generate_paired_network()], [generate_kinase_network()],
#' [build_network_grid()] and [run_network_grid()]. Built by
#' `data-raw/string_12.5_direct.R`.
#'
#' @format A tibble with 863,892 rows and 3 columns:
#' \describe{
#'   \item{head}{Gene symbol, interaction partner 1.}
#'   \item{tail}{Gene symbol, interaction partner 2.}
#'   \item{cost}{Edge cost for PCSF -- lower means a stronger/more confident interaction.}
#' }
#' @source <https://github.com/pamgene/STRING_download>, <https://string-db.org/help/faq/>.
"string_12.5_direct"

#' Previous reference PPI network (STRING v12.5, before 2026-11-11)
#'
#' The default `ppi_network` until 2026-11-11, kept to reproduce earlier
#' results. An edge is kept if at least one of the `experiments`,
#' `database`, `textmining` or `coexpression` channels has direct evidence
#' (> 0) and STRING's full `combined_score` is `>= 0.5`. That combined
#' score also contains evidence transferred from other species and genomic
#' context, so those still decide which edges pass and set their cost.
#' `cost = 1 - (combined score/1000)`. Superseded by [string_12.5_direct].
#'
#' @format A tibble with 1,218,528 rows and 3 columns:
#' \describe{
#'   \item{head}{Gene symbol, interaction partner 1.}
#'   \item{tail}{Gene symbol, interaction partner 2.}
#'   \item{cost}{Edge cost for PCSF -- lower means a stronger/more confident interaction.}
#' }
#' @source <https://github.com/pamgene/STRING_download>.
"string_12.5"

#' Reference PPI network, kinase-only subset (502 kinases)
#'
#' A restriction of an earlier STRING v12 build of the reference network
#' (see [string_12.5_direct]) to interactions involving a curated set of 502
#' kinases, 107,508 rows. Same `head`/`tail`/`cost` schema. Not the
#' default -- pass explicitly as `ppi_network` when working with kinases
#' only (e.g. faster iteration during development). Not yet regenerated
#' from the current v12.5 default network.
#'
#' @format A tibble with 107,508 rows and 3 columns (`head`, `tail`, `cost`).
#' @source STRING v12, filtered the same way as
#'   [string_12.5] (an earlier build).
"ppi_networkv12_502_kins"
