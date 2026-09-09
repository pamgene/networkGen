#' Reference human protein-protein interaction network (STRING v12)
#'
#' The PPI network used as `ppi_network` throughout this project's analyses
#' -- a filtered, high-confidence subset of STRING v12 human interactions.
#' Built by the STRING-download pipeline (see `DevOpti/STRING_download` /
#' `150-110 Golden data/Dev2024-17-GDS_Dora/generate_string_ppi_for_kinograte.r`),
#' not regenerated here; bundled as-is so `networkGen` ships a working
#' reference network out of the box instead of requiring every caller to
#' independently rebuild one. This is the default `ppi_network` for
#' [generate_paired_network()] and [generate_kinase_network()].
#'
#' @format A tibble with 1,138,434 rows and 3 columns:
#' \describe{
#'   \item{head}{Gene symbol, interaction partner 1.}
#'   \item{tail}{Gene symbol, interaction partner 2.}
#'   \item{cost}{Edge cost for PCSF -- lower means a stronger/more confident interaction.}
#' }
#' @source STRING v12 (<https://string-db.org>), filtered to high-reliability
#'   edges (excludes transferred, textmining-only, and coexpression-only evidence).
"ppi_networkv12"

#' Reference PPI network, filtered subset
#'
#' A further-filtered subset of [ppi_networkv12], 855,626 rows. Same
#' `head`/`tail`/`cost` schema and source. Not the default -- pass explicitly
#' as `ppi_network` when this narrower filtering is wanted instead of the
#' full reference network.
#'
#' @format A tibble with 855,626 rows and 3 columns (`head`, `tail`, `cost`) -- see [ppi_networkv12].
#' @source Same as [ppi_networkv12].
"ppi_networkv12_filt"

#' Reference PPI network, kinase-only subset (502 kinases)
#'
#' [ppi_networkv12] restricted to interactions involving a curated set of 502
#' kinases, 107,508 rows. Same `head`/`tail`/`cost` schema and source. Not
#' the default -- pass explicitly as `ppi_network` when working with kinases
#' only (e.g. faster iteration during development).
#'
#' @format A tibble with 107,508 rows and 3 columns (`head`, `tail`, `cost`) -- see [ppi_networkv12].
#' @source Same as [ppi_networkv12].
"ppi_networkv12_502_kins"
