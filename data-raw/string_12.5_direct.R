# Builds data/string_12.5_direct.rda: STRING v12.5 (human), direct evidence only -- experiments, database and text
# mining; no coexpression, no evidence transferred from other species, no genomic context. The edge score is
# recomputed from those channels with STRING's own combination rule (STRING FAQ "How are the scores computed?",
# combine_subscores.v2.py), via generate_string_ppi() from https://github.com/pamgene/STRING_download (R/helper.R,
# pinned commit below). Kept: score >= 0.5; cost = 1 - score.
#
# Run from the package root: Rscript data-raw/string_12.5_direct.R
# Inputs are downloaded from STRING (links.full ~117 MB, protein.info ~2 MB) unless local copies are given with the
# environment variables STRING_LINKS_FULL and STRING_INFO.
suppressMessages({ library(dplyr); library(readr) })

STRING_VERSION  <- "12.5"
SPECIES         <- "9606"
HELPER_COMMIT   <- "8aa9b342125a0815f59d98dcbcdbbfcb2a31ed2d"
HELPER_URL      <- sprintf("https://raw.githubusercontent.com/pamgene/STRING_download/%s/R/helper.R", HELPER_COMMIT)
LINKS_FULL_URL  <- sprintf("https://stringdb-downloads.org/download/protein.links.full.v%s/%s.protein.links.full.v%s.txt.gz",
                           STRING_VERSION, SPECIES, STRING_VERSION)
INFO_URL        <- sprintf("https://stringdb-downloads.org/download/protein.info.v%s/%s.protein.info.v%s.txt.gz",
                           STRING_VERSION, SPECIES, STRING_VERSION)

# local copy if given, else download into a temporary folder
input_file <- function(env_var, url) {
  local <- Sys.getenv(env_var)
  if (nzchar(local)) { stopifnot(file.exists(local)); return(local) }
  dest <- file.path(tempdir(), basename(url))
  if (!file.exists(dest)) { op <- options(timeout = 3600); on.exit(options(op)); utils::download.file(url, dest, mode = "wb") }
  dest
}

source(HELPER_URL)   # generate_string_ppi(), combine_string_channels()
links_full <- read_delim(input_file("STRING_LINKS_FULL", LINKS_FULL_URL), delim = " ", show_col_types = FALSE,
                         progress = FALSE, col_types = cols(protein1 = "c", protein2 = "c", .default = "i"))
info <- read_tsv(input_file("STRING_INFO", INFO_URL), show_col_types = FALSE, progress = FALSE)

string_12.5_direct <- generate_string_ppi(NULL, info, links_full, best_evidence = c("experiments", "database", "textmining"),
                                          transferred = FALSE, genomic = FALSE, min_score = 0.5)
stopifnot(!anyNA(string_12.5_direct$head), !anyNA(string_12.5_direct$tail), nrow(string_12.5_direct) == 863892)
save(string_12.5_direct, file = "data/string_12.5_direct.rda", compress = "xz")
cat(nrow(string_12.5_direct), "rows,", length(unique(c(string_12.5_direct$head, string_12.5_direct$tail))), "proteins\n")
