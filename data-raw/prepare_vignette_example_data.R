# Prepares the small, real example datasets used by vignettes/networkGen.Rmd
# (vignettes/example_data/). Not run as part of the package build -- a
# provenance record of how those files were produced, to re-run by hand if
# the source data changes.
#
# Source files are real UKA/sensitivity exports living in the sibling
# Network_generation project, not part of this repo:
#   - single-condition example:  Network_generation/uka/GDSi/csUKA_NUDUL_CDKi_PHA-793887.csv
#   - kinase-only batch example: Network_generation/uka/GDS/UKA_cs_GDS_STK_ivPNKL_wKL1_log10.csv
#   - paired batch example:      same CDKi file above (all 3 doses, not just one)
#     + Network_generation/data/sensitivity_gds_median_stk_cells_vs_RL_lessna.csv
#   - multi-dataset (raw, for load_uka_dataset_files()) example:
#     Network_generation/uka/<real study folder>/*.csv (3 real timepoints,
#     PTK+STK pairs each -- copied verbatim, not cleaned, since
#     load_uka_dataset_files() needs the real
#     UKA_<PTK|STK>_<number>_<dataset name>.csv file names). Adjust
#     `multidataset_src` below to match the actual source folder name.
#
# Adjust `network_generation_repo` below if that project lives elsewhere.
#
# Anonymization (this repo is public): the underlying values are real
# experimental results, and the multi-dataset example's real source folder
# name was a real internal study/collaborator identifier -- dropped on copy
# (files land directly under "multidataset", not a named subfolder; the
# real folder name is deliberately not recorded here either). Every genuine
# score/statistic column across all example files
# additionally gets independent random noise added (add_privacy_noise()
# below) before writing -- large enough that the real values aren't
# recoverable, small enough (relative to each column's own spread) that the
# vignette's qualitative behavior (which rows clear spec_cutoff/perc_cutoff,
# roughly how big the resulting networks are) is unaffected. Identifier/
# structural columns (kinase names, UniProt/Entrez IDs, rank, peptide set
# size) are left untouched -- they're not experimental results, and Entrez/
# UniProt IDs are public gene identifiers, not perturbing them keeps the
# example biologically coherent.

library(networkGen)

set.seed(20260910) # reproducible noise -- re-running this script regenerates identical files

#' Add independent random noise to a data frame's score columns
#'
#' @param df Data frame to perturb.
#' @param cols Names of numeric columns to add noise to (others untouched).
#' @param sd_frac Noise SD as a fraction of each column's own SD -- large
#'   enough that individual real values aren't recoverable, small enough
#'   that the column's overall range/shape stays representative.
add_privacy_noise <- function(df, cols, sd_frac = 0.3) {
  for (col in cols) {
    if (!col %in% colnames(df)) next
    x <- df[[col]]
    noise_sd <- sd_frac * stats::sd(x, na.rm = TRUE)
    df[[col]] <- x + stats::rnorm(length(x), mean = 0, sd = noise_sd)
  }
  df
}

network_generation_repo <- Sys.getenv("NETWORK_GENERATION_REPO", "../Network_generation")
out_dir <- file.path("vignettes", "example_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Single-condition example ------------------------------------------
# One drug (a CDK inhibitor, PHA-793887), filtered down to one dose/timepoint
# (4 uM vs DMSO control) so the vignette can go straight from reading the
# file to uka_top() -- cleaning and condition selection both happen here,
# not in the vignette.
raw_single <- read.csv(
  file.path(network_generation_repo, "uka/GDSi/csUKA_NUDUL_CDKi_PHA-793887.csv"),
  check.names = FALSE
)
cleaned_single <- clean_uka_to_kinograte(raw_single, cs = TRUE)
cleaned_single <- add_privacy_noise(cleaned_single, c("LogFC", "fscore"))
single_condition <- cleaned_single[cleaned_single$Sgroup_contrast == "4_T01 vs DMSO", ]

write.csv(
  single_condition,
  file.path(out_dir, "cdk_inhibitor_PHA-793887_4uM.csv"),
  row.names = FALSE
)
cat("single-condition example:", nrow(single_condition), "rows, condition:",
  unique(single_condition$Sgroup_contrast), "\n")

# --- Paired batch example (same drug, all 3 doses + matching sensitivity) --
# All 3 dose/timepoint conditions from the same CDKi-NUDUL export (not
# filtered to one, unlike the single-condition example above), for a paired
# batch demo: one sensitivity profile (below, cell line NUDUL1) pairs with
# several kinase-activity conditions from that same cell line.
all_doses <- cleaned_single
write.csv(
  all_doses,
  file.path(out_dir, "cdk_inhibitor_PHA-793887_all_doses.csv"),
  row.names = FALSE
)
cat("paired-batch UKA example:", nrow(all_doses), "rows,",
  length(unique(all_doses$Sgroup_contrast)), "conditions:",
  paste(unique(all_doses$Sgroup_contrast), collapse = ", "), "\n")

raw_sens <- read.csv(
  file.path(network_generation_repo, "data/sensitivity_gds_median_stk_cells_vs_RL_lessna.csv"),
  check.names = FALSE
)
cleaned_sens <- clean_sens_to_kinograte(raw_sens, control = NULL, zscore = FALSE)
cleaned_sens <- add_privacy_noise(cleaned_sens, "LogFC")
write.csv(
  cleaned_sens,
  file.path(out_dir, "sensitivity_example.csv"),
  row.names = FALSE
)
cat("sensitivity example:", nrow(cleaned_sens), "rows,",
  length(unique(cleaned_sens$cell_line)), "cell lines\n")

# --- Kinase-only batch (multi-condition) example -------------------------
# A larger kinase-activity dataset (66 conditions), trimmed to 3 conditions
# so the batch example stays fast, then cleaned.
raw_full <- read.csv(
  file.path(network_generation_repo, "uka/GDS/UKA_cs_GDS_STK_ivPNKL_wKL1_log10.csv"),
  check.names = FALSE
)
keep_conditions <- c("_A4-Fuk vs A3-KAW", "_HT vs A4-Fuk", "_MC-116 vs A4-Fuk")
raw_trimmed <- raw_full[raw_full$Sgroup_contrast %in% keep_conditions, ]
cleaned_batch <- clean_uka_to_kinograte(raw_trimmed, cs = TRUE)
cleaned_batch <- add_privacy_noise(cleaned_batch, c("LogFC", "fscore"))

write.csv(
  cleaned_batch,
  file.path(out_dir, "batch_example_3conditions.csv"),
  row.names = FALSE
)
cat("batch example:", nrow(cleaned_batch), "rows,",
  length(unique(cleaned_batch$Sgroup_contrast)), "conditions:",
  paste(unique(cleaned_batch$Sgroup_contrast), collapse = ", "), "\n")

# --- Multi-dataset raw example (for load_uka_dataset_files()) ------------
# 3 real timepoints (01/02/03), each a real PTK+STK pair -- kept in the raw
# Tercen wide-export shape (not cleaned) so the vignette can demonstrate
# load_uka_dataset_files() parsing the real UKA_<PTK|STK>_<number>_<dataset
# name>.csv naming convention and merging each pair, not just reading an
# already-cleaned file. Landed directly under "multidataset" -- not a named
# subfolder -- since load_uka_dataset_files() uses the immediate parent
# folder name as part of the dataset identifier (previously the real
# internal study/collaborator name; see the anonymization note above).
# Spaces in the real file names are replaced with underscores on copy --
# R CMD check flags space-containing file names as non-portable, and
# load_uka_dataset_files()'s parsing is agnostic to the difference.
# Not hardcoded here -- the real folder name is an internal study/
# collaborator identifier and this script's own source is public.
multidataset_src <- file.path(network_generation_repo, "uka", Sys.getenv("MULTIDATASET_SRC_FOLDER"))
multidataset_out <- file.path(out_dir, "multidataset")
dir.create(multidataset_out, showWarnings = FALSE, recursive = TRUE)
src_files <- list.files(multidataset_src, full.names = TRUE)
dest_files <- file.path(multidataset_out, gsub(" ", "_", basename(src_files)))

score_cols <- c("Median Final score", "Mean Significance Score", "Mean Specificity Score", "Median Kinase Statistic")
for (i in seq_along(src_files)) {
  df <- read.csv(src_files[i], check.names = FALSE)
  df <- add_privacy_noise(df, score_cols)
  write.csv(df, dest_files[i], row.names = FALSE)
}
cat("multi-dataset raw example:", length(dest_files), "files written\n")
