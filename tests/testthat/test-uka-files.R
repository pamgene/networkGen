test_that("load_uka_dataset_files merges PTK/STK pairs sharing a folder + number_dataset suffix", {
  dir1 <- file.path(tempdir(), "uka_files_test", "folderA")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)

  write.csv(data.frame(x = 1:2), file.path(dir1, "UKA_PTK_01_TvsC.csv"), row.names = FALSE)
  write.csv(data.frame(x = 3:4), file.path(dir1, "UKA_STK_01_TvsC.csv"), row.names = FALSE)

  out <- load_uka_dataset_files(list.files(dir1, full.names = TRUE))

  expect_equal(nrow(out), 4)
  expect_equal(length(unique(out$dataset)), 1)
  expect_setequal(out$kinase_type, c("PTK", "STK"))
})

test_that("load_uka_dataset_files treats the same number_dataset suffix in different folders as different datasets", {
  base <- file.path(tempdir(), "uka_files_test2")
  dirA <- file.path(base, "folderA")
  dirB <- file.path(base, "folderB")
  dir.create(dirA, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirB, recursive = TRUE, showWarnings = FALSE)

  write.csv(data.frame(x = 1), file.path(dirA, "UKA_PTK_01_TvsC.csv"), row.names = FALSE)
  write.csv(data.frame(x = 2), file.path(dirB, "UKA_PTK_01_TvsC.csv"), row.names = FALSE)

  out <- load_uka_dataset_files(c(
    file.path(dirA, "UKA_PTK_01_TvsC.csv"),
    file.path(dirB, "UKA_PTK_01_TvsC.csv")
  ))

  expect_equal(length(unique(out$dataset)), 2)
})

test_that("load_uka_dataset_files falls back to the file name for a non-matching path, no kinase_type", {
  dir1 <- file.path(tempdir(), "uka_files_test3")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir1, "UKA_NUDUL_AKTi_MK-2206.csv")
  write.csv(data.frame(x = 1), path, row.names = FALSE)

  out <- load_uka_dataset_files(path)

  expect_equal(out$dataset, "UKA_NUDUL_AKTi_MK-2206")
  expect_true(is.na(out$kinase_type))
})

test_that("load_uka_dataset_files handles a dataset name containing underscores", {
  dir1 <- file.path(tempdir(), "uka_files_test4")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir1, "UKA_PTK_03_24h_Pr2_TvsC.csv")
  write.csv(data.frame(x = 1), path, row.names = FALSE)

  out <- load_uka_dataset_files(path)

  expect_true(grepl("03_24h_Pr2_TvsC$", out$dataset))
  expect_equal(out$kinase_type, "PTK")
})

# Regression test: real csUKA exports carry a Tercen app-instance number in
# their column prefixes that differs file to file (e.g.
# "csUKA_app0.UKA0.Sgroup_contrast" vs "csUKA_app1.UKA00.Sgroup_contrast").
# Binding the raw, still-prefixed frames together (the old behavior) unions
# both prefix variants as separate columns -- real for the file it came
# from, all-NA for every other file's rows. A downstream clean_fn's
# clean_tercen_columns() then collapses both variants onto the same short
# name, and select(all_of(name)) could silently grab the all-NA one instead
# of the real one for a given dataset -- producing an NA condition column
# for that dataset even though the raw file had a real value. Cleaning each
# file's columns before binding (this test) means every file already
# shares the one canonical name, so no duplicate/all-NA column is ever
# created.
test_that("load_uka_dataset_files cleans each file's Tercen prefix before binding, so differing app-instance numbers don't create duplicate/all-NA columns", {
  dir1 <- file.path(tempdir(), "uka_files_test5")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)

  path_a <- file.path(dir1, "UKA_PTK_01_DatasetA.csv")
  path_b <- file.path(dir1, "UKA_PTK_02_DatasetB.csv")

  write.csv(
    data.frame(
      `csUKA_app0.UKA0.Sgroup_contrast` = "DrugA_vs_DMSO",
      `csUKA_app0.UKA0.Kinase Name` = "AKT1",
      check.names = FALSE
    ),
    path_a, row.names = FALSE
  )
  write.csv(
    data.frame(
      `csUKA_app1.UKA00.Sgroup_contrast` = "DrugB_vs_DMSO",
      `csUKA_app1.UKA00.Kinase Name` = "MTOR",
      check.names = FALSE
    ),
    path_b, row.names = FALSE
  )

  out <- load_uka_dataset_files(c(path_a, path_b))

  # exactly one Sgroup_contrast column -- not two prefix variants collapsed
  # onto the same name
  expect_equal(sum(colnames(out) == "Sgroup_contrast"), 1)
  expect_equal(sum(colnames(out) == "Kinase Name"), 1)
  # and it's populated for every row, not NA for the "other" file's rows
  expect_false(anyNA(out$Sgroup_contrast))
  expect_setequal(out$Sgroup_contrast, c("DrugA_vs_DMSO", "DrugB_vs_DMSO"))
})
