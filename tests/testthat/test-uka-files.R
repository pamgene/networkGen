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
