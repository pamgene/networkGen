test_that("uka_top filters by spec_cutoff then percentile ranks and renames prize", {
  uka <- data.frame(
    uniprotname = c("K1", "K2", "K3"),
    LogFC = c(3.0, -1.0, 0.1),
    fscore = c(2.0, 2.0, 0.5)
  )
  res <- uka_top(uka, spec_cutoff = 1.0, perc_cutoff = 0, rank_uka_abs = TRUE)
  # K3 dropped by spec_cutoff (fscore 0.5 < 1.0)
  expect_equal(sort(res$name), c("K1", "K2"))
  expect_true(all(c("name", "prize", "type", "LogFC") %in% colnames(res)))
  expect_equal(unique(res$type), "Kinase")
})

test_that("sens_top lowers the effective cutoff when balance = TRUE", {
  sens <- data.frame(uniprotname = c("S1", "S2", "S3"), LogFC = c(-3.0, -1.0, 2.0))
  res_balanced <- sens_top(sens, perc_cutoff = 0.8, balance = TRUE)
  res_unbalanced <- sens_top(sens, perc_cutoff = 0.8, balance = FALSE)
  expect_true(nrow(res_balanced) >= nrow(res_unbalanced))
  expect_equal(unique(res_balanced$type), "Sensitivity")
})

test_that("overlap_uka_sens computes relative overlap", {
  expect_equal(overlap_uka_sens(c("A", "B", "C"), c("B", "C", "D")), 2 / 3)
  expect_equal(overlap_uka_sens(character(0), c("A")), 0)
  expect_equal(overlap_uka_sens(c("A"), character(0)), 0) # max(length(sens), 1) guards divide-by-zero
})

test_that("clean_tercen_columns keeps only the last dot-segment of each column name", {
  df <- data.frame(`path.to.Sample` = 1, `other.path.Kinase Name` = 2, check.names = FALSE)
  cleaned <- clean_tercen_columns(df)
  expect_equal(colnames(cleaned), c("Sample", "Kinase Name"))
})

test_that("clean_uka_to_kinograte selects and renames the mean/median columns", {
  uka <- data.frame(
    `x.Sgroup_contrast` = c("A_vs_B"),
    `x.Kinase Name` = c("KIN1"),
    `x.Median Kinase Statistic` = c(1.5),
    `x.Mean Specificity Score` = c(2.0),
    check.names = FALSE
  )
  res <- clean_uka_to_kinograte(uka, cs = FALSE)
  expect_equal(colnames(res), c("Sgroup_contrast", "uniprotname", "LogFC", "fscore"))
  expect_equal(res$uniprotname, "KIN1")
  expect_equal(res$LogFC, 1.5)
  expect_equal(res$fscore, 2.0)
})
