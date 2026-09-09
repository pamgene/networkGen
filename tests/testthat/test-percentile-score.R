test_that("percentile_score ranks by absolute value, highest first", {
  df <- data.frame(uniprotname = c("P1", "P2", "P3"), LogFC = c(1.0, -3.0, 0.5))
  res <- percentile_score(df, uniprotname, LogFC)

  expect_true(all(c("name", "percentile_score") %in% colnames(res)))
  expect_true(all(res$percentile_score >= 0 & res$percentile_score <= 1))
  # P2 has the largest |LogFC| (3.0) -> highest percentile_score
  expect_equal(res$name[which.max(res$percentile_score)], "P2")
})

test_that("percentile_score deduplicates by name and drops NA metric", {
  df <- data.frame(uniprotname = c("P1", "P1", "P2"), LogFC = c(1.0, 2.0, NA))
  res <- percentile_score(df, uniprotname, LogFC)
  expect_equal(nrow(res), 1)
  expect_equal(res$name, "P1")
})

test_that("percentile_score_noabs ranks by signed value, sign matters", {
  df <- data.frame(uniprotname = c("P1", "P2"), LogFC = c(1.0, -2.0))

  res_hi <- percentile_score_noabs(df, uniprotname, LogFC, rank_lowest_highest = FALSE)
  expect_true(res_hi$percentile_score[res_hi$name == "P1"] > res_hi$percentile_score[res_hi$name == "P2"])

  res_lo <- percentile_score_noabs(df, uniprotname, LogFC, rank_lowest_highest = TRUE)
  expect_true(res_lo$percentile_score[res_lo$name == "P2"] > res_lo$percentile_score[res_lo$name == "P1"])
})

test_that("percentile_score_fast matches percentile_score_noabs on abs ranking direction", {
  df <- data.frame(uniprotname = c("P1", "P2", "P3"), LogFC = c(1.0, -3.0, 0.5))
  res <- percentile_score_fast(df, uniprotname, LogFC)
  expect_true("percentile_score" %in% colnames(res))
  expect_equal(res$name[which.max(res$percentile_score)], "P2")
})

test_that("percentile_score_fast handles edge cases: empty, single row, tied values", {
  empty_df <- data.frame(uniprotname = character(), LogFC = numeric())
  expect_equal(nrow(percentile_score_fast(empty_df, uniprotname, LogFC)), 0)

  single_df <- data.frame(uniprotname = "P1", LogFC = 1.0)
  res_single <- percentile_score_fast(single_df, uniprotname, LogFC)
  expect_equal(nrow(res_single), 1)

  tied_df <- data.frame(uniprotname = c("P1", "P2", "P3"), LogFC = c(1.0, 1.0, 1.0))
  res_tied <- percentile_score_fast(tied_df, uniprotname, LogFC)
  expect_true(all(res_tied$percentile_score == res_tied$percentile_score[1]))
})
