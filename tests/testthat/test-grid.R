test_that("build_network_grid crosses condition x spec_cutoff x perc_cutoff", {
  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("condA", "condB"), each = 3),
    uniprotname = rep(c("K1", "K2", "K3"), 2),
    LogFC = c(3, 2, 1, 3, 2, 1),
    fscore = 2
  )

  grid <- build_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast",
    spec_cutoff = c(0, 0.5), perc_cutoff = 0, b = 1
  )

  # 2 conditions x 2 spec_cutoff x 1 perc_cutoff = 4 cells
  expect_length(grid, 4)
  combos <- purrr::map(grid, ~ list(condition = .x$condition, spec_cutoff = .x$spec_cutoff))
  expect_setequal(
    purrr::map_chr(combos, ~ paste(.x$condition, .x$spec_cutoff)),
    c("condA 0", "condA 0.5", "condB 0", "condB 0.5")
  )
})

test_that("build_network_grid splits by dataset when a dataset column is present, before cleaning", {
  raw <- data.frame(
    dataset = rep(c("ds1", "ds2"), each = 6),
    Sgroup_contrast = rep(c("condA", "condB"), times = 2, each = 3),
    uniprotname = rep(c("K1", "K2", "K3"), 4),
    LogFC = rep(c(3, 2, 1), 4),
    fscore = 2
  )

  # clean_fn here is just a pass-through that drops the dataset column,
  # mirroring how a real clean_uka_to_kinograte*() select() would drop it
  drop_dataset <- function(x) x[, setdiff(colnames(x), "dataset"), drop = FALSE]

  grid <- build_network_grid(
    raw, clean_fn = drop_dataset, condition_col = "Sgroup_contrast",
    spec_cutoff = 0, perc_cutoff = 0, b = 1
  )

  # 2 datasets x 2 conditions x 1 spec_cutoff x 1 perc_cutoff = 4 cells
  expect_length(grid, 4)
  expect_setequal(purrr::map_chr(grid, "dataset"), c("ds1", "ds2"))
  expect_setequal(purrr::map_chr(grid, "condition"), c("condA", "condB"))
})

test_that("run_network_grid builds one task per grid cell, kinase-only when sens is NULL", {
  future::plan(future::sequential)
  local_mocked_bindings(
    generate_kinase_network = function(uka, condition, spec_cutoff, b, w, ppi_network, write) {
      list(condition = condition, spec_cutoff = spec_cutoff, n_terminals = nrow(uka))
    }
  )

  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("condA", "condB"), each = 3),
    uniprotname = rep(c("K1", "K2", "K3"), 2),
    LogFC = c(3, 2, 1, 3, 2, 1),
    fscore = 2
  )

  out <- run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
    ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), write = FALSE
  )

  expect_length(out, 2)
  expect_equal(sort(purrr::map_chr(out, ~ .x$meta$condition)), c("condA", "condB"))
  expect_equal(out[[1]]$result$condition, out[[1]]$meta$condition)
})

test_that("run_network_grid reuses one sens profile across every condition, per sens_perc_cutoff", {
  future::plan(future::sequential)
  captured_sens <- list()
  local_mocked_bindings(
    generate_paired_network = function(uka, sens, condition, spec_cutoff, b, w, ppi_network, write) {
      captured_sens[[condition]] <<- sens
      list(condition = condition)
    }
  )

  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("dose1", "dose2"), each = 2),
    uniprotname = rep(c("K1", "K2"), 2),
    LogFC = c(3, 2, 3, 2),
    fscore = 2
  )
  raw_sens <- data.frame(uniprotname = paste0("T", 1:5), LogFC = c(-5, -4, -3, -2, -1))

  out <- run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
    sens = raw_sens, sens_perc_cutoff = 0.7,
    ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), write = FALSE
  )

  expect_length(out, 2)
  expect_identical(captured_sens[["dose1"]], captured_sens[["dose2"]])
  expect_setequal(captured_sens[["dose1"]]$name, c("T1", "T2", "T3"))
})

test_that("run_network_grid grids sens_perc_cutoff independently of every uka-side dimension", {
  future::plan(future::sequential)
  captured_sens_names <- list()
  local_mocked_bindings(
    generate_paired_network = function(uka, sens, condition, spec_cutoff, b, w, ppi_network, write) {
      captured_sens_names[[length(captured_sens_names) + 1]] <<- sens$name
      list(condition = condition)
    }
  )

  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("dose1", "dose2"), each = 2),
    uniprotname = rep(c("K1", "K2"), 2),
    LogFC = c(3, 2, 3, 2),
    fscore = 2
  )
  raw_sens <- data.frame(uniprotname = paste0("T", 1:5), LogFC = c(-5, -4, -3, -2, -1))

  out <- run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
    sens = raw_sens, sens_perc_cutoff = c(0.7, 0.9),
    ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), write = FALSE
  )

  # 2 conditions x 2 sens_perc_cutoff values -- 4 tasks, not 2.
  expect_length(out, 4)

  # sens_perc_cutoff = 0.7 (threshold 0.5 after the default balance shift)
  # keeps T1/T2/T3 (3 hits); 0.9 (threshold 0.7) keeps T1/T2 (2 hits) --
  # confirms each value actually reaches sens_top() and changes what's
  # built, not just carried through unused.
  n_hits_by_cutoff <- vapply(seq_along(out), function(i) length(captured_sens_names[[i]]), integer(1))
  cutoffs <- vapply(out, function(x) x$meta$sens_perc_cutoff, numeric(1))
  expect_true(all(n_hits_by_cutoff[cutoffs == 0.7] == 3))
  expect_true(all(n_hits_by_cutoff[cutoffs == 0.9] == 2))
})

test_that("run_network_grid passes a single-network grid's ppi_network once to the batch, not per-task; a multi-network grid keeps it per-task", {
  future::plan(future::sequential)
  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("condA", "condB"), each = 2),
    uniprotname = rep(c("K1", "K2"), 2), LogFC = c(3, 2, 3, 2), fscore = 2
  )
  ppi_one <- data.frame(head = "A", tail = "B", cost = 0.1)
  ppi_two <- data.frame(head = "C", tail = "D", cost = 0.2)

  capture <- NULL
  local_mocked_bindings(
    generate_networks_batch = function(tasks, generate_fn, ppi_network = NULL, extra_args = list(), progress = FALSE) {
      capture <<- list(
        n = length(tasks),
        per_task_ppi = vapply(tasks, function(t) "ppi_network" %in% names(t$args), logical(1)),
        shared = ppi_network
      )
      lapply(tasks, function(t) list(result = list(), meta = t$meta))
    }
  )

  run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
    ppi_network = ppi_one, write = FALSE
  )
  expect_false(any(capture$per_task_ppi))          # not embedded in each task
  expect_identical(capture$shared, ppi_one)        # passed once as the shared network

  run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
    ppi_network = list(one = ppi_one, two = ppi_two), write = FALSE
  )
  expect_true(all(capture$per_task_ppi))           # genuine multi-network grid: kept per-task
  expect_null(capture$shared)
})

test_that("run_network_grid requires sens_perc_cutoff when sens is given", {
  cleaned_uka <- data.frame(
    Sgroup_contrast = "dose1", uniprotname = c("K1", "K2"), LogFC = c(3, 2), fscore = 2
  )
  raw_sens <- data.frame(uniprotname = "T1", LogFC = -5)

  expect_error(
    run_network_grid(
      cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
      sens = raw_sens, ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), write = FALSE
    ),
    "sens_perc_cutoff"
  )
})

test_that("run_network_grid refuses to proceed when the grid exceeds max_tasks", {
  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(paste0("cond", 1:5), each = 2),
    uniprotname = rep(c("K1", "K2"), 5), LogFC = 1, fscore = 2
  )

  expect_error(
    run_network_grid(
      cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 1,
      ppi_network = data.frame(head = "A", tail = "B", cost = 0.1), write = FALSE, max_tasks = 3
    ),
    "max_tasks"
  )
})

test_that("run_network_grid works end-to-end with real PCSF, one folder per (spec_cutoff, perc_cutoff) combination", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  local_path_ppi <- data.frame(
    head = c("A", "B", "D"), tail = c("B", "D", "C"), cost = c(0.1, 0.1, 0.1)
  )
  cleaned_uka <- data.frame(
    Sgroup_contrast = rep(c("cond1", "cond2"), each = 3),
    uniprotname = rep(c("A", "B", "C"), 2),
    LogFC = rep(1, 6),
    fscore = 2
  )

  respath <- file.path(tempdir(), "run_network_grid_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  out <- run_network_grid(
    cleaned_uka, condition_col = "Sgroup_contrast", spec_cutoff = 0, perc_cutoff = 0, b = 10,
    ppi_network = local_path_ppi, respath = respath, write = TRUE
  )

  expect_length(out, 2)
  for (item in out) {
    expect_s3_class(item$result, "networkGen_result")
    expect_true("D" %in% item$result$nodes$Protein)
  }

  files <- list.files(respath, recursive = TRUE)
  expect_true(any(grepl("nodes_cond1", files)))
  expect_true(any(grepl("nodes_cond2", files)))
  expect_true(any(grepl("params\\.csv$", files)))
  # both conditions share ONE folder (same spec_cutoff/perc_cutoff), not two
  expect_equal(length(unique(dirname(files))), 1)
})
