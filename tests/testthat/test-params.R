test_that("capture_params generically captures arbitrary parameters, not just UKA-shaped ones", {
  params <- capture_params(x = 1, y = "hello", some_var <- 42)

  expect_equal(params$x, 1)
  expect_equal(params$y, "hello")
  expect_true("symbol_names" %in% names(attributes(params)))
  expect_equal(unname(attr(params, "symbol_names")[["x"]]), "1")
})

test_that("capture_params uses each unnamed argument's own expression as its name", {
  respath <- "some/path"
  params <- capture_params(respath)

  expect_true("respath" %in% names(params))
  expect_equal(params$respath, "some/path")
})

test_that("build_param_folder encodes recognized parameters into a folder name, one folder per combination", {
  params <- capture_params(
    uka = "uka_data", ppi_network = "ppi_networkv12",
    spec_cutoff = 0.7, rank_uka_abs = TRUE, b = 1.5
  )

  folder <- build_param_folder(params)
  expect_true(grepl("thr_0.7", folder))
  expect_true(grepl("ppi_networkv12", folder))
  expect_true(grepl("ukaabs1_b1.5", folder))
})

test_that("build_param_folder handles a spec_cutoff vector (a sweep captured for one folder) the same as before", {
  params <- capture_params(
    uka = "uka_data", ppi_network = "ppi_networkv12",
    spec_cutoff = c(0.7, 0.9), rank_uka_abs = TRUE, b = 1.5
  )

  expect_true(grepl("thr_0.7-0.9", build_param_folder(params)))
})

test_that("build_param_folder omits spec_cutoff/perc_cutoff from the name when they're 0", {
  params <- capture_params(
    uka = "uka_data", ppi_network = "ppi_networkv12",
    spec_cutoff = 0, perc_cutoff = 0, rank_uka_abs = TRUE, b = 1.5
  )

  folder <- build_param_folder(params)
  expect_false(grepl("thr_", folder))
  expect_false(grepl("perc", folder))
  expect_true(grepl("ukaabs1_b1.5", folder))
})

test_that("build_param_folder still shows spec_cutoff/perc_cutoff when they're not 0", {
  params <- capture_params(
    uka = "uka_data", ppi_network = "ppi_networkv12",
    spec_cutoff = 0.7, perc_cutoff = 0.5, rank_uka_abs = TRUE, b = 1.5
  )

  folder <- build_param_folder(params)
  expect_true(grepl("thr_0.7", folder))
  expect_true(grepl("perc0.5", folder))
})

test_that("build_param_folder nests under a uka-named subfolder only when sens is present (paired analyses)", {
  # Real callers always pass data frames for uka/sens/ppi_network (never a
  # scalar), so capture_params() records their variable name, not their
  # value -- mirror that with data.frame placeholders here rather than bare
  # strings (which are themselves scalars and would trigger the "show the
  # value" path instead).
  uka_data <- data.frame(name = "K1")
  sens_data <- data.frame(name = "S1")
  ppi <- data.frame(head = "A", tail = "B", cost = 0.1)

  unpaired <- capture_params(uka = uka_data, ppi_network = ppi, rank_uka_abs = TRUE, b = 1)
  paired <- capture_params(uka = uka_data, sens = sens_data, ppi_network = ppi, rank_uka_abs = TRUE, b = 1)

  expect_false(grepl("/", build_param_folder(unpaired)))
  expect_true(startsWith(build_param_folder(paired), "uka_data/"))
})

test_that("prepare_run_params captures params, creates a parameter-encoded folder, and writes params.csv", {
  respath <- file.path(tempdir(), "prepare_run_params_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  uka <- data.frame(name = "K1")
  ppi_networkv12 <- data.frame(head = "A", tail = "B", cost = 0.1)

  params <- prepare_run_params(
    respath = respath,
    uka = uka,
    ppi_network = ppi_networkv12,
    spec_cutoff = c(0.7, 0.9),
    rank_uka_abs = TRUE,
    b = 1.5
  )

  expect_true(dir.exists(params$respath))
  expect_true(file.exists(file.path(params$respath, "params.csv")))
  expect_true(grepl("thr_0.7-0.9", params$respath))
  expect_true(grepl("ukaabs1_b1.5", params$respath))
  expect_identical(params$uka, uka)
})

test_that("prepare_run_params gives each distinct spec_cutoff its own folder -- conditions share it via per-file naming, not per-folder", {
  respath <- file.path(tempdir(), "prepare_run_params_flatten_test")
  unlink(respath, recursive = TRUE)
  dir.create(respath)

  uka <- data.frame(name = "K1")
  ppi_networkv12 <- data.frame(head = "A", tail = "B", cost = 0.1)
  conditions <- paste0("cond", 1:4)
  spec_cutoffs <- c(0.5, 0.7)

  # One prepare_run_params() call per spec_cutoff (2 folders); each folder is
  # then reused for all 4 conditions' generate_*_network(res.path=..., write=TRUE)
  # calls, which distinguish conditions by filename (nodes_<condition>_spec*.csv),
  # not by folder -- matching kinograte_pg_pcsf()'s existing per-file convention.
  folders <- purrr::map_chr(spec_cutoffs, function(sc) {
    prepare_run_params(
      respath = respath, uka = uka, ppi_network = ppi_networkv12,
      spec_cutoff = sc, rank_uka_abs = TRUE, b = 1
    )$respath
  })

  expect_equal(length(unique(folders)), length(spec_cutoffs))
})

test_that("spec_suffix omits the suffix for 0 and is vectorized", {
  expect_equal(spec_suffix(0), "")
  expect_equal(spec_suffix(0.7), "_spec0.7")
  expect_equal(spec_suffix(c(0, 0.7)), c("", "_spec0.7"))
})

test_that("save_params writes a clean parameter/value table and drops uka_fam/respath", {
  respath <- file.path(tempdir(), "save_params_test")
  dir.create(respath, showWarnings = FALSE)

  params <- list(a = 1, uka_fam = "shouldbedropped", respath = respath)
  attr(params, "symbol_names") <- c(a = "1", uka_fam = "fam_df", respath = "\"some/path\"")

  df <- save_params(params, respath = respath)
  expect_false("uka_fam" %in% df$parameter)
  expect_false("respath" %in% df$parameter)
  expect_true(file.exists(file.path(respath, "params.csv")))
})
