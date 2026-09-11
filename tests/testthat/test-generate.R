toy_ppi <- data.frame(
  head = c("K1", "K2", "K3", "S1"),
  tail = c("K2", "K3", "S1", "S2"),
  cost = c(0.1, 0.1, 0.2, 0.1)
)

toy_uka <- data.frame(name = c("K1", "K2"), prize = c(0.9, 0.8), type = "Kinase", LogFC = c(1.5, -2.0))
toy_sens <- data.frame(name = c("S1", "S2"), prize = c(0.7, 0.6), type = "Sensitivity", LogFC = c(-1.5, -2.0))

# A hand-verifiable topology: terminals A, B, C. Two candidate paths connect
# B to C -- via D (cheap, total cost 0.2) or via E-F (expensive, total cost
# 0.3). PCSF should always prefer the cheap path: PCSF_rand_pg()'s random
# noise only ever *increases* edge cost (base_cost * (1 + noise), noise in
# [0, r], default r = 0.1), never decreases it -- so the cheap path's worst
# case (0.1*1.1 + 0.1*1.1 = 0.22) is always below the expensive path's best
# case (exactly 0.3, no noise applied) -- deterministic across every
# randomized run, not just "usually right." Expect D in the result, E and F
# absent. At b = 10 (matches the value already validated above), collecting
# all three terminals' prizes only pays off if they're merged into one
# connected component (sharing a single costly dummy-root edge, w = 10,
# instead of paying it three times) -- so the solver is economically forced
# to connect them, not just inclined to.
path_ppi <- data.frame(
  head = c("A", "B", "D", "B", "E", "F"),
  tail = c("B", "D", "C", "E", "F", "C"),
  cost = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
)
path_uka <- data.frame(name = c("A", "B", "C"), prize = c(0.9, 0.9, 0.9), type = "Kinase", LogFC = c(1, 1, 1))

test_that("generate_kinase_network picks the cheaper of two connecting paths, not just any connected subnetwork (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  result <- generate_kinase_network(
    uka = path_uka, ppi_network = path_ppi, spec_cutoff = 0, b = 10,
    condition = "path_test", write = FALSE, seed = 12345
  )

  expect_s3_class(result, "networkGen_result")

  nodes_found <- result$nodes$Protein
  expect_true(all(c("A", "B", "C", "D") %in% nodes_found))
  expect_false(any(c("E", "F") %in% nodes_found))

  has_edge <- function(edges, x, y) {
    any((edges$from == x & edges$to == y) | (edges$from == y & edges$to == x))
  }
  expect_true(has_edge(result$edges, "A", "B"))
  expect_true(has_edge(result$edges, "B", "D"))
  expect_true(has_edge(result$edges, "D", "C"))
  expect_false(has_edge(result$edges, "B", "E"))
  expect_false(has_edge(result$edges, "E", "F"))
  expect_false(has_edge(result$edges, "F", "C"))
})

test_that("generate_kinase_network degrades to NULL, not an error, when PCSF's compiled code is unavailable", {
  skip_if(pcsf_functional(), "real PCSF is installed -- this test only covers the stub-PCSF failure path")

  result <- generate_kinase_network(
    uka = toy_uka, ppi_network = toy_ppi, spec_cutoff = 0, b = 1.5,
    condition = "test", write = FALSE
  )
  expect_null(result)
})

test_that("generate_kinase_network builds a real network end-to-end (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  result <- generate_kinase_network(
    uka = toy_uka, ppi_network = toy_ppi, spec_cutoff = 0, b = 10,
    condition = "test", write = FALSE
  )
  expect_s3_class(result, "networkGen_result")
  expect_s3_class(result$network, "igraph")
  expect_true(all(c("nodes", "edges", "wc_df", "maintitle", "params") %in% names(result)))
})

test_that("result edges carry the reference PPI's cost, matched regardless of edge orientation (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  # Distinct per-edge costs so a wrong match would show up as a wrong value,
  # not just coincidentally equal. PCSF's output edge orientation is
  # arbitrary, so some result edges will be (tail, head) relative to this
  # table -- the cost lookup must still find them.
  ppi <- data.frame(
    head = c("A", "B", "D", "B", "E", "F"),
    tail = c("B", "D", "C", "E", "F", "C"),
    cost = c(0.11, 0.22, 0.33, 0.44, 0.55, 0.66)
  )
  uka <- data.frame(name = c("A", "B", "C"), prize = c(0.9, 0.9, 0.9), type = "Kinase", LogFC = c(1, 1, 1))

  result <- generate_kinase_network(
    uka = uka, ppi_network = ppi, spec_cutoff = 0, b = 10,
    condition = "cost_test", write = FALSE, seed = 12345
  )

  expect_true("cost" %in% colnames(result$edges))
  expect_false(anyNA(result$edges$cost))

  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  want <- stats::setNames(ppi$cost, key(ppi$head, ppi$tail))
  expect_equal(unname(result$edges$cost), unname(want[key(result$edges$from, result$edges$to)]))
})

test_that("generate_paired_network builds a real network end-to-end (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  result <- generate_paired_network(
    uka = toy_uka, sens = toy_sens, ppi_network = toy_ppi, spec_cutoff = 0, b = 10,
    condition = "test", write = FALSE
  )
  expect_s3_class(result, "networkGen_result")
})

test_that("written nodes/edges/wc_df filenames omit _spec when spec_cutoff is 0, include it otherwise (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  respath_zero <- file.path(tempdir(), "spec_suffix_zero_test")
  unlink(respath_zero, recursive = TRUE)
  dir.create(respath_zero)
  generate_kinase_network(
    uka = path_uka, ppi_network = path_ppi, spec_cutoff = 0, b = 10,
    condition = "cond", write = TRUE, res.path = respath_zero, seed = 12345
  )
  expect_true(file.exists(file.path(respath_zero, "nodes_cond.csv")))
  expect_true(file.exists(file.path(respath_zero, "edges_cond.csv")))
  expect_true(file.exists(file.path(respath_zero, "wc_df_cond.csv")))
  expect_false(file.exists(file.path(respath_zero, "nodes_cond_spec0.csv")))

  respath_nonzero <- file.path(tempdir(), "spec_suffix_nonzero_test")
  unlink(respath_nonzero, recursive = TRUE)
  dir.create(respath_nonzero)
  generate_kinase_network(
    uka = path_uka, ppi_network = path_ppi, spec_cutoff = 0.7, b = 10,
    condition = "cond", write = TRUE, res.path = respath_nonzero, seed = 12345
  )
  expect_true(file.exists(file.path(respath_nonzero, "nodes_cond_spec0.7.csv")))
  expect_true(file.exists(file.path(respath_nonzero, "edges_cond_spec0.7.csv")))
  expect_true(file.exists(file.path(respath_nonzero, "wc_df_cond_spec0.7.csv")))
})

test_that("ppi_network defaults to the bundled ppi_network_human_filtered_v12.5", {
  expect_identical(formals(generate_kinase_network)$ppi_network, as.symbol("ppi_network_human_filtered_v12.5"))
  expect_identical(formals(generate_paired_network)$ppi_network, as.symbol("ppi_network_human_filtered_v12.5"))
})

# This is the one deliberately "heavy" test in the suite -- it runs PCSF
# over the full ~1.2M-edge bundled ppi_network_human_filtered_v12.5 (a few
# seconds), rather than a tiny hand-built graph. Everything else in this
# file uses small toy networks precisely so the rest of the suite stays
# fast; this test exists only to confirm the default-argument wiring and
# real-scale behavior work, not to be a template for how every PCSF test
# should be written.
test_that("generate_kinase_network works against the real bundled reference network with real gene symbols, using the default ppi_network (needs real PCSF)", {
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  set.seed(42)
  real_genes <- sample(unique(c(ppi_network_human_filtered_v12.5$head, ppi_network_human_filtered_v12.5$tail)), 15)
  uka_real <- data.frame(name = real_genes, prize = stats::runif(15, 0.5, 1), type = "Kinase", LogFC = stats::rnorm(15))

  # 15 random genes aren't a biologically connected set: at the default
  # w = 2 most get parked on the artificial root and drop out of the node
  # list. This test is about real-scale PCSF wiring, not the defaults, so
  # it uses b = 10, w = 10 -- high enough w that terminals must route
  # through the network to be kept.
  result <- generate_kinase_network(uka = uka_real, spec_cutoff = 0, b = 10, w = 10, condition = "real_network_test", write = FALSE)

  expect_s3_class(result, "networkGen_result")
  expect_true(all(real_genes %in% result$nodes$Protein))
})
