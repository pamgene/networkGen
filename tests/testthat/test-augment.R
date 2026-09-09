test_that("augment_by_threshold_steps expands from a core subgraph (steps = 1)", {
  # Chain A-B-C-D-E, all edges cost 0.1 except D-E which is expensive (0.9)
  full_edges <- data.frame(
    from = c("A", "B", "C", "D"),
    to = c("B", "C", "D", "E"),
    cost = c(0.1, 0.1, 0.1, 0.9)
  )
  full_graph <- igraph::graph_from_data_frame(full_edges, directed = FALSE)

  core_graph <- igraph::graph_from_data_frame(
    data.frame(from = "B", to = "C", cost = 0.1),
    directed = FALSE, vertices = c("B", "C")
  )

  aug <- augment_by_threshold_steps(core_graph, full_graph, cost_threshold = 0.1, steps = 1)
  aug_nodes <- igraph::V(aug)$name

  # From {B, C}, one cheap-edge step should pull in A (via A-B) and D (via C-D), not E (expensive edge)
  expect_true(all(c("A", "B", "C", "D") %in% aug_nodes))
  expect_false("E" %in% aug_nodes)
})

test_that("augment_by_threshold_steps works across multiple steps (regression test for the incident_edges name-shadowing bug)", {
  # Chain A-B-C-D-E-F, all cheap edges
  full_edges <- data.frame(
    from = c("A", "B", "C", "D", "E"),
    to = c("B", "C", "D", "E", "F"),
    cost = rep(0.1, 5)
  )
  full_graph <- igraph::graph_from_data_frame(full_edges, directed = FALSE)
  core_graph <- igraph::graph_from_data_frame(
    data.frame(from = "C", to = "D", cost = 0.1),
    directed = FALSE, vertices = c("C", "D")
  )

  # steps = 2 previously errored ("attempt to apply non-function") because the
  # loop reassigned `incident_edges` (a local variable) over the function of
  # the same name on the first iteration.
  expect_no_error({
    aug <- augment_by_threshold_steps(core_graph, full_graph, cost_threshold = 0.1, steps = 2)
  })
  aug_nodes <- igraph::V(aug)$name
  expect_true(all(c("A", "B", "C", "D", "E", "F") %in% aug_nodes))
})
