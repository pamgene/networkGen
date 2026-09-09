#' Augment a PCSF result with additional threshold-connected neighbors
#'
#' Starting from a solved PCSF subgraph's nodes, iteratively expands the node
#' set by following edges in `full_graph` whose cost is at or below
#' `cost_threshold`, for `steps` iterations, then returns the induced
#' subgraph of `full_graph` on the final node set.
#'
#' Two bugs fixed while migrating this from `Network_generation`'s
#' `R/kinograte_PG.R` (unused there, so never actually exercised): the
#' starting node set referenced an undefined `core_nodes0` (typo for
#' `core_nodes`), and the loop reassigned `incident_edges` (a local variable)
#' over the function of the same name, which would break on any `steps > 1`
#' call.
#'
#' @param pcsf_graph An `igraph` object -- the solved PCSF subgraph whose
#'   node set is the starting point.
#' @param full_graph An `igraph` object -- the full PPI network to expand into,
#'   with an edge attribute `cost`.
#' @param cost_threshold Maximum edge cost to follow when expanding. Default 0.1.
#' @param steps Number of expansion iterations. Default 1.
#'
#' @return The induced subgraph of `full_graph` on the expanded node set.
#' @export
augment_by_threshold_steps <- function(pcsf_graph, full_graph, cost_threshold = 0.1, steps = 1) {
  core_nodes <- igraph::V(pcsf_graph)$name
  print(paste0("number of input nodes: ", length(unique(core_nodes))))

  current_nodes <- core_nodes

  for (step in 1:steps) {
    print(paste0("Step ", step, " of ", steps))

    current_vids <- which(igraph::V(full_graph)$name %in% current_nodes)

    incident <- igraph::incident_edges(full_graph, current_vids, mode = "all")
    incident <- unique(unlist(incident))

    strong_edges <- incident[igraph::E(full_graph)[incident]$cost <= cost_threshold]

    strong_verts <- igraph::ends(full_graph, strong_edges, names = TRUE)
    strong_verts <- unique(as.vector(strong_verts))

    current_nodes <- unique(c(current_nodes, strong_verts))

    print(paste0("After step ", step, ", nodes count: ", length(current_nodes)))
  }

  augmented_graph <- igraph::induced_subgraph(full_graph, which(igraph::V(full_graph)$name %in% current_nodes))
  print(paste0("Number of nodes in final augmented graph: ", length(unique(igraph::V(augmented_graph)$name))))

  augmented_graph
}
