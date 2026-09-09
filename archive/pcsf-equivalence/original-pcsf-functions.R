# Archived original (igraph-based) PCSF path.
#
# Removed from R/pcsf-core.R -- kept here, not deleted outright, so the
# equivalence proof in verify_equivalence.R stays runnable. See
# docs/adr/0005-remove-original-pcsf-path.md for why this was
# removed from the package, and verify_equivalence.R for how to re-derive
# the equivalence result (rather than trust this comment) if the fast path
# ever changes.
#
# Verbatim copies of kinograte_pg() and PCSF_rand() as they existed in
# networkGen before removal -- unchanged except calling networkGen:::call_sr()
# explicitly, since call_sr() is an unexported package internal and this file
# is sourced standalone, outside the package namespace.

PCSF_rand <- function(ppi, terminals, n = 8, r = 0.1, w = 2, b = 1, mu = 0.0005, dummies) {
  if (missing(ppi)) {
    stop("Need to specify an interaction network \"ppi\".")
  }
  if (!methods::is(ppi, "igraph")) {
    stop("The interaction network \"ppi\" must be an igraph object.")
  }
  if (missing(terminals)) {
    stop("  Need to provide terminal nodes as a named numeric vector, \n",
      "  where node names must be same as in the interaction network.")
  }
  if (is.null(names(terminals))) {
    stop("  The terminal nodes must be provided as a named numeric vector, \n",
      "  where node names must be same as in the interaction network.")
  }

  terminal_names <- names(terminals)
  terminal_values <- as.numeric(terminals)

  node_names <- igraph::V(ppi)$name
  node_prz <- vector(mode = "numeric", length = length(node_names))
  index <- match(terminal_names, node_names)
  percent <- signif((length(index) - sum(is.na(index))) / length(index) * 100, 4)
  if (percent < 5) {
    stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
  }
  cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
  terminal_names <- terminal_names[!is.na(index)]
  terminal_values <- terminal_values[!is.na(index)]
  index <- index[!is.na(index)]
  node_prz[index] <- terminal_values

  if (missing(dummies) || is.null(dummies) || is.na(dummies)) {
    dummies <- terminal_names
  }

  cat("  Solving the PCSF by adding random noise to the edge costs...\n")

  node_degrees <- igraph::degree(ppi)
  hub_penalization <- -mu * node_degrees

  node_prizes <- b * node_prz
  index <- which(node_prizes == 0)
  node_prizes[index] <- hub_penalization[index]

  edges <- igraph::ends(ppi, es = igraph::E(ppi))
  from <- c(rep("DUMMY", length(dummies)), edges[, 1])
  to <- c(dummies, edges[, 2])

  all_nodes <- NULL

  for (i in 1:n) {
    cat("Run", i, "of", n, "\n")
    random_vals <- stats::runif(length(igraph::E(ppi)), 0, r)

    cost <- c(rep(w, length(dummies)), igraph::E(ppi)$weight + igraph::E(ppi)$weight * random_vals)

    output <- networkGen:::call_sr(from, to, cost, node_names, node_prizes)

    edge <- data.frame(as.character(output[[1]]), as.character(output[[2]]))
    colnames(edge) <- c("source", "target")
    edge <- edge[which(edge[, 1] != "DUMMY"), ]
    edge <- edge[which(edge[, 2] != "DUMMY"), ]

    graph <- igraph::graph.data.frame(edge, directed = FALSE)
    assign(paste0("graph_", i), graph)
    all_nodes <- c(all_nodes, igraph::V(graph)$name)
  }

  node_frequency <- table(all_nodes)
  print(utils::head(node_frequency, 10))

  node_names <- names(node_frequency)
  node_prizes <- as.numeric(node_frequency)

  adj_matrix <- matrix(0, length(node_names), length(node_names))
  colnames(adj_matrix) <- node_names
  rownames(adj_matrix) <- node_names
  for (i in 1:n) {
    graph <- get(paste0("graph_", i))
    edges <- igraph::ends(graph, es = igraph::E(graph))
    x <- match(edges[, 1], node_names)
    y <- match(edges[, 2], node_names)
    if (length(x) > 0 & length(y) > 0) {
      for (j in seq_along(x)) {
        if (x[j] >= y[j]) {
          k <- x[j]
          l <- y[j]
        } else {
          k <- y[j]
          l <- x[j]
        }
        adj_matrix[k, l] <- adj_matrix[k, l] + 1
      }
    }
  }

  if (sum(adj_matrix) != 0) {
    subnet <- igraph::graph_from_adjacency_matrix(adj_matrix, weighted = TRUE, mode = "undirected")
    index <- match(igraph::V(subnet)$name, node_names)
    igraph::V(subnet)$prize <- node_prizes[index]

    igraph::V(subnet)$type <- "Steiner"
    index <- match(terminal_names, igraph::V(subnet)$name)
    index <- index[!is.na(index)]
    igraph::V(subnet)$type[index] <- "Terminal"

    class(subnet) <- c("PCSF", "igraph")

    return(subnet)
  } else {
    stop("  Subnetwork can not be identified for a given parameter set.\n",
      "  Provide a compatible b or mu value with your terminal prize list...\n\n")
  }
}

kinograte_pg <- function(df, ppi_network, maintitle, n = 8, w = 10, r = 0.1, b = 1.5,
                          mu = 0.005, cluster = TRUE, seed = NULL, res.path, spec_cutoff, condition, write) {
  maintitle <- paste0(condition, " - Network with specificity cutoff = ", spec_cutoff, ", Number of nodes = ", nrow(df))

  terms <- df$prize
  names(terms) <- df$name
  ppi_net <- base::as.data.frame(ppi_network)
  ppi <- PCSF::construct_interactome(ppi_net)
  if (!is.null(seed)) {
    base::set.seed(seed)
  }

  subnet <- tryCatch(
    {
      PCSF_rand(ppi, terms, n = n, w = w, r = r, b = b, mu = mu)
    },
    error = function(e) {
      message("PCSF_rand error: ", e$message)
      return(NULL)
    }
  )

  if (is.null(subnet) || igraph::vcount(subnet) == 0 || igraph::ecount(subnet) == 0) {
    message("PCSF_rand returned empty network.")
    return(NULL)
  }
  edges <- igraph::as_data_frame(subnet, what = "edges")

  nodes <- dplyr::tibble(Protein = c(
    base::unique(edges$from),
    base::unique(edges$to)
  )) %>% dplyr::distinct()

  nodes <- dplyr::left_join(nodes, df, by = c(Protein = "name")) %>%
    dplyr::mutate(
      type = ifelse(is.na(.data$type), "Hidden", .data$type),
      prize = ifelse(is.na(.data$prize), 2, .data$prize)
    )

  if (!"degree" %in% colnames(nodes)) {
    node_degrees <- dplyr::tibble(
      Protein = igraph::V(subnet)$name,
      degree = igraph::degree(subnet)
    )
    nodes <- dplyr::left_join(nodes, node_degrees, by = "Protein")
  }

  subnet_nodes <- igraph::V(subnet)$name
  missing_nodes <- setdiff(names(terms), subnet_nodes)
  if (length(missing_nodes) > 0 && write) {
    missing_df <- df %>% dplyr::filter(.data$name %in% missing_nodes)
    readr::write_csv(missing_df, paste0(res.path, "/missing_nodes_", condition, networkGen:::spec_suffix(spec_cutoff), ".csv"))
  } else {
    missing_df <- NULL
  }

  kinograte_nodes <- nodes %>%
    as.data.frame() %>%
    dplyr::select(-dplyr::any_of("title"))

  if (write) {
    readr::write_csv(
      kinograte_nodes,
      paste0(res.path, "/nodes_", condition, networkGen:::spec_suffix(spec_cutoff), ".csv")
    )
  }

  if (!"pathway" %in% colnames(nodes)) {
    nodes$title <- if ("LogFC_all" %in% colnames(nodes)) {
      paste0("LFC = ", nodes$LogFC_all)
    } else {
      paste0("LFC = ", round(nodes$LogFC, 2))
    }
  } else {
    nodes$title <- nodes$pathway
  }

  if (cluster) {
    wc <- suppressWarnings(igraph::edge.betweenness.community(subnet,
      weights = igraph::E(subnet)$weight, directed = FALSE,
      bridges = TRUE
    ))
    wc_df <- dplyr::tibble(id = wc$names, cluster = wc$membership)

    return(list(
      network = subnet, nodes = nodes, edges = edges, missing = missing_df,
      wc_df = wc_df, maintitle = maintitle
    ))
  } else {
    return(list(network = subnet, nodes = nodes, edges = edges, missing_nodes = missing_df, maintitle = maintitle))
  }
}
