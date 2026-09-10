#' Call PCSF's compiled Steiner-forest solver
#'
#' Thin wrapper around the compiled `PCSF` package's `_PCSF_call_sr` C++
#' entry point. Internal -- not exported.
#'
#' `.Call(..., PACKAGE = "PCSF")` only works if `PCSF`'s namespace (and so
#' its compiled DLL) is already loaded -- unlike calling an R-level function
#' via `PCSF::something()`, it does *not* trigger lazy loading itself. Listing
#' `PCSF` in `Imports:` alone is not enough to guarantee that either: R only
#' eagerly loads an Imports package's namespace when this package's own
#' `NAMESPACE` has a real `import()`/`importFrom()` entry for it. The
#' `@importFrom PCSF construct_interactome` below exists specifically to
#' create that entry, so `library(networkGen)` alone loads `PCSF` before any
#' `call_sr()` call can run -- this was missed originally (every test passed
#' under `pkgload::load_all()`, which eagerly loads all `Imports:` regardless
#' of `NAMESPACE` directives, masking the gap) and only caught by testing a
#' real installed-and-`library()`d session.
#'
#' @param from,to,cost Edge list (character/character/numeric), including
#'   `"DUMMY"` root edges.
#' @param node_names Character vector of all node names.
#' @param node_prizes Numeric vector of node prizes, same order as `node_names`.
#'
#' @return List with the solved edge list (as returned by the C++ solver).
#' @importFrom PCSF construct_interactome
#' @keywords internal
call_sr <- function(from, to, cost, node_names, node_prizes) {
  .Call("_PCSF_call_sr", PACKAGE = "PCSF", from, to, cost, node_names, node_prizes)
}

#' Solve PCSF with randomized edge costs
#'
#' Runs the prize-collecting Steiner forest algorithm `n` times, each time
#' adding random noise to edge costs, and returns the union of all runs as a
#' consensus subnetwork. Dataframe-based -- skips `PCSF::construct_interactome()`.
#' A slower igraph-based implementation existed alongside this one and was
#' removed after confirming the two were bit-exact equivalent (including the
#' raw random values driving the solver) -- see
#' `archive/pcsf-equivalence/` and `docs/adr/0005-remove-original-pcsf-path.md`.
#'
#' @param edges_df Data frame with columns `head`, `tail`, `cost` (the PPI network).
#' @param terminals Named numeric vector: names are terminal node names
#'   (matching `edges_df`), values are their prizes.
#' @param n Number of randomized runs to solve and combine into a consensus
#'   network -- each run adds different random noise to edge costs (see `r`),
#'   and the final result is the union across all `n` runs. Higher `n` gives
#'   a more stable result at the cost of speed; the default (8) is often
#'   lowered for faster iteration while developing/testing. Default 8.
#' @param r Fraction of random noise added to each edge cost per run (0-1),
#'   used to build the `n`-run consensus (see `n`). Default 0.1.
#' @param w The cost of connecting a terminal node straight to an artificial
#'   root, i.e. the price of "parking" a terminal on its own rather than
#'   reaching it through the network. Low `w` lets isolated terminals in
#'   cheaply (more, smaller fragments); high `w` means a terminal is only
#'   worth keeping if the network already passes near it (fewer, more
#'   connected components -- isolated terminals get dropped). A lone
#'   terminal survives when `b * prize > w`. Default 2 (PCSF's own default).
#' @param b Multiplier applied to every node prize. Raising `b` makes prizes
#'   outweigh edge costs, so more input nodes -- and more connecting
#'   (Steiner) nodes -- get pulled in: a bigger network. It also rescues
#'   lone terminals, since `b * prize > w` is their break-even. If PCSF
#'   finds no subnetwork worth building, raise `b` (or lower `w`). Default 2.
#' @param mu Hub-penalization strength: non-terminal (Steiner) nodes are
#'   penalized by `mu * their degree in the interactome`, discouraging the
#'   algorithm from routing through very high-degree "hub" genes just because
#'   they happen to connect to many things. Default 0.0005.
#' @param dummies Optional character vector of dummy-root target names;
#'   defaults to `names(terminals)`.
#'
#' @return A list of class `c("PCSF_df", "list")` with `edges` (a data
#'   frame with columns `from`, `to`, `weight` -- the number of the `n`
#'   runs the edge appeared in -- and `cost` -- the reference PPI's
#'   interaction cost for that edge, carried over from `edges_df`),
#'   `nodes` (name/prize/type data frame), and `n_runs`.
#' @export
PCSF_rand_pg <- function(edges_df, terminals, n = 8, r = 0.1,
                          w = 2, b = 2, mu = 0.0005, dummies = NULL) {
  if (missing(edges_df)) {
    stop("Need to specify edge dataframe with columns: head, tail, cost")
  }
  if (!all(c("head", "tail", "cost") %in% colnames(edges_df))) {
    stop("edges_df must have columns: head, tail, cost")
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

  # === REPLICATE construct_interactome EXACTLY ===
  node_names <- unique(c(as.character(edges_df[, 1]), as.character(edges_df[, 2])))
  temp_graph <- igraph::graph.data.frame(
    edges_df[, 1:2],
    vertices = node_names,
    directed = FALSE
  )
  igraph::E(temp_graph)$weight <- as.numeric(edges_df[, 3])
  temp_graph <- igraph::simplify(temp_graph)

  node_degrees <- igraph::degree(temp_graph)
  node_degrees <- node_degrees[node_names]
  edges_simplified <- igraph::ends(temp_graph, es = igraph::E(temp_graph))
  edge_weights_simplified <- igraph::E(temp_graph)$weight

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

  hub_penalization <- -mu * node_degrees

  node_prizes <- b * node_prz
  index <- which(node_prizes == 0)
  node_prizes[index] <- hub_penalization[index]

  from <- c(rep("DUMMY", length(dummies)), edges_simplified[, 1])
  to <- c(dummies, edges_simplified[, 2])
  base_costs <- edge_weights_simplified

  all_nodes <- NULL

  for (i in 1:n) {
    cat("Run", i, "of", n, "\n")
    random_noise <- stats::runif(length(base_costs), 0, r)

    cost <- c(
      rep(w, length(dummies)),
      base_costs + (base_costs * random_noise)
    )

    output <- call_sr(from, to, cost, node_names, node_prizes)

    edge <- data.frame(
      as.character(output[[1]]),
      as.character(output[[2]]),
      stringsAsFactors = FALSE
    )
    colnames(edge) <- c("source", "target")

    edge <- edge[which(edge[, 1] != "DUMMY"), ]
    edge <- edge[which(edge[, 2] != "DUMMY"), ]

    assign(paste0("graph_", i), edge)
    all_nodes <- c(all_nodes, edge$source, edge$target)
  }

  node_frequency <- table(all_nodes)
  node_names_out <- names(node_frequency)
  node_prizes_out <- as.numeric(node_frequency)

  adj_matrix <- matrix(0, length(node_names_out), length(node_names_out))
  colnames(adj_matrix) <- node_names_out
  rownames(adj_matrix) <- node_names_out

  for (i in 1:n) {
    graph <- get(paste0("graph_", i))
    edges <- graph

    x <- match(edges[, 1], node_names_out)
    y <- match(edges[, 2], node_names_out)

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
    edge_indices <- which(adj_matrix > 0, arr.ind = TRUE)
    result_edges <- data.frame(
      from = node_names_out[edge_indices[, 1]],
      to = node_names_out[edge_indices[, 2]],
      weight = adj_matrix[edge_indices],
      stringsAsFactors = FALSE
    )

    # Reattach the reference PPI's interaction cost to each selected edge.
    # The solver returns only endpoint-name pairs, so match them back to the
    # *simplified* interactome's costs -- the exact values the solver
    # searched over, so any parallel-edge collapsing igraph::simplify() did
    # is already accounted for. The match is direction-agnostic: an
    # undirected edge is stored (A,B) here but may be (B,A) in the PPI, so
    # both sides are keyed on the pair sorted into a fixed order
    # (pmin = the name that sorts first, pmax = the name that sorts last).
    ppi_key <- paste(
      pmin(edges_simplified[, 1], edges_simplified[, 2]),
      pmax(edges_simplified[, 1], edges_simplified[, 2]),
      sep = "||"
    )
    ppi_cost <- edge_weights_simplified[!duplicated(ppi_key)]
    names(ppi_cost) <- ppi_key[!duplicated(ppi_key)]
    result_key <- paste(
      pmin(result_edges$from, result_edges$to),
      pmax(result_edges$from, result_edges$to),
      sep = "||"
    )
    result_edges$cost <- unname(ppi_cost[result_key])

    index <- match(node_names_out, node_names_out)
    result_nodes <- data.frame(
      name = node_names_out,
      prize = node_prizes_out[index],
      type = ifelse(node_names_out %in% terminal_names, "Terminal", "Steiner"),
      stringsAsFactors = FALSE
    )
    result <- list(
      edges = result_edges,
      nodes = result_nodes,
      n_runs = n
    )

    class(result) <- c("PCSF_df", "list")
    return(result)
  } else {
    stop("  Subnetwork can not be identified for a given parameter set.\n",
      "  Provide a compatible b or mu value with your terminal prize list...\n\n")
  }
}

#' Build a PCSF network from prize data and a PPI network
#'
#' Calls [PCSF_rand_pg()] directly on `ppi_network` without going through
#' `PCSF::construct_interactome()`. Used by [generate_paired_network()] /
#' [generate_kinase_network()].
#'
#' @param df Data frame of terminal nodes with columns `name`, `prize`
#'   (and typically `type`, `LogFC`), e.g. from [uka_top()] / [sens_top()].
#' @param ppi_network Data frame with columns `head`, `tail`, `cost` -- the
#'   PPI network to search over. Always required explicitly; never defaults
#'   to a global.
#' @param maintitle Ignored if given; recomputed from `condition`/`spec_cutoff`.
#' @param n,w,r,b,mu Passed to [PCSF_rand_pg()].
#' @param cluster If `TRUE` (default), cluster the resulting network with
#'   `igraph::edge.betweenness.community()` and include `wc_df`.
#' @param seed Optional RNG seed for reproducibility.
#' @param res.path Output folder, used only when `write = TRUE`.
#' @param spec_cutoff Specificity cutoff, recorded in `maintitle` and (when
#'   `write = TRUE`) output file names.
#' @param condition Condition/comparison label, recorded in `maintitle` and
#'   (when `write = TRUE`) output file names.
#' @param write If `TRUE`, write `nodes_<condition>_spec<spec_cutoff>.csv`,
#'   `edges_<condition>_spec<spec_cutoff>.csv`, `wc_df_<condition>_spec<spec_cutoff>.csv`
#'   (if `cluster = TRUE`), and (if any) `missing_nodes_<condition>_spec<spec_cutoff>.csv`
#'   to `res.path`. Pass [prepare_run_params()]'s resulting `respath` as
#'   `res.path` to get a parameter-encoded, provenance-tracked folder rather
#'   than an arbitrary caller-chosen one.
#'
#' @return A list: `network` (igraph), `nodes`, `edges`, `missing`/`missing_nodes`
#'   (data frames), `wc_df` (if `cluster = TRUE`), `maintitle`. `NULL` if PCSF
#'   fails to find a subnetwork.
#' @export
kinograte_pg_pcsf <- function(df, ppi_network, maintitle, n = 8, w = 2, r = 0.1, b = 2,
                               mu = 0.005, cluster = TRUE, seed = NULL, res.path, spec_cutoff,
                               condition, write) {
  print("Building PCSF network...")
  maintitle <- paste0(
    condition, " - Network with specificity cutoff = ",
    spec_cutoff, ", Number of nodes = ", nrow(df)
  )

  terms <- df$prize
  names(terms) <- df$name

  if (!is.null(seed)) {
    set.seed(seed)
  }

  ppi_network <- base::as.data.frame(ppi_network)

  subnet <- tryCatch(
    {
      PCSF_rand_pg(
        edges_df = ppi_network,
        terminals = terms,
        n = n,
        w = w,
        r = r,
        b = b,
        mu = mu
      )
    },
    error = function(e) {
      message("PCSF_rand_pg error: ", e$message)
      return(NULL)
    }
  )

  if (is.null(subnet)) {
    return(NULL)
  }

  edges <- subnet$edges

  nodes <- dplyr::tibble(
    Protein = c(base::unique(edges$from), base::unique(edges$to))
  ) %>%
    dplyr::distinct()

  nodes <- dplyr::left_join(nodes, df, by = c("Protein" = "name")) %>%
    dplyr::mutate(
      type = ifelse(is.na(.data$type), "Hidden", .data$type),
      prize = ifelse(is.na(.data$prize), 2, .data$prize)
    )

  if (!"degree" %in% colnames(nodes)) {
    degree_from <- table(edges$from)
    degree_to <- table(edges$to)
    all_degrees <- c(degree_from, degree_to)
    degree_sum <- tapply(all_degrees, names(all_degrees), sum)

    node_degrees <- dplyr::tibble(
      Protein = names(degree_sum),
      degree = as.numeric(degree_sum)
    )

    nodes <- dplyr::left_join(nodes, node_degrees, by = "Protein")
  }

  subnet_nodes <- unique(c(edges$from, edges$to))
  missing_nodes <- setdiff(names(terms), subnet_nodes)

  if (length(missing_nodes) > 0 && write) {
    missing_df <- df %>% dplyr::filter(.data$name %in% missing_nodes)
    readr::write_csv(
      missing_df,
      paste0(res.path, "/missing_nodes_", condition, spec_suffix(spec_cutoff), ".csv")
    )
  } else {
    missing_df <- NULL
  }

  kinograte_nodes <- nodes %>%
    as.data.frame() %>%
    dplyr::select(-dplyr::any_of("title"))

  if (write) {
    readr::write_csv(
      kinograte_nodes,
      paste0(res.path, "/nodes_", condition, spec_suffix(spec_cutoff), ".csv")
    )
    readr::write_csv(
      edges,
      paste0(res.path, "/edges_", condition, spec_suffix(spec_cutoff), ".csv")
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
    subnet_igraph <- igraph::graph_from_data_frame(
      d = edges[, c("from", "to", "weight")],
      directed = FALSE
    )

    wc <- suppressWarnings(
      igraph::edge.betweenness.community(
        subnet_igraph,
        weights = igraph::E(subnet_igraph)$weight,
        directed = FALSE,
        bridges = TRUE
      )
    )

    wc_df <- dplyr::tibble(id = wc$names, cluster = wc$membership)

    if (write) {
      readr::write_csv(
        wc_df,
        paste0(res.path, "/wc_df_", condition, spec_suffix(spec_cutoff), ".csv")
      )
    }

    return(list(
      network = subnet_igraph,
      nodes = nodes,
      edges = edges,
      missing = missing_df,
      wc_df = wc_df,
      maintitle = maintitle
    ))
  } else {
    subnet_igraph <- igraph::graph_from_data_frame(
      d = edges[, c("from", "to", "weight")],
      directed = FALSE
    )

    return(list(
      network = subnet_igraph,
      nodes = nodes,
      edges = edges,
      missing_nodes = missing_df,
      maintitle = maintitle
    ))
  }
}

