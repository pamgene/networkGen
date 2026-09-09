# Re-derive the equivalence proof behind removing the original PCSF path.
#
# Not part of the package: not in tests/testthat (so not run by
# devtools::test() / R CMD check), not sourced by anything in R/. Run it by
# hand only if there's real reason to re-question whether the fast path
# (kinograte_pg_pcsf() / PCSF_rand_pg(), still in R/pcsf-core.R) is safe to
# keep as the sole implementation -- e.g. after changing PCSF_rand_pg()'s
# random-noise logic, or upgrading the PCSF package.
#
#   Rscript archive/pcsf-equivalence/verify_equivalence.R
#
# What this checks: both paths solve PCSF via the same compiled entry point,
# call_sr(). This compares the raw `cost` vector each path builds
# (base_cost + base_cost * random_noise) and passes into call_sr(), per
# randomized run -- i.e. whether the two paths consume stats::runif()
# identically, not just whether their final output agrees. See
# docs/adr/0005-remove-original-pcsf-path.md for the decision this
# supports, and git history (pre-removal) for the version of this check that
# ran as a testthat test while both paths still existed in the package.

library(networkGen)

# Run from the package root: Rscript archive/pcsf-equivalence/verify_equivalence.R
original_fns_path <- if (file.exists("archive/pcsf-equivalence/original-pcsf-functions.R")) {
  "archive/pcsf-equivalence/original-pcsf-functions.R"
} else {
  "original-pcsf-functions.R"
}
source(original_fns_path)

capture <- new.env()
capture$fast <- list()
capture$orig <- list()
capture$which <- "fast"
trace(networkGen:::call_sr, tracer = quote({
  bucket <- if (capture$which == "fast") "fast" else "orig"
  capture[[bucket]][[length(capture[[bucket]]) + 1]] <- cost
}), print = FALSE)

run_check <- function(label, ppi_network, uka, b, n) {
  cat("\n===", label, "===\n")
  capture$fast <<- list()
  capture$orig <<- list()

  capture$which <<- "fast"
  fast <- networkGen:::kinograte_pg_pcsf(
    df = uka, ppi_network = ppi_network, spec_cutoff = 0, b = b,
    condition = "verify", write = FALSE, seed = 12345, n = n
  )
  capture$which <<- "orig"
  original <- kinograte_pg(
    df = uka, ppi_network = ppi_network, spec_cutoff = 0, b = b,
    condition = "verify", write = FALSE, seed = 12345, n = n
  )

  costs_ok <- length(capture$fast) == length(capture$orig) &&
    all(mapply(identical, capture$fast, capture$orig))
  nodes_ok <- identical(sort(fast$nodes$Protein), sort(original$nodes$Protein))
  edges_ok <- nrow(fast$edges) == nrow(original$edges)

  cat("  raw cost vectors bit-exact identical across", length(capture$fast), "runs:", costs_ok, "\n")
  cat("  node sets identical:", nodes_ok, "\n")
  cat("  edge counts identical:", edges_ok, "\n")
  all(costs_ok, nodes_ok, edges_ok)
}

toy_ppi <- data.frame(
  head = c("A", "B", "D", "B", "E", "F"),
  tail = c("B", "D", "C", "E", "F", "C"),
  cost = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
)
toy_uka <- data.frame(name = c("A", "B", "C"), prize = c(0.9, 0.9, 0.9), type = "Kinase", LogFC = c(1, 1, 1))

set.seed(42)
real_genes <- sample(unique(c(ppi_networkv12$head, ppi_networkv12$tail)), 15)
real_uka <- data.frame(name = real_genes, prize = stats::runif(15, 0.5, 1), type = "Kinase", LogFC = stats::rnorm(15))

ok1 <- run_check("toy topology (A-B-C-D-E-F, n=3)", toy_ppi, toy_uka, b = 10, n = 3)
ok2 <- run_check("real ppi_networkv12, 15 real genes (n=8)", ppi_networkv12, real_uka, b = 50, n = 8)

untrace(networkGen:::call_sr)

cat("\n=== RESULT:", if (ok1 && ok2) "PASS -- fast and original paths remain bit-exact equivalent" else "FAIL -- investigate before trusting the fast path alone", "===\n")
