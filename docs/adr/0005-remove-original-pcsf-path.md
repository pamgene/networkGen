# Remove the original (igraph-based) PCSF path

`networkGen` carried two implementations of the same PCSF build:
`kinograte_pg_pcsf()`/`PCSF_rand_pg()` (dataframe-based, skips
`PCSF::construct_interactome()`) and `kinograte_pg()`/`PCSF_rand()`
(igraph-based, goes through it). The fast path was already the default
everywhere; the original was carried alongside it, unused in practice, with
only a weak automated test (matching node sets and edge counts, not exact
output) standing between "probably equivalent" and "definitely equivalent."

Verified before removing: intercepted `call_sr()`, the compiled solver
entry point both paths call, and compared the raw `cost` vector
(`base_cost + base_cost * random_noise`) each path builds per randomized
run -- i.e. whether the two paths consume `stats::runif()` identically, not
just whether their final output agrees. Bit-exact identical at two scales:
the package's toy A–B–C–D–E–F test topology (3 runs), and the real
`ppi_networkv12` (1.1M edges) with real gene data (8 runs, 569,232-length
cost vectors, character-for-character identical).

Decided: remove `kinograte_pg()`/`PCSF_rand()` and the `use_fast` parameter
from `generate_paired_network()`/`generate_kinase_network()` (both now
always use the fast path). The removed functions and the equivalence check
are not deleted outright -- they're kept as plain, non-package files under
`archive/pcsf-equivalence/` (excluded from the package build via
`.Rbuildignore`, not discovered by `testthat`), including a standalone
script that re-derives the bit-exact result on demand. This keeps the proof
a real, browsable file in the repo rather than something only recoverable
from git history, without keeping dead code in the package or a
now-untestable check in the routine test suite.

Why now: the equivalence test as it stood only checked derived output, not
the random values themselves -- not what "properly tested, with the seed and
all" requires. Strengthening it to bit-exact first, then removing the
original path only once that passed, keeps the removal justified by a
concrete result rather than an assumption carried over from when the fast
path was first written.
