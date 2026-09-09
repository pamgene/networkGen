# Regression test for a real bug: library(networkGen) alone did not load
# PCSF's namespace (and so not its compiled DLL), because call_sr()'s
# .Call(..., PACKAGE = "PCSF") doesn't trigger lazy loading the way
# PCSF::something() does, and PCSF being listed in DESCRIPTION's Imports:
# isn't enough on its own -- R only eagerly loads an Imports package when
# this package's NAMESPACE has a real import()/importFrom() entry for it.
# Fixed by adding @importFrom PCSF construct_interactome near call_sr().
#
# This must run in a genuinely fresh subprocess (callr::r()), not just a
# fresh R6/environment within this session -- every previous test of this
# passed under pkgload::load_all(), which eagerly loads all Imports:
# regardless of NAMESPACE directives and so masked the bug completely. Only
# testing a real library()-after-install session catches this class of
# issue.
test_that("library(networkGen) alone loads PCSF's namespace (fresh subprocess)", {
  skip_if_not_installed("callr")
  skip_if_not(pcsf_functional(), "real PCSF compiled package not installed in this dev environment")

  result <- callr::r(function() {
    library(networkGen)
    list(
      pcsf_loaded = isNamespaceLoaded("PCSF"),
      symbol_available = tryCatch(
        {
          getNativeSymbolInfo("_PCSF_call_sr", PACKAGE = "PCSF")
          TRUE
        },
        error = function(e) FALSE
      )
    )
  })

  expect_true(result$pcsf_loaded)
  expect_true(result$symbol_available)
})
