#' Whether the real, compiled PCSF package (not the roxygen-doc-generation
#' stub used in this dev environment) is installed and usable.
pcsf_functional <- function() {
  tryCatch(
    {
      getNativeSymbolInfo("_PCSF_call_sr", PACKAGE = "PCSF")
      TRUE
    },
    error = function(e) FALSE
  )
}
