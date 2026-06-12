# Helpers shared across test files.
#
# The built tables ship as package data (data/*.rda) produced by
# inst/scripts/make-data.R. Before that pipeline has been run (e.g. on a fresh
# skeleton), the data objects do not exist yet. These helpers load a table and
# skip() the test if it is not available, so the package checks cleanly at every
# stage of development rather than erroring on missing data.

pe_try <- function(fn, ...) {
  out <- tryCatch(fn(...), error = function(e) e)
  if (inherits(out, "error")) {
    testthat::skip(paste0("dataset not built yet: ", conditionMessage(out)))
  }
  out
}

pe_kinases  <- function(...) pe_try(PhosphoEnzymes::get_kinases, ...)
pe_phosphat <- function(...) pe_try(PhosphoEnzymes::get_phosphatases, ...)
pe_unified  <- function()    pe_try(PhosphoEnzymes::get_phosphoenzymes)

# Look up a single gene by symbol; returns a one-row data frame or NULL.
pe_row <- function(df, sym) {
  hit <- df[df$symbol == sym, , drop = FALSE]
  if (nrow(hit) == 0L) return(NULL)
  hit
}

# Myotubularin family symbols (all lipid phosphatases; several are
# pseudophosphatases). Used by the category invariants.
PE_MTM_FAMILY <- c(
  "MTM1",
  paste0("MTMR", 1:14),
  "SBF1", "SBF2"   # MTMR5 / MTMR13 aliases
)

# --- Kinase family symbols (used by the kinase-side category invariants) ------
# Each family is unambiguously non-protein end-to-end, so a family invariant can
# assert no member acts on protein -- EXCEPT the diacylglycerol kinases, which
# include DGKQ, a genuine dual lipid+protein kinase. The DGK invariant therefore
# allows the dual exception; the others are strict. (PI3K is deliberately not a
# family here: class I p110s are dual, so PIK3CA is guarded as a trap case.)
PE_DGK_FAMILY  <- paste0("DGK", c("A", "B", "D", "E", "G", "H", "I", "K", "Q", "Z"))
PE_HK_FAMILY   <- c("HK1", "HK2", "HK3", "GCK", "HKDC1")
PE_AK_FAMILY   <- paste0("AK", 1:9)
PE_PI4K_FAMILY <- c("PI4KA", "PI4KB", "PI4K2A", "PI4K2B")
