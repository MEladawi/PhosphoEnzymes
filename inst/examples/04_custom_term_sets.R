# The EC/GO typing rules ship as cited data, and are overridable.
#
# get_term_set() reads one of the four rule tables; validate_term_set() lints
# them; and get_kinases(term_sets = ...) / get_phosphatases(term_sets = ...)
# recompute the substrate/rigor/provenance columns under your own rules WITHOUT
# rebuilding the package. Passing the defaults back reproduces the shipped table
# exactly -- the round-trip identity that anchors every custom recompute.
#
# The recompute path is the one place that needs the Suggests packages
# (dplyr/purrr/stringr/readr); this script skips that part if they are absent.

library(PhosphoEnzymes)

# Each row of a term set is a cited rule: a term id, the substrate it implies,
# and whether it contributes to rigor, to substrate typing, or both.
kinase_go <- get_term_set("kinase", "go")
cat("First rows of the kinase GO rule table:\n")
print(head(kinase_go[, c("term_id", "substrate", "substrate_subtype", "role")]))

cat("\nLinting the four shipped defaults (0 issues = clean):\n")
print(nrow(validate_term_set()))

suggests <- c("dplyr", "purrr", "stringr", "readr")
have_suggests <- all(vapply(suggests, requireNamespace, logical(1),
                            quietly = TRUE))

if (!have_suggests) {
  message("Skipping the term_sets= recompute: needs ",
          paste(suggests, collapse = ", "))
} else {
  # Round-trip: feeding the default term set back reproduces substrate_call.
  default_ts <- list(
    kinase_ec = get_term_set("kinase", "ec"),
    kinase_go = get_term_set("kinase", "go")
  )
  recomputed <- get_kinases(term_sets = default_ts)
  cat("\nRound-trip identity on substrate_call:\n")
  print(identical(
    as.data.frame(recomputed)[, "substrate_call"],
    as.data.frame(get_kinases())[, "substrate_call"]
  ))

  # A custom edit: drop the lipid-kinase GO rules. Those genes stop being typed
  # as lipid kinases; their rigor is untouched (the two axes are independent).
  ts2 <- get_term_set("kinase", "go")
  keep <- ts2$substrate_subtype != "lipid" | is.na(ts2$substrate_subtype)
  retyped <- get_kinases(term_sets = list(kinase_go = ts2[keep, ]))
  cat("\nLipid-typed kinases, shipped vs lipid-GO-dropped:\n")
  print(c(
    shipped = sum(grepl("lipid", get_kinases()$nonprotein_substrate_type)),
    retyped = sum(grepl("lipid", retyped$nonprotein_substrate_type))
  ))
}
