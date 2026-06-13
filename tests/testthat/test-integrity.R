# Cross-class and cross-table integrity.
#
# These are pure identifier / derivation checks. The cross-class overlap is on the PROTEIN
# subsets: no human gene should be a protein kinase AND a protein phosphatase (that would be a
# parser bug). Non-protein bifunctional enzymes (e.g. PFKFB 6-phosphofructo-2-kinase/
# fructose-2,6-bisphosphatases) legitimately appear in both masters as non-protein members and
# must NOT trip this. A versioned allow-list (accepted_cross_class.csv, currently empty) covers a
# genuine bifunctional protein case.

test_that("no Ensembl gene ID is a protein enzyme in both masters", {
  k <- pe_kinases()
  p <- pe_phosphat()
  # Curated-core protein sets: a real bifunctional protein enzyme would be curated_core in both.
  # A comprehensive-only (Provisional, 0-axis) protein-phosphatase GO annotation on a real kinase
  # (e.g. an IBA-propagated term on LCK) is GO noise, not a pipeline collision, and is quarantined
  # in strict mode -- so the hard-fail keys on the curated core, not the comprehensive union.
  protein_in_both <- intersect(k$ensembl_gene_id[k$acts_on_protein & k$curated_core],
                               p$ensembl_gene_id[p$acts_on_protein & p$curated_core])
  allow_path <- system.file("extdata", "accepted_cross_class.csv", package = "PhosphoEnzymes")
  allow <- if (nzchar(allow_path)) utils::read.csv(allow_path)$ensembl_gene_id else character(0)
  unexpected <- setdiff(protein_in_both, allow)
  expect_identical(
    unexpected, character(0),
    info = paste("Genes that are a protein enzyme in BOTH masters (expected none):",
                 paste(unexpected, collapse = ", "))
  )
})

test_that("the unified summary is exactly the union of the two masters", {
  k  <- pe_kinases()
  p  <- pe_phosphat()
  pe <- pe_unified()

  expect_equal(nrow(pe), nrow(k) + nrow(p))
  expect_setequal(pe$ensembl_gene_id,
                  c(k$ensembl_gene_id, p$ensembl_gene_id))
  expect_setequal(
    pe$ensembl_gene_id[pe$regulator_class == "kinase"], k$ensembl_gene_id)
  expect_setequal(
    pe$ensembl_gene_id[pe$regulator_class == "phosphatase"], p$ensembl_gene_id)
})

test_that("shared columns in the summary match the masters row-for-row", {
  k  <- pe_kinases()
  p  <- pe_phosphat()
  pe <- pe_unified()
  shared <- c("acts_on_protein", "acts_on_nonprotein", "nonprotein_substrate_type",
              "n_evidence_dimensions", "evidence_tier", "curated_core")

  # Join class-routed: a non-protein bifunctional gene can be in both masters, so ensembl_gene_id
  # alone is ambiguous -- match on (ensembl_gene_id, regulator_class).
  master <- rbind(
    cbind(k[, c("ensembl_gene_id", shared)], regulator_class = "kinase"),
    cbind(p[, c("ensembl_gene_id", shared)], regulator_class = "phosphatase")
  )
  m <- merge(master, pe[, c("ensembl_gene_id", "regulator_class", shared)],
             by = c("ensembl_gene_id", "regulator_class"), suffixes = c(".master", ".pe"))
  for (col in shared) {
    expect_equal(m[[paste0(col, ".master")]], m[[paste0(col, ".pe")]],
                 info = paste("summary column drifted from master:", col))
  }
})

test_that("regulatory subunits do not leak into the catalytic phosphatase master", {
  p <- pe_phosphat()
  reg <- pe_try(PhosphoEnzymes::get_phosphatase_regulators)
  skip_if(is.null(reg) || nrow(reg) == 0L, "no regulator table")
  expect_identical(
    intersect(p$ensembl_gene_id, reg$ensembl_gene_id), character(0)
  )
})
