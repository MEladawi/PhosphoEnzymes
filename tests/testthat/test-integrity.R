# Cross-class and cross-table integrity.
#
# These are pure identifier / derivation checks - uncontroversial, no biology to
# argue with. The cross-class overlap is the CI-hard form of the build-time
# warning described in the plan (no whitelist; expected empty, since no known
# human gene is both a protein kinase and a protein phosphatase).

test_that("no Ensembl gene ID appears in both the kinase and phosphatase masters", {
  k <- pe_kinases()
  p <- pe_phosphat()
  overlap <- intersect(k$ensembl_gene_id, p$ensembl_gene_id)
  expect_identical(
    overlap, character(0),
    info = paste("Cross-class overlap (expected none):",
                 paste(overlap, collapse = ", "))
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

  master <- rbind(
    k[, c("ensembl_gene_id", shared)],
    p[, c("ensembl_gene_id", shared)]
  )
  m <- merge(master, pe[, c("ensembl_gene_id", shared)],
             by = "ensembl_gene_id", suffixes = c(".master", ".pe"))
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
