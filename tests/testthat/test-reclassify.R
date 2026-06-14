# The term_sets= reclassification override re-types the catalog from raw per-gene evidence under a
# user-supplied EC/GO term set. Its correctness rests on the runtime gate twin agreeing with the
# build gate: reclassifying with the default term set must reproduce the shipped tables exactly.

# Assemble the four-table default term-set list. get_term_set() yields one table at a time, so the
# recompute path's named-list shape is built from the four calls (the same shape validate_term_set
# loads internally).
default_term_sets <- function() {
  list(
    kinase_ec      = PhosphoEnzymes::get_term_set("kinase", "ec"),
    kinase_go      = PhosphoEnzymes::get_term_set("kinase", "go"),
    phosphatase_ec = PhosphoEnzymes::get_term_set("phosphatase", "ec"),
    phosphatase_go = PhosphoEnzymes::get_term_set("phosphatase", "go"))
}

# Columns the override recomputes; these are exactly the typing/rigor/provenance fields that must be
# reproduced under the default term set.
.reclassify_round_trip_cols <- c(
  "acts_on_protein", "acts_on_nonprotein", "nonprotein_substrate_type", "substrate_call",
  "n_evidence_dimensions", "evidence_tier", "curated_core", "dual_protein_nonprotein")

test_that("get_kinases(term_sets = defaults) reproduces the shipped kinase master", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  base <- PhosphoEnzymes::get_kinases()
  recl <- PhosphoEnzymes::get_kinases(term_sets = default_term_sets())
  base <- base[order(base$ensembl_gene_id), , drop = FALSE]
  recl <- recl[order(recl$ensembl_gene_id), , drop = FALSE]

  expect_identical(base$ensembl_gene_id, recl$ensembl_gene_id)
  for (cc in .reclassify_round_trip_cols) {
    expect_equal(recl[[cc]], base[[cc]], info = cc)
  }
})

test_that("get_phosphatases(term_sets = defaults) reproduces the shipped phosphatase master", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  base <- PhosphoEnzymes::get_phosphatases()
  recl <- PhosphoEnzymes::get_phosphatases(term_sets = default_term_sets())
  base <- base[order(base$ensembl_gene_id), , drop = FALSE]
  recl <- recl[order(recl$ensembl_gene_id), , drop = FALSE]

  expect_identical(base$ensembl_gene_id, recl$ensembl_gene_id)
  for (cc in .reclassify_round_trip_cols) {
    expect_equal(recl[[cc]], base[[cc]], info = cc)
  }
})

test_that("dropping a non-protein lipid GO rule reduces lipid calls for that class only", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  ts  <- default_term_sets()
  kgo <- ts$kinase_go
  lipid_rows <- kgo$role == "rigor+substrate" & kgo$substrate == "nonprotein" &
    grepl("lipid", kgo$substrate_subtype)
  expect_gt(sum(lipid_rows), 0)
  ts$kinase_go <- kgo[!lipid_rows, , drop = FALSE]

  n_lipid_before <- sum(grepl("lipid", PhosphoEnzymes::get_kinases()$nonprotein_substrate_type))
  modified <- PhosphoEnzymes::get_kinases(term_sets = ts)
  n_lipid_after <- sum(grepl("lipid", modified$nonprotein_substrate_type))
  expect_lt(n_lipid_after, n_lipid_before)

  # The untouched class (phosphatases) is unaffected when only the kinase GO table changes.
  p_base <- PhosphoEnzymes::get_phosphatases()
  p_mod  <- PhosphoEnzymes::get_phosphatases(term_sets = ts)
  p_base <- p_base[order(p_base$ensembl_gene_id), , drop = FALSE]
  p_mod  <- p_mod[order(p_mod$ensembl_gene_id), , drop = FALSE]
  expect_equal(p_mod$nonprotein_substrate_type, p_base$nonprotein_substrate_type)
})

test_that("the default get_kinases() path runs and a term set with error rows warns, not stops", {
  # Default path: no term_sets=, must work regardless of the Suggests being reachable.
  expect_s3_class(PhosphoEnzymes::get_kinases(), "data.frame")

  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  # A missing citation is an error-severity validation issue; for a USER term set it degrades to a
  # warning and the reclassification proceeds (only the default set is held to a hard contract).
  bad_ec <- PhosphoEnzymes::get_term_set("kinase", "ec")
  bad_ec$citation[1] <- ""
  expect_warning(
    out <- PhosphoEnzymes::get_kinases(term_sets = list(kinase_ec = bad_ec)),
    "validation error")
  expect_s3_class(out, "data.frame")
})
