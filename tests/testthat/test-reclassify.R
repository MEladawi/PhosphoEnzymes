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

test_that("a non-default GO accession warns instead of a silent no-op", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  # The GMT-free recompute resolves GO membership from each gene's `go_terms` in
  # the sidecar, which records only the build's candidate accessions (the union
  # of the default GO term tables). A GO accession a user ADDS beyond those
  # defaults is absent from every gene's record and matches no gene -- even one
  # genuinely annotated with it in the full ontology -- so it would silently
  # under-call substrate. The accessor must instead warn, naming the bad id.
  # GO:0004683 = Ca2+/calmodulin-dependent kinase activity (not a default term).
  novel_go <- "GO:0004683"
  ts <- default_term_sets()
  expect_false(novel_go %in% ts$kinase_go$term_id)
  ts$kinase_go <- dplyr::bind_rows(
    ts$kinase_go,
    dplyr::tibble(
      term_id = novel_go, class = "kinase", substrate = "protein",
      substrate_subtype = "", role = "rigor+substrate", scope = "exact",
      citation = "AmiGO", note = "seeded non-default protein GO term"))

  expect_warning(out <- PhosphoEnzymes::get_kinases(term_sets = ts), novel_go)
  expect_s3_class(out, "data.frame")

  # The clean default set must NOT trip the guard. A benign vroom read warning
  # may still surface, so assert specifically that no guard warning is raised,
  # not that there is none at all.
  default_warnings <- testthat::capture_warnings(
    PhosphoEnzymes::get_kinases(term_sets = default_term_sets()))
  expect_false(any(grepl("silently inert", default_warnings)))
})

test_that("term_sets= aborts on a correctness fault but only warns on a hygiene fault", {
  # Default path: no term_sets=, must work regardless of the Suggests being reachable.
  expect_s3_class(PhosphoEnzymes::get_kinases(), "data.frame")

  skip_if_not_installed("dplyr")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")
  skip_if_not_installed("readr")

  # Hygiene fault: a missing citation does not change which genes a rule matches,
  # so for an inline USER term set it degrades to a warning and reclassification
  # proceeds (only a shipped or persisted set is held to the cited-data contract).
  bad_ec <- PhosphoEnzymes::get_term_set("kinase", "ec")
  bad_ec$citation[1] <- ""
  expect_warning(
    out <- PhosphoEnzymes::get_kinases(term_sets = list(kinase_ec = bad_ec)),
    "validation")
  expect_s3_class(out, "data.frame")

  # Correctness fault: a term_id tagged both protein and nonprotein would call a
  # gene "dual" off a contradiction, so the inline path must abort, not warn.
  overlap_ec <- PhosphoEnzymes::get_term_set("kinase", "ec")
  dup <- overlap_ec[overlap_ec$term_id == "2.7.10.-", , drop = FALSE]
  dup$substrate <- "nonprotein"
  overlap_ec <- rbind(overlap_ec, dup)
  expect_error(
    PhosphoEnzymes::get_kinases(term_sets = list(kinase_ec = overlap_ec)),
    "validation error")
})
