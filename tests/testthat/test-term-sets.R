# Engine-function unit tests for the term-set machinery. Shared helpers (build_file,
# extdata_dir, tidyverse attaches) are provided by helper-build.R (testthat auto-sources
# helper-*.R before tests).

test_that("read_gmt_accession_map returns Ensembl-id vectors keyed by GO accession", {
  source(build_file("utils.R"), local = TRUE)
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c(
    "PROTEIN KINASE ACTIVITY%GOMF%GO:0004672\tna\tENSG1\tENSG2",
    "LIPID KINASE ACTIVITY%GOMF%GO:0001727\tna\tENSG3"), tmp)
  m <- read_gmt_accession_map(tmp)
  expect_setequal(names(m), c("GO:0004672", "GO:0001727"))
  expect_setequal(m[["GO:0004672"]], c("ENSG1", "ENSG2"))
  expect_identical(m[["GO:0001727"]], "ENSG3")
})

test_that("resolve_term_sets builds EC matchers and GO id sets per class", {
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts <- load_term_sets(extdata_dir())
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c(
    "PROTEIN KINASE ACTIVITY%GOMF%GO:0004672\tna\tENSGP",
    "LIPID KINASE ACTIVITY%GOMF%GO:0001727\tna\tENSGL",
    "PHOSPHOPROTEIN PHOSPHATASE ACTIVITY%GOMF%GO:0004721\tna\tENSGPP"), tmp)
  res <- resolve_term_sets(ts, tmp)
  expect_true(matches_ec_rules("2.7.10.1", res$kinase$ec_protein))
  expect_true(matches_ec_rules("2.7.1.67", res$kinase$ec_rigor))
  expect_false(matches_ec_rules("2.7.1.67", res$kinase$ec_protein))
  expect_true(matches_ec_rules("2.7.1.67", res$kinase$ec_nonprotein))
  expect_true(matches_ec_rules("3.1.3.16", res$phosphatase$ec_protein))
  expect_false(matches_ec_rules("3.1.3.67", res$phosphatase$ec_protein))
  expect_true(matches_ec_rules("3.9.1.3", res$phosphatase$ec_protein))
  expect_true(matches_ec_rules("3.6.1.1", res$phosphatase$ec_nonprotein))
  expect_true("ENSGP"  %in% res$kinase$go_protein_ids)
  expect_true("ENSGL"  %in% res$kinase$go_nonprotein_ids)
  expect_true("ENSGPP" %in% res$phosphatase$go_protein_ids)
})

test_that("validate_term_set passes the default set and catches seeded faults", {
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts <- load_term_sets(extdata_dir())
  issues <- validate_term_set(ts)
  expect_false(any(issues$severity == "error"))
  # seed an overlap fault: add a duplicate of a protein EC row tagged nonprotein
  bad <- ts
  ke  <- ts$tables$kinase_ec
  dup <- dplyr::mutate(dplyr::filter(ke, term_id == "2.7.10.-"), substrate = "nonprotein")
  bad$tables$kinase_ec <- dplyr::bind_rows(ke, dup)
  bad_issues <- validate_term_set(bad)
  expect_true(any(bad_issues$severity == "error" & grepl("overlap", bad_issues$message)))
})

test_that("term-set GO selection covers the frozen non-protein kinase panel", {
  gmt <- file.path(extdata_dir(), "go_mf_genesets_with_iea_ensembl.gmt")
  skip_if_not(file.exists(gmt), "GMT snapshot absent")
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts  <- load_term_sets(extdata_dir())
  res <- resolve_term_sets(ts, gmt)
  h <- readr::read_tsv(file.path(extdata_dir(), "hgnc_complete_set.txt"), show_col_types = FALSE)
  ensg <- function(sym) h$ensembl_gene_id[match(sym, h$symbol)]
  for (sym in c("PIP4K2A","PIP4K2B","PIP4K2C","PIKFYVE","PI4KA","SPHK1","DGKA")) {
    expect_true(ensg(sym) %in% res$kinase$go_nonprotein_ids,
                info = paste(sym, "should be in go_nonprotein_ids"))
  }
})

test_that("ec_axis_flags types protein vs nonprotein EC by the term set", {
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts  <- load_term_sets(extdata_dir())
  tmp <- tempfile(fileext = ".gmt"); writeLines("X%GOMF%GO:0004672\tna\tENSGP", tmp)
  res <- resolve_term_sets(ts, tmp)
  f_pi4ka <- ec_axis_flags("2.7.1.67", res$kinase)
  expect_true(f_pi4ka$ec_rigor); expect_false(f_pi4ka$ec_protein); expect_true(f_pi4ka$ec_nonprotein)
  expect_identical(f_pi4ka$nonprotein_subtypes, "lipid")
  f_phpt1 <- ec_axis_flags("3.9.1.3", res$phosphatase)
  expect_true(f_phpt1$ec_protein); expect_true(f_phpt1$ec_rigor)
  # a code in no curated set fires nothing
  f_none <- ec_axis_flags("9.9.9.9", res$kinase)
  expect_false(f_none$ec_rigor); expect_false(f_none$ec_protein); expect_false(f_none$ec_nonprotein)
})

test_that("GO:0052866 descendant set is pinned (broad-parent drift tripwire)", {
  # The phosphoinositide-phosphatase lipid subtype rides on this one broad subtree parent.
  # A broad parent is safe only while bounded: a future GO release could silently add a
  # non-lipid descendant to its propagated set without the term id changing. Pin the member
  # set; on drift, confirm every new member is a genuine phosphoinositide phosphatase, then
  # update this snapshot (and re-verify INPP5A is still excluded).
  skip_if_not(file.exists("../../inst/extdata/go_mf_genesets_with_iea_ensembl.gmt"),
              "GMT snapshot absent")
  source("../../inst/scripts/build/utils.R", local = TRUE)
  m <- read_gmt_accession_map("../../inst/extdata/go_mf_genesets_with_iea_ensembl.gmt")
  members <- sort(unique(m[["GO:0052866"]]))
  members <- members[nzchar(members)]
  pinned <- c(
    "ENSG00000003987", "ENSG00000040933", "ENSG00000063601", "ENSG00000078269",
    "ENSG00000087053", "ENSG00000100330", "ENSG00000102043", "ENSG00000108389",
    "ENSG00000109452", "ENSG00000110536", "ENSG00000112367", "ENSG00000115020",
    "ENSG00000122126", "ENSG00000132376", "ENSG00000132958", "ENSG00000139505",
    "ENSG00000148384", "ENSG00000155099", "ENSG00000159082", "ENSG00000163719",
    "ENSG00000165458", "ENSG00000171100", "ENSG00000185133", "ENSG00000198825",
    "ENSG00000204084", "ENSG00000211456", "ENSG00000274391", "ENSG00000281614",
    "ENSG00000284792", "ENSG00000291802")
  expect_setequal(members, pinned)
})

test_that("validate_term_set hard-errors on an obsolete GO term_id", {
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts <- load_term_sets(extdata_dir())
  # The shipped default set is clean (the obsolete GO:0004437 was removed).
  expect_false(any(validate_term_set(ts)$severity == "error" &
                     grepl("obsolete", validate_term_set(ts)$message)))
  # Re-introduce the known-obsolete id; the denylist must fail it.
  bad <- ts
  bad$tables$phosphatase_go <- dplyr::bind_rows(
    ts$tables$phosphatase_go,
    dplyr::tibble(term_id = "GO:0004437", class = "phosphatase", substrate = "nonprotein",
                  substrate_subtype = "lipid", role = "rigor+substrate", scope = "exact",
                  citation = "AmiGO", note = "seeded obsolete term"))
  issues <- validate_term_set(bad)
  expect_true(any(issues$severity == "error" & grepl("obsolete", issues$message)))
})
