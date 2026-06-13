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
