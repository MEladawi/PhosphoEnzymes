# Engine-function unit tests for the term-set machinery. The build engine is not in the
# package namespace, so locate it via system.file (installed / load_all) with a source-tree
# fallback, and attach the tidyverse namespaces the engine assumes are present.
suppressPackageStartupMessages({
  library(stringr); library(readr); library(purrr); library(tibble); library(dplyr)
})

build_file <- function(name) {
  p <- system.file("scripts", "build", name, package = "PhosphoEnzymes")
  if (nzchar(p) && file.exists(p)) return(p)
  alt <- file.path("..", "..", "inst", "scripts", "build", name)
  if (file.exists(alt)) return(alt)
  testthat::skip(paste("build engine not locatable:", name))
}

extdata_dir <- function() {
  p <- system.file("extdata", package = "PhosphoEnzymes")
  if (nzchar(p) && dir.exists(p)) return(p)
  alt <- file.path("..", "..", "inst", "extdata")
  if (dir.exists(alt)) return(alt)
  testthat::skip("inst/extdata not locatable")
}

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
