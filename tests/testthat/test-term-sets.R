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
