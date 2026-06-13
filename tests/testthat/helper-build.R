# Shared helpers for engine-function unit tests: locate the build engine (works under
# devtools::test() and an installed check) and attach the tidyverse namespaces the engine
# functions assume the caller has loaded.
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
