# Maintainer-only regeneration of the shipped reference tables.
#
# This is NOT run at install or load time. It rebuilds data/*.rda from the pinned
# source snapshots in inst/extdata/ with no live network access, so the shipped
# tables are deterministic and citable. Run it from the package root on a release:
#
#   Rscript inst/scripts/make-data.R
#
# then commit the regenerated data/*.rda and inst/build_manifest.yaml.

stopifnot("run from the package root (DESCRIPTION not found)" = file.exists("DESCRIPTION"))

# The build engine lives outside the installed namespace; source it here.
build_dir <- file.path("inst", "scripts", "build")
invisible(lapply(list.files(build_dir, pattern = "[.]R$", full.names = TRUE), source))

# --- Kinases -----------------------------------------------------------------
# Build offline from the committed snapshots; assert the functional-gate QC gate.
result <- build_kinase_list(
  refresh_data = FALSE,
  data_in_dir  = "inst/extdata",
  output_dir   = tempfile("kinase_build_"),   # legacy file outputs not needed here
  write_files  = FALSE,
  quiet        = FALSE)

if (!isTRUE(result$sanity_passed))
  stop("QC sanity genes failed; aborting data build.")

human_kinases <- harmonize_kinases_to_package_schema(result$kinases)
message(sprintf("human_kinases: %d rows x %d cols", nrow(human_kinases), ncol(human_kinases)))

if (!requireNamespace("usethis", quietly = TRUE))
  stop("Package 'usethis' is required to write data/*.rda.")
usethis::use_data(human_kinases, compress = "xz", overwrite = TRUE)

# --- Build manifest ----------------------------------------------------------
# One machine-readable provenance file: build date, package version, per-table
# row count + content hash, and the recorded version of every input snapshot.
package_version <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
rda_path <- file.path("data", "human_kinases.rda")
src <- result$manifest

manifest <- list(
  build_date      = as.character(Sys.Date()),
  package_version = unname(package_version),
  tables = list(
    human_kinases = list(
      n_rows = nrow(human_kinases),
      n_cols = ncol(human_kinases),
      md5    = unname(tools::md5sum(rda_path)))),
  sources = lapply(seq_len(nrow(src)), function(i) list(
    source  = src$source[i],
    file    = src$file[i],
    version = src$version[i],
    fetched = as.character(src$fetched[i]))))

if (requireNamespace("yaml", quietly = TRUE)) {
  yaml::write_yaml(manifest, file.path("inst", "build_manifest.yaml"))
  message("Wrote inst/build_manifest.yaml")
} else {
  warning("Package 'yaml' not available; skipped inst/build_manifest.yaml")
}

message("make-data.R complete.")
