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

# --- Phosphatases ------------------------------------------------------------
phos <- build_phosphatase_list(refresh_data = FALSE, data_in_dir = "inst/extdata", quiet = FALSE)
human_phosphatases <- harmonize_phosphatases_to_package_schema(phos$phosphatases)
message(sprintf("human_phosphatases: %d rows x %d cols (protein subset %d)",
                nrow(human_phosphatases), ncol(human_phosphatases), sum(human_phosphatases$acts_on_protein)))

# --- Unified summary + provenance (derived from the two masters) -------------
human_phosphoenzymes  <- build_unified_summary(human_kinases, human_phosphatases)
membership_provenance <- build_membership_provenance(human_kinases, human_phosphatases)
message(sprintf("human_phosphoenzymes: %d rows; membership_provenance: %d rows",
                nrow(human_phosphoenzymes), nrow(membership_provenance)))

if (!requireNamespace("usethis", quietly = TRUE))
  stop("Package 'usethis' is required to write data/*.rda.")
usethis::use_data(human_kinases, human_phosphatases, human_phosphoenzymes,
                  compress = "xz", overwrite = TRUE)

# Static CSVs for non-R (Python/CLI/HPC) consumers + the provenance sidecar, shipped in extdata.
readr::write_csv(human_phosphoenzymes,  file.path("inst", "extdata", "human_phosphoenzymes.csv"))
readr::write_csv(membership_provenance, file.path("inst", "extdata", "membership_provenance.csv"))

# --- Build manifest ----------------------------------------------------------
# One machine-readable provenance file: build date, package version, per-table
# row count + content hash, and the recorded version of every input snapshot.
package_version <- read.dcf("DESCRIPTION", fields = "Version")[1, 1]
src <- result$manifest
table_entry <- function(name, df) list(
  n_rows = nrow(df), n_cols = ncol(df),
  md5    = unname(tools::md5sum(file.path("data", paste0(name, ".rda")))))

manifest <- list(
  build_date      = as.character(Sys.Date()),
  package_version = unname(package_version),
  tables = list(
    human_kinases        = table_entry("human_kinases", human_kinases),
    human_phosphatases   = table_entry("human_phosphatases", human_phosphatases),
    human_phosphoenzymes = table_entry("human_phosphoenzymes", human_phosphoenzymes)),
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
