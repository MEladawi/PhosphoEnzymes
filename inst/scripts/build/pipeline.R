# helper_code/pipeline.R
# The single entry-point function. build_kinase_list() runs the whole pipeline end to
# end -- refresh sources, map every source to Ensembl via HGNC, classify by function,
# and (optionally) write outputs -- and returns the result invisibly. All progress
# messaging lives here and in the helpers, keeping build_kinases.R a one-line call.

# Install any missing packages, then attach the ones the pipeline needs.
# `attached` are used directly (data wrangling, I/O, workbook). `installed_only` are used
# via `::` by the HGNC source fetch (hgnc / its dependency lubridate) and are deliberately
# NOT attached -- attaching lubridate would mask base intersect()/union()/setdiff().
ensure_packages <- function() {
  attached      <- c("readr","dplyr","tidyr","stringr","purrr","tibble","rvest","readxl","openxlsx")
  installed_only <- c("hgnc","lubridate")
  missing <- c(attached, installed_only)
  missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(invisible(lapply(attached, library, character.only = TRUE)))
}

# Single source of truth for the build version: the `version:` field of CITATION.cff. Keeping
# the tag here means a citation and the recorded build version cannot disagree.
read_pipeline_version <- function(citation_path = "CITATION.cff") {
  if (!file.exists(citation_path)) {
    "unversioned"
  } else {
    version_line <- str_subset(readLines(citation_path, warn = FALSE), "^version:")
    if (length(version_line)) str_squish(str_remove(version_line[1], "^version:")) else "unversioned"
  }
}

#' Build the comprehensive human kinase reference table.
#'
#' @param refresh_data  TRUE re-fetches auto-updatable sources (at most once/day);
#'                      FALSE reuses the cached files in `data_in_dir` (offline, reproducible).
#' @param data_in_dir   Folder holding (and caching) the source files.
#' @param output_dir    Folder for the generated outputs (created if needed).
#' @param write_files   TRUE writes the CSV/XLSX/list/manifest files; FALSE builds in memory only.
#' @param quiet         TRUE suppresses progress messages and the QC printout.
#' @param go_include_iea TRUE (default) uses the IEA-inclusive GO MF GMT; FALSE uses the
#'                      no-IEA (manual/experimental-only) variant. The choice is recorded in
#'                      the manifest.
#' @param hgnc_archive_url Optional exact HGNC monthly-archive URL to pin the identifier
#'                      bridge to a specific release; NULL (default) fetches the latest monthly
#'                      and records the resolved archive URL in the manifest.
#' @return (invisibly) a list with: `kinases` (the table), `unmapped`, `manifest`
#'         (source versions), and `sanity_passed` (TRUE if all QC sanity genes passed).
build_kinase_list <- function(refresh_data    = TRUE,
                              data_in_dir     = "data_in",
                              output_dir      = "data_out",
                              write_files     = TRUE,
                              quiet           = FALSE,
                              go_include_iea  = TRUE,
                              hgnc_archive_url = NULL) {

  run_pipeline <- function() {
    ensure_packages()
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    build_date       <- Sys.Date()
    pipeline_version <- read_pipeline_version()

    # Finalise the source registry for this run: GO IEA variant + optional HGNC pin.
    registry <- configure_registry(SOURCE_REGISTRY,
                                   go_include_iea   = go_include_iea,
                                   hgnc_archive_url = hgnc_archive_url)

    message("Checking and refreshing source files ...")
    update_sources(registry, data_in_dir, refresh = refresh_data)
    source_manifest <- build_source_manifest(registry, data_in_dir)
    path_for <- function(source_key) source_file_path(source_key, data_in_dir, registry)

    message("Reading HGNC complete set ...")
    hgnc_bridge <- build_hgnc_bridge(path_for("hgnc_complete_set"))

    message("Loading sources ...")
    go_sets           <- load_go_functional_sets(path_for("go_mf_genesets"))
    # Experimental-GO proxy for the evidence tier: membership in the NO-IEA GO kinase-activity
    # set (non-electronic annotations only). Loaded from the no-IEA GMT separately from whichever
    # GO variant feeds the gate, so the signal is available either way (file-level proxy).
    go_experimental_ids <- load_go_experimental_ids(
      file.path(data_in_dir, GO_GENESET_VARIANTS$no_iea$filename))
    pkinfam_source    <- load_pkinfam_kinome(path_for("uniprot_pkinfam"), hgnc_bridge)
    manning_source    <- load_manning_kinome(path_for("manning_kinome"), hgnc_bridge)
    kinhub_source     <- load_kinhub_kinome(path_for("kinhub"), hgnc_bridge)
    uniprot_kw_source <- load_uniprot_keyword_kinome(path_for("uniprot_keyword_kinase"), hgnc_bridge)
    idg_source        <- load_idg_dark_kinome(path_for("idg_dark_kinome"), hgnc_bridge)
    ec_source         <- load_ec_kinome(hgnc_bridge$gene_metadata)
    # Curated pseudokinase symbols -> base Ensembl IDs; sets catalytic_status in classify.
    pseudokinase_ensembl_ids <- read_csv(path_for("pseudokinases"), show_col_types = FALSE)$symbol |>
      map_chr(~ hgnc_bridge$resolve_to_ensembl(source_symbols = .x)) |>
      (\(ids) unique(ids[!is.na(ids)]))()

    membership <- list(
      pkinfam         = pkinfam_source$ensembl_ids,
      manning         = manning_source$ensembl_ids,
      kinhub          = kinhub_source$ensembl_ids,
      go_umbrella     = go_sets$kinase_activity_umbrella,
      ec              = ec_source$ensembl_ids,
      uniprot_keyword = uniprot_kw_source$ensembl_ids,
      idg_dark        = idg_source$ensembl_ids)
    universe_ensembl_ids <- reduce(membership, union) |>
      intersect(hgnc_bridge$gene_metadata$ensembl_gene_id) |> sort()
    message(sprintf("UNIVERSE: %d genes", length(universe_ensembl_ids)))

    kinase_taxonomy  <- build_kinase_taxonomy(uniprot_kw_source$taxonomy_table,
                                              kinhub_source$taxonomy_table,
                                              manning_source$taxonomy_table)
    unmapped_records <- bind_rows(pkinfam_source$unmapped, manning_source$unmapped,
                                  kinhub_source$unmapped, uniprot_kw_source$unmapped, idg_source$unmapped)

    message("Classifying and assembling ...")
    kinases_table <- classify_kinases(universe_ensembl_ids, hgnc_bridge, go_sets,
                                      ec_source, membership, kinase_taxonomy,
                                      go_experimental_ids = go_experimental_ids,
                                      pseudokinase_ensembl_ids = pseudokinase_ensembl_ids)

    if (write_files) {
      message("Writing outputs to ", output_dir, "/ ...")
      write_outputs(kinases_table, unmapped_records, output_dir, source_manifest,
                    build_date, pipeline_version)
    }
    sanity_passed <- qc_report(kinases_table, unmapped_records, verbose = !quiet,
                               go_protein_kinase_n = length(go_sets$protein_kinase_activity))

    invisible(list(kinases       = kinases_table,
                   unmapped      = unmapped_records,
                   manifest      = source_manifest,
                   sanity_passed = sanity_passed))
  }

  if (quiet) suppressMessages(run_pipeline()) else run_pipeline()
}
