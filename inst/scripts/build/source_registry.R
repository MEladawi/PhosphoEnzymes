# helper_code/source_registry.R
# Single registry describing every input source: its local (version-less) filename,
# where to fetch the latest copy, whether it can be auto-fetched, and how to read the
# version string that ends up in the README. update_sources() refreshes the auto
# sources; build_source_manifest() records exactly what was used in this run.

# --- version readers --------------------------------------------------------

# Pull the "Release: 2026_01 of 28-Jan-2026" line out of the pkinfam header.
read_pkinfam_release <- function(file_path) {
  header_lines <- readLines(file_path, n = 40, warn = FALSE)
  release_line <- str_subset(header_lines, "^Release:")
  if (length(release_line)) str_squish(str_remove(release_line[1], "^Release:")) else NA_character_
}

# Files without an internal version string are tracked by their file date (local time, so
# the day boundary matches Sys.Date()).
file_modification_date <- function(file_path) as.character(as.Date(file.mtime(file_path), tz = Sys.timezone()))

# --- the UniProt KW-0418 REST query (reviewed human proteins, keyword "Kinase") ----
uniprot_keyword_kinase_url <- paste0(
  "https://rest.uniprot.org/uniprotkb/stream?",
  "query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0418%29",
  "&format=tsv&fields=accession,gene_primary,protein_name,ec,protein_families")

# --- GO molecular-function GMT variants (Bader Lab EM_Genesets) ---------------
# Two distributions of the same GO MF gene sets differ only in whether IEA (Inferred from
# Electronic Annotation) evidence is included. The IEA-inclusive file is the default: the
# functional gate is a RECALL-over-precision decision (we would rather include a real kinase
# annotated only electronically than miss it), and GO-umbrella-only singletons are already
# flagged low-confidence / non-curated downstream. Switch with build_kinase_list(go_include_iea=).
# Distinct local filenames so both variants can be cached side by side for offline reuse.
go_geneset_base <- "https://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/GO/"
GO_GENESET_VARIANTS <- list(
  with_iea = list(
    filename = "go_mf_genesets_with_iea_ensembl.gmt",
    url      = paste0(go_geneset_base, "Human_GO_mf_with_GO_iea_ensembl.gmt"),
    label    = "IEA-inclusive (electronic annotations included)"),
  no_iea = list(
    filename = "go_mf_genesets_no_iea_ensembl.gmt",
    url      = paste0(go_geneset_base, "Human_GO_mf_no_GO_iea_ensembl.gmt"),
    label    = "no IEA (manual/experimental annotations only)"))
# read_version factory for whichever GO variant is in use (one format, no duplication).
go_read_version <- function(variant_label) {
  function(file_path) paste0("Bader Lab EM_Genesets current_release, ", variant_label,
                             ", fetched ", file_modification_date(file_path))
}

# --- HGNC fetch (optionally pinned to an exact monthly archive) ---------------
# HGNC's plain download URL is not stable (the archive filename changes each release), so the
# hgnc package resolves an archive URL for us. A pinned URL (build_kinase_list(hgnc_archive_url=))
# fetches that exact release; otherwise the latest monthly is used. Either way the resolved
# archive URL is written to a sidecar so build_source_manifest() can record the EXACT input
# that produced this build (re-fetchable later), even on offline reruns.
HGNC_SOURCE_URL_SIDECAR <- "hgnc_source_url.txt"
make_hgnc_fetch <- function(pinned_archive_url = NULL) {
  function(destination) {
    walk(c("lubridate", "hgnc"),
         \(pkg) if (!requireNamespace(pkg, quietly = TRUE))
           install.packages(pkg, repos = "https://cloud.r-project.org"))
    resolved_url <- pinned_archive_url %||% hgnc::latest_monthly_url()
    hgnc::download_hgnc_dataset(url = resolved_url,
                                path = dirname(destination), filename = basename(destination))
    writeLines(resolved_url, file.path(dirname(destination), HGNC_SOURCE_URL_SIDECAR))
  }
}
read_hgnc_version <- function(file_path) {
  url_file <- file.path(dirname(file_path), HGNC_SOURCE_URL_SIDECAR)
  archive  <- if (file.exists(url_file)) readLines(url_file, n = 1, warn = FALSE) else "archive URL not recorded"
  paste0("HGNC monthly | archive: ", archive, " | file dated ", file_modification_date(file_path))
}

# --- the registry -----------------------------------------------------------
SOURCE_REGISTRY <- list(
  hgnc_complete_set = list(
    description    = "HGNC complete set (identifier bridge)",
    local_filename = "hgnc_complete_set.txt",
    download_url   = "https://www.genenames.org (monthly archive via the hgnc package)",
    auto_updatable = TRUE,
    # Default: latest monthly. configure_registry() can swap in a pinned-archive fetch.
    fetch          = make_hgnc_fetch(NULL),
    read_version   = read_hgnc_version),

  uniprot_pkinfam = list(
    description    = "UniProt pkinfam (curated protein kinome)",
    local_filename = "pkinfam.txt",
    download_url   = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pkinfam.txt",
    auto_updatable = TRUE,
    read_version   = read_pkinfam_release),

  # GO variant (filename/url/read_version) is finalised by configure_registry(); these are
  # the IEA-inclusive defaults.
  go_mf_genesets = list(
    description    = "GO molecular-function gene sets (Bader Lab EM_Genesets, Ensembl-keyed)",
    local_filename = GO_GENESET_VARIANTS$with_iea$filename,
    download_url   = GO_GENESET_VARIANTS$with_iea$url,
    auto_updatable = TRUE,
    read_version   = go_read_version(GO_GENESET_VARIANTS$with_iea$label)),

  manning_kinome = list(
    description    = "Manning kinome (kinase.com)",
    local_filename = "kinase.com_manning_list.xls",
    download_url   = "https://raw.githubusercontent.com/IDG-Kinase/DarkKinaseTools/master/data-raw/dark_kinases/kinase.com_list.xls",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("Manning et al. 2002 (static), fetched", file_modification_date(file_path))),

  # KinHub republishes the Manning kinome (Eid et al. 2017) and has not changed since. The
  # bundled input is a reconstructed facts table -- HGNC-normalised gene memberships and the
  # Manning group/family/subfamily labels, one row per gene -- not a copy of the source web
  # page; it is pinned (auto_updatable = FALSE) so offline reruns use exactly this file. It is
  # committed to the repo (see .gitignore exception) so it is available to fresh clones.
  kinhub = list(
    description    = "KinHub human kinase list (reconstructed facts)",
    local_filename = "kinhub_facts.tsv",
    download_url   = "http://www.kinhub.org/kinases.html",
    auto_updatable = FALSE,
    read_version   = function(file_path) paste("kinhub.org (Eid et al. 2017), reconstructed membership facts dated", file_modification_date(file_path))),

  uniprot_keyword_kinase = list(
    description    = "UniProtKB reviewed human, keyword KW-0418 (Kinase)",
    local_filename = "uniprot_kinase_KW-0418_human.tsv",
    download_url   = uniprot_keyword_kinase_url,
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("UniProtKB reviewed, fetched", file_modification_date(file_path))),

  idg_dark_kinome = list(
    description    = "IDG understudied (dark) kinome",
    local_filename = "IDG_dark_kinase_list.csv",
    download_url   = "https://raw.githubusercontent.com/IDG-Kinase/DarkKinaseTools/master/data-raw/dark_kinases/Dark%20Kinase%20List.csv",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("IDG DarkKinaseTools, fetched", file_modification_date(file_path))))

# --- runtime configuration of the registry ----------------------------------

# Finalise the registry for a run: select the GO MF GMT variant (IEA-inclusive by default;
# sets filename, URL, and a variant-recording read_version) and pin HGNC to an exact monthly
# archive URL (or keep latest-monthly when NULL). Single entry point so per-run knobs live in
# one place.
configure_registry <- function(registry, go_include_iea = TRUE, hgnc_archive_url = NULL) {
  variant <- if (isTRUE(go_include_iea)) GO_GENESET_VARIANTS$with_iea else GO_GENESET_VARIANTS$no_iea
  registry$go_mf_genesets$local_filename <- variant$filename
  registry$go_mf_genesets$download_url   <- variant$url
  registry$go_mf_genesets$read_version   <- go_read_version(variant$label)
  registry$hgnc_complete_set$fetch       <- make_hgnc_fetch(hgnc_archive_url)
  registry
}

# Convenience: absolute path of a registered source's local file.
source_file_path <- function(source_key, data_in_dir, registry = SOURCE_REGISTRY) {
  file.path(data_in_dir, registry[[source_key]]$local_filename)
}

# Refresh the auto-updatable sources. A file is re-downloaded only if it is missing or
# was last fetched on an earlier day (so same-day reruns reuse the copy on disk); set
# refresh = FALSE to force fully offline reuse, or force_refresh = TRUE to always fetch.
update_sources <- function(registry, data_in_dir, refresh = TRUE, force_refresh = FALSE) {
  for (source_key in names(registry)) {
    source_entry <- registry[[source_key]]
    file_path    <- file.path(data_in_dir, source_entry$local_filename)

    if (!isTRUE(source_entry$auto_updatable)) {
      if (!file.exists(file_path))
        warning(sprintf("Manual source missing (download it yourself): %s", source_entry$local_filename))
      else
        message(sprintf("  [manual] %s", source_entry$local_filename))
      next
    }

    already_fetched_today <- file.exists(file_path) &&
      as.Date(file.mtime(file_path), tz = Sys.timezone()) >= Sys.Date()
    if (!refresh || (already_fetched_today && !force_refresh)) {
      message(sprintf("  [cached] %s", source_entry$local_filename))
      next
    }

    message(sprintf("  [fetch ] %s", source_entry$local_filename))
    fetch_into <- if (is.function(source_entry$fetch)) {
      source_entry$fetch
    } else {
      function(destination) utils::download.file(source_entry$download_url, destination,
                                                 quiet = TRUE, mode = "wb")
    }
    # Fetch into a temporary file in the same folder and only replace the cached copy on
    # success, so a failed refresh can never destroy the existing file.
    temp_path <- tempfile(pattern = paste0(source_key, "_"), tmpdir = data_in_dir)
    fetched_ok <- tryCatch({
      fetch_into(temp_path)
      file.exists(temp_path) && file.info(temp_path)$size > 0L
    }, error = function(e) {
      message("           fetch failed: ", conditionMessage(e)); FALSE
    })

    if (fetched_ok) {
      file.rename(temp_path, file_path)                    # atomic replace (same filesystem)
    } else {
      unlink(temp_path)
      if (file.exists(file_path)) message("           keeping cached copy")
      else stop(sprintf("Could not fetch %s and no cached copy exists.", source_entry$local_filename))
    }
  }
}

# Record what was actually used this run (source, file, version, url, fetch date).
build_source_manifest <- function(registry, data_in_dir) {
  names(registry) |>
    map(\(source_key) {
      source_entry <- registry[[source_key]]
      file_path    <- file.path(data_in_dir, source_entry$local_filename)
      present      <- file.exists(file_path)
      tibble(
        source  = source_entry$description,
        file    = source_entry$local_filename,
        version = if (present) source_entry$read_version(file_path) %||% NA_character_ else "MISSING",
        fetched = if (present) file_modification_date(file_path) else NA_character_,
        url     = source_entry$download_url)
    }) |>
    list_rbind()
}
