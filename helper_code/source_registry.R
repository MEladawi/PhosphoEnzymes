# helper_code/source_registry.R
# Single registry describing every input source: its local (version-less) filename,
# where to fetch the latest copy, whether it can be auto-fetched, and how to read the
# version string that ends up in the README. update_sources() refreshes the auto
# sources; build_source_manifest() records exactly what was used in this run.

# --- version readers --------------------------------------------------------

# Pull the "Release: 2026_01 of 28-Jan-2026" line out of the pkinfam header.
read_pkinfam_release <- function(file_path) {
  header_lines <- readLines(file_path, n = 40, warn = FALSE)
  release_line <- grep("^Release:", header_lines, value = TRUE)
  if (length(release_line)) str_squish(sub("^Release:", "", release_line[1])) else NA_character_
}

# Files without an internal version string are tracked by their file date (local time, so
# the day boundary matches Sys.Date()).
file_modification_date <- function(file_path) as.character(as.Date(file.mtime(file_path), tz = Sys.timezone()))

# --- the UniProt KW-0418 REST query (reviewed human proteins, keyword "Kinase") ----
uniprot_keyword_kinase_url <- paste0(
  "https://rest.uniprot.org/uniprotkb/stream?",
  "query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0418%29",
  "&format=tsv&fields=accession,gene_primary,protein_name,ec,protein_families")

# --- the registry -----------------------------------------------------------
SOURCE_REGISTRY <- list(
  hgnc_complete_set = list(
    description    = "HGNC complete set (identifier bridge)",
    local_filename = "hgnc_complete_set.txt",
    download_url   = "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("HGNC monthly, file dated", file_modification_date(file_path))),

  uniprot_pkinfam = list(
    description    = "UniProt pkinfam (curated protein kinome)",
    local_filename = "pkinfam.txt",
    download_url   = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pkinfam.txt",
    auto_updatable = TRUE,
    read_version   = read_pkinfam_release),

  go_mf_genesets = list(
    description    = "GO molecular-function gene sets (Bader Lab EM_Genesets, Ensembl-keyed)",
    local_filename = "go_mf_genesets_ensembl.gmt",
    download_url   = "https://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/GO/Human_GO_mf_with_GO_iea_ensembl.gmt",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("Bader Lab EM_Genesets current_release, fetched", file_modification_date(file_path))),

  manning_kinome = list(
    description    = "Manning kinome (kinase.com)",
    local_filename = "kinase.com_manning_list.xls",
    download_url   = "https://raw.githubusercontent.com/IDG-Kinase/DarkKinaseTools/master/data-raw/dark_kinases/kinase.com_list.xls",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("Manning et al. 2002 (static), fetched", file_modification_date(file_path))),

  kinhub = list(
    description    = "KinHub human kinase list",
    local_filename = "kinhub_kinases.html",
    download_url   = "http://www.kinhub.org/kinases.html",
    auto_updatable = TRUE,
    read_version   = function(file_path) paste("kinhub.org (Eid et al. 2017, static), fetched", file_modification_date(file_path))),

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

# Convenience: absolute path of a registered source's local file.
source_file_path <- function(source_key, data_in_dir) {
  file.path(data_in_dir, SOURCE_REGISTRY[[source_key]]$local_filename)
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
    tryCatch(
      utils::download.file(source_entry$download_url, file_path, quiet = TRUE, mode = "wb"),
      error = function(e) {
        if (file.exists(file_path)) message("           download failed; keeping cached copy")
        else stop(sprintf("Could not download %s and no cached copy exists: %s",
                          source_entry$local_filename, conditionMessage(e)))
      })
  }
}

# Record what was actually used this run (source, file, version, url, fetch date).
build_source_manifest <- function(registry, data_in_dir) {
  map_dfr(names(registry), function(source_key) {
    source_entry <- registry[[source_key]]
    file_path    <- file.path(data_in_dir, source_entry$local_filename)
    present      <- file.exists(file_path)
    tibble(
      source  = source_entry$description,
      file    = source_entry$local_filename,
      version = if (present) source_entry$read_version(file_path) %||% NA_character_ else "MISSING",
      fetched = if (present) file_modification_date(file_path) else NA_character_,
      url     = if (is.na(source_entry$download_url)) "manual download (registration required)"
                else source_entry$download_url)
  })
}
