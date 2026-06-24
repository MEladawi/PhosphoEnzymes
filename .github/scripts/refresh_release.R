#!/usr/bin/env Rscript
# Scheduled maintenance: re-fetch the auto-updatable sources, rebuild the tables
# through the QC gate, and -- only if the table CONTENT actually changed -- bump
# the version, write a NEWS entry, and emit a PR summary. The calling workflow
# turns a real change into a pull request for human review.
#
# The change signal is a content fingerprint of the three masters AND the
# substrate-evidence sidecar, stable to serialization (and the masters stable to
# the term_set_md5 stamp), so a no-op rebuild does not look like a change. The
# sidecar is tracked because it is a shipped runtime input the masters cannot
# stand in for (it carries the raw GO/EC accessions the term_sets= recompute
# reads); the build manifest is deliberately excluded (volatile dates, and its
# only content-bearing fields are the per-table md5 stamps already covered).
# The QC gate inside make-data.R aborts the run on any sanity-gene or invariant
# failure, so a bad upstream release never reaches a PR.
#
# Two side findings are surfaced for the reviewer:
#   * which refreshable SOURCE FILES actually changed (by content md5, so the
#     fetch-date in the version string does not create noise, and pinned sources
#     -- which are never refreshed -- never appear);
#   * any GO term that went obsolete upstream, with the build's own message. An
#     obsolete term need not change membership (a redundant grouping term does
#     not), so when content is unchanged but a term went obsolete we still raise
#     an issue rather than bury the curation reminder in the run log.

stopifnot("run from the package root" = file.exists("DESCRIPTION"))

extdata <- "inst/extdata"
data_files <- file.path("data", c("human_kinases.rda", "human_phosphatases.rda",
                                  "human_phosphoenzymes.rda"))
sidecar_path <- file.path(extdata, "substrate_evidence.csv")
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Canonical, serialization-stable fingerprint of one shipped table.
canonical_md5 <- function(df) {
  attr(df, "term_set_md5") <- NULL
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(col)
    if (is.list(col)) vapply(col, paste, character(1), collapse = "|")
    else as.character(col))
  df <- df[order(df$ensembl_gene_id), sort(names(df)), drop = FALSE]
  tmp <- tempfile(); on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE, na = "")
  unname(tools::md5sum(tmp))
}

# Canonical, serialization-stable fingerprint of the sidecar CSV.
# The sidecar (substrate_evidence.csv) is a runtime input for the public
# term_sets= accessor override, which re-evaluates user-supplied term sets
# against raw GO/EC accessions stored in the sidecar. Unlike the build
# manifest (which carries volatile fetch/build dates and is covered by master
# fingerprints via its term_set_md5 stamp), the sidecar has no volatile fields
# and must be tracked directly: a refresh that swaps one GO accession for
# another in the same category changes the sidecar but leaves all master
# columns byte-identical.
canonical_md5_sidecar <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  # Natural key is (ensembl_gene_id, regulator_class); order by both.
  df <- df[do.call(order, df[c("ensembl_gene_id", "regulator_class")]),
           sort(names(df)), drop = FALSE]
  tmp <- tempfile(); on.exit(unlink(tmp))
  utils::write.csv(df, tmp, row.names = FALSE, na = "")
  unname(tools::md5sum(tmp))
}

load_one <- function(path) {
  e <- new.env(); load(path, envir = e); get(ls(e)[1], envir = e)
}

summarise_tables <- function(paths) {
  k <- load_one(paths[1]); p <- load_one(paths[2])
  list(
    fingerprint = vapply(paths, function(x) canonical_md5(load_one(x)), character(1)),
    sidecar_fingerprint = canonical_md5_sidecar(sidecar_path),
    n_kinases = nrow(k), n_phosphatases = nrow(p),
    protein_kinome = sum(k$acts_on_protein),
    protein_phosphatome = sum(p$acts_on_protein))
}

emit <- function(key, value) {
  out <- Sys.getenv("GITHUB_OUTPUT", "")
  line <- sprintf("%s=%s", key, value)
  if (nzchar(out)) cat(line, "\n", file = out, sep = "", append = TRUE)
  else message("[output] ", line)
}

# Content md5 of each REFRESHABLE source file (pinned sources are never listed).
source("inst/scripts/build/source_registry.R")
auto_sources <- Filter(function(x) isTRUE(x$auto_updatable), SOURCE_REGISTRY)
source_file_md5 <- function() {
  vapply(auto_sources, function(s) {
    p <- file.path(extdata, s$local_filename)
    if (file.exists(p)) unname(tools::md5sum(p)) else NA_character_
  }, character(1))
}
source_labels <- vapply(auto_sources, function(s) s$description, character(1))

# --- 1. Snapshot the committed state BEFORE the rebuild overwrites it ---------
before <- summarise_tables(data_files)
before_source_md5 <- source_file_md5()

# --- 2. Refresh-rebuild through the QC gate (stops the run on QC failure),
#        capturing any obsolete-GO warning the build raises. -------------------
Sys.setenv(PHOSPHOENZYMES_REFRESH = "true")
obsolete_msgs <- character(0)
withCallingHandlers(
  source("inst/scripts/make-data.R"),
  warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("obsolete", msg, ignore.case = TRUE))
      obsolete_msgs <<- c(obsolete_msgs, trimws(msg))
    # do not muffle -- let the warning still surface in the run log
  })
obsolete_msgs <- unique(obsolete_msgs)

# --- 3. Snapshot the rebuilt state -------------------------------------------
after <- summarise_tables(data_files)
after_source_md5 <- source_file_md5()

moved <- source_labels[which(before_source_md5 != after_source_md5 |
                             is.na(before_source_md5) != is.na(after_source_md5))]
source_lines <- if (length(moved)) {
  vapply(unname(moved), function(s) sprintf("- %s", s), character(1))
} else {
  "- (no refreshable source file changed content)"
}

obsolete_block <- if (length(obsolete_msgs)) {
  c("", "## :warning: Newly-obsolete GO terms (curation needed)",
    "An upstream GO release retired a term still cited in the term sets. Remove or",
    "replace it in the term-set CSV and add it to the obsolete-GO denylist:", "",
    paste0("> ", unlist(strsplit(obsolete_msgs, "\n"))))
} else {
  character(0)
}

# --- 4. Surface an obsolete term even when the content did not change ---------
emit("obsolete", if (length(obsolete_msgs)) "true" else "false")
if (length(obsolete_msgs))
  writeLines(c("A scheduled refresh found a GO term that went obsolete upstream.",
               "Membership did not necessarily change, so no release PR was opened,",
               "but the term set needs a curation edit:", obsolete_block),
             ".github/obsolete_body.md")

# --- 5. Did the table content change? ----------------------------------------
# Check BOTH master tables and sidecar: either masters differ OR sidecar differs.
masters_changed <- !identical(before$fingerprint, after$fingerprint)
sidecar_changed <- !identical(before$sidecar_fingerprint, after$sidecar_fingerprint)
if (!masters_changed && !sidecar_changed) {
  message("No content change after refresh -- nothing to release.")
  emit("changed", "false")
  quit(save = "no", status = 0)
}
message("Table content changed after refresh -- preparing a release.")

# --- 6. Bump the patch version -----------------------------------------------
desc <- readLines("DESCRIPTION")
vline <- grep("^Version:", desc)
old_version <- trimws(sub("^Version:", "", desc[vline]))
parts <- as.integer(strsplit(old_version, ".", fixed = TRUE)[[1]])
parts[length(parts)] <- parts[length(parts)] + 1L
new_version <- paste(parts, collapse = ".")
desc[vline] <- paste0("Version: ", new_version)
writeLines(desc, "DESCRIPTION")

# --- 7. Change summary (counts + which source files moved) -------------------
delta <- function(label, a, b)
  sprintf("- %s: %d -> %d (%+d)", label, a, b, b - a)
count_lines <- c(
  delta("kinases", before$n_kinases, after$n_kinases),
  delta("protein kinome", before$protein_kinome, after$protein_kinome),
  delta("phosphatases", before$n_phosphatases, after$n_phosphatases),
  delta("protein phosphatome", before$protein_phosphatome,
        after$protein_phosphatome))

# --- 8. NEWS entry -----------------------------------------------------------
news_entry <- c(
  sprintf("# PhosphoEnzymes %s", new_version), "",
  "* Data refresh from updated upstream sources (automated maintenance). The",
  "  rebuild passed the full QC sanity gate. Table changes:", "",
  count_lines, "")
news <- if (file.exists("NEWS.md")) readLines("NEWS.md") else character(0)
writeLines(c(news_entry, news), "NEWS.md")

# --- 9. PR body --------------------------------------------------------------
pr_body <- c(
  sprintf("Automated source refresh: version %s -> %s.", old_version, new_version),
  "", "The refreshable upstream sources were re-fetched and the tables rebuilt;",
  "the full QC sanity gate passed. This PR is for review before release.",
  "", "## Table changes", count_lines,
  "", "## Source files that changed", source_lines,
  obsolete_block)
writeLines(pr_body, ".github/pr_body.md")

emit("changed", "true")
emit("version", new_version)
message(sprintf("Prepared release %s.", new_version))
