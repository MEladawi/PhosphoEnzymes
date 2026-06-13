#' Access the human kinase reference table
#'
#' Two orthogonal knobs. `mode` filters rigor (`curated_core`) and `substrate`
#' filters the substrate booleans; the two never couple, so `mode = "strict"`
#' is not a protein filter -- a strictly-curated lipid or sugar kinase is kept
#' unless you also ask for `substrate = "protein"`.
#'
#' @param mode `"comprehensive"` (the full membership union, default) or
#'   `"strict"` (`curated_core == TRUE`: genes carrying at least one class-
#'   evidence kind -- an expert structural/sequence catalog and/or a class-
#'   specific EC).
#' @param substrate `"any"` (default), `"protein"` (keep `acts_on_protein`), or
#'   `"nonprotein"` (keep `acts_on_nonprotein`). A dual enzyme satisfies both.
#' @return A [tibble][tibble::tibble] of human kinases, one row per gene.
#' @examples
#' k <- get_kinases()
#' nrow(k)
#' get_kinases(mode = "strict", substrate = "protein")
#' @export
get_kinases <- function(mode = c("comprehensive", "strict"),
                        substrate = c("any", "protein", "nonprotein")) {
  mode <- match.arg(mode)
  substrate <- match.arg(substrate)
  .pe_filter(.pe_load("human_kinases"), mode, substrate)
}

#' Access the human protein-phosphatase reference table
#'
#' @inheritParams get_kinases
#' @return A [tibble][tibble::tibble] of human phosphatases, one row per gene.
#' @examples
#' p <- get_phosphatases()
#' nrow(p)
#' get_phosphatases(mode = "strict", substrate = "protein")
#' @export
get_phosphatases <- function(mode = c("comprehensive", "strict"),
                             substrate = c("any", "protein", "nonprotein")) {
  mode <- match.arg(mode)
  substrate <- match.arg(substrate)
  .pe_filter(.pe_load("human_phosphatases"), mode, substrate)
}

#' Access the unified human phospho-enzyme summary
#'
#' A thin, class-agnostic table spanning both kinases and phosphatases (one row per gene,
#' `regulator_class` distinguishing them), carrying only the shared substrate/evidence columns.
#' Join to [get_kinases()] / [get_phosphatases()] by `ensembl_gene_id` to recover class-specific
#' taxonomy. Ships all rows including Provisional; the summary has no `mode` knob (filter on
#' `evidence_tier` / `curated_core` yourself), but it carries `acts_on_protein` /
#' `acts_on_nonprotein`, so the `substrate` knob applies.
#'
#' @param substrate `"any"` (default), `"protein"`, or `"nonprotein"`.
#' @return A [tibble][tibble::tibble] spanning both enzyme classes.
#' @examples
#' pe <- get_phosphoenzymes()
#' table(pe$regulator_class)
#' get_phosphoenzymes(substrate = "protein")
#' @export
get_phosphoenzymes <- function(substrate = c("any", "protein", "nonprotein")) {
  substrate <- match.arg(substrate)
  .pe_filter(.pe_load("human_phosphoenzymes"), mode = "comprehensive", substrate = substrate)
}

#' Report the term-set provenance carried by a reference table
#'
#' The masters and the unified summary carry a `term_set_md5` attribute (the md5
#' fingerprints of the four default EC/GO term-set CSVs the build typed them
#' against). `provenance()` surfaces those fingerprints plus, when the `yaml`
#' package is available, the pinned source release strings from the shipped
#' build manifest.
#'
#' @param x A table returned by [get_kinases()], [get_phosphatases()], or
#'   [get_phosphoenzymes()].
#' @return A list with `term_set_md5` (the named md5 vector) and, when readable,
#'   `releases` (the source release strings from the build manifest).
#' @examples
#' provenance(get_kinases())
#' @export
provenance <- function(x) {
  md5 <- attr(x, "term_set_md5")
  out <- list(term_set_md5 = md5)
  manifest_path <- system.file("build_manifest.yaml", package = "PhosphoEnzymes")
  if (nzchar(manifest_path) && requireNamespace("yaml", quietly = TRUE)) {
    manifest <- yaml::read_yaml(manifest_path)
    if (!is.null(manifest$sources)) {
      out$releases <- stats::setNames(
        vapply(manifest$sources, function(s) as.character(s$version %||% NA_character_), character(1)),
        vapply(manifest$sources, function(s) as.character(s$source %||% NA_character_), character(1)))
    }
  }
  out
}

#' Read a default EC/GO term-set table
#'
#' Returns one of the four shipped declarative term-set CSVs (8 columns:
#' term_id, class, substrate, substrate_subtype, role, scope, citation, note).
#' These are the cited rules the build gate uses to type each gene by class and
#' substrate.
#'
#' @param class `"kinase"` or `"phosphatase"`.
#' @param type `"ec"` or `"go"`.
#' @return A [tibble][tibble::tibble] (or a data frame if tibble is unavailable),
#'   one row per term.
#' @examples
#' get_term_set("kinase", "ec")
#' @export
get_term_set <- function(class = c("kinase", "phosphatase"),
                         type = c("ec", "go")) {
  class <- match.arg(class)
  type <- match.arg(type)
  path <- system.file("extdata", paste0(class, "_", type, "_terms.csv"),
                      package = "PhosphoEnzymes")
  if (!nzchar(path)) stop("term-set CSV not found: ", class, "_", type, "_terms.csv")
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses = "character")
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

#' Write a custom term-set table to disk
#'
#' Structurally validates `table` (via [validate_term_set()]) and, if it carries
#' no error-severity issues, writes it as a CSV named `<class>_<type>_terms.csv`
#' into `dir`. The returned path can later feed a reclassification override.
#'
#' @param class `"kinase"` or `"phosphatase"`.
#' @param type `"ec"` or `"go"`.
#' @param table A data frame in the term-set schema.
#' @param dir Destination directory. Defaults to a per-session location under
#'   `tools::R_user_dir("PhosphoEnzymes", "data")`.
#' @return The written file path, invisibly.
#' @examples
#' \dontrun{
#' set_term_set("kinase", "ec", get_term_set("kinase", "ec"))
#' }
#' @export
set_term_set <- function(class = c("kinase", "phosphatase"),
                         type = c("ec", "go"), table, dir = NULL) {
  class <- match.arg(class)
  type <- match.arg(type)
  key <- paste0(class, "_", type)
  issues <- validate_term_set(stats::setNames(list(table), key))
  if (any(issues$severity == "error")) {
    msgs <- paste0("  - ", issues$message[issues$severity == "error"], collapse = "\n")
    stop("term-set validation failed:\n", msgs)
  }
  if (is.null(dir)) dir <- tools::R_user_dir("PhosphoEnzymes", "data")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(key, "_terms.csv"))
  utils::write.csv(table, path, row.names = FALSE, quote = TRUE)
  invisible(path)
}

#' Structurally validate one or more term-set tables
#'
#' Mirrors the build-side structural checks (a subset -- no gene-set resolution
#' or canary checks) so a user can lint a hand-edited term-set table before
#' feeding it back. With no argument it loads and validates the four shipped
#' default CSVs; otherwise pass a named list (any subset of
#' `kinase_ec` / `kinase_go` / `phosphatase_ec` / `phosphatase_go`).
#'
#' @param term_sets `NULL` (validate the four shipped defaults) or a named list
#'   of term-set data frames to validate.
#' @return A `data.frame(severity, table, term_id, message)` -- one row per
#'   issue, empty when clean. `severity == "error"` rows are fatal.
#' @examples
#' validate_term_set()
#' @export
validate_term_set <- function(term_sets = NULL) {
  if (is.null(term_sets)) {
    keys <- c("kinase_ec", "kinase_go", "phosphatase_ec", "phosphatase_go")
    term_sets <- stats::setNames(lapply(keys, function(k) {
      parts <- strsplit(k, "_", fixed = TRUE)[[1]]
      get_term_set(parts[1], parts[2])
    }), keys)
  }
  need <- c("term_id", "class", "substrate", "substrate_subtype",
            "role", "scope", "citation", "note")
  issues <- list()
  add <- function(severity, table, term_id, message)
    issues[[length(issues) + 1L]] <<- data.frame(
      severity = severity, table = table, term_id = term_id,
      message = message, stringsAsFactors = FALSE)

  for (nm in names(term_sets)) {
    tbl <- as.data.frame(term_sets[[nm]], stringsAsFactors = FALSE)
    miss <- setdiff(need, names(tbl))
    if (length(miss)) {
      add("error", nm, NA_character_,
          paste("missing columns:", paste(miss, collapse = ",")))
      next
    }
    blank <- function(v) is.na(v) | !nzchar(v)
    # Every row must carry a citation.
    no_cite <- tbl$term_id[blank(tbl$citation)]
    for (t in no_cite) add("error", nm, t, "missing citation")
    # Non-protein rigor rows must name a substrate_subtype.
    is_rig <- tbl$role == "rigor+substrate"
    bad_np <- tbl$term_id[is_rig & tbl$substrate == "nonprotein" &
                            blank(tbl$substrate_subtype)]
    for (t in bad_np) add("error", nm, t, "nonprotein rigor row lacks substrate_subtype")
    # Umbrella rows must be substrate-agnostic.
    bad_umb <- tbl$term_id[tbl$role == "rigor_umbrella" & tbl$substrate != "na"]
    for (t in bad_umb) add("error", nm, t, "umbrella row must have substrate == na")
    # No term_id may carry both protein and nonprotein among rigor+substrate rows.
    rig <- tbl[is_rig, , drop = FALSE]
    if (nrow(rig)) {
      n_sub <- vapply(split(rig$substrate, rig$term_id),
                      function(s) length(unique(s)), integer(1))
      for (t in names(n_sub)[n_sub > 1L])
        add("error", nm, t, "substrate overlap: term_id tagged both protein and nonprotein")
    }
  }
  if (length(issues)) Reduce(rbind, issues) else
    data.frame(severity = character(), table = character(),
               term_id = character(), message = character(),
               stringsAsFactors = FALSE)
}

# ---- internal helpers -------------------------------------------------------

#' @keywords internal
#' @noRd
.pe_filter <- function(df, mode = "comprehensive", substrate = "any") {
  keep_attr <- attr(df, "term_set_md5")
  if (mode == "strict") df <- df[!is.na(df$curated_core) & df$curated_core, , drop = FALSE]
  if (substrate == "protein")    df <- df[df$acts_on_protein, , drop = FALSE]
  if (substrate == "nonprotein") df <- df[df$acts_on_nonprotein, , drop = FALSE]
  if (is.null(attr(df, "term_set_md5")) && !is.null(keep_attr))
    attr(df, "term_set_md5") <- keep_attr
  df
}

#' @keywords internal
#' @noRd
.pe_load <- function(dataset) {
  e <- new.env(parent = emptyenv())
  utils::data(list = dataset, package = "PhosphoEnzymes", envir = e)
  obj <- get(dataset, envir = e)
  keep_attr <- attr(obj, "term_set_md5")
  out <- if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(obj) else obj
  if (!is.null(keep_attr)) attr(out, "term_set_md5") <- keep_attr
  out
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
