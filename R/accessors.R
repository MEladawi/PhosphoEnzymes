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
#' @param term_sets `NULL` (default) returns the shipped master, typed against the
#'   default EC/GO term sets. Otherwise a named list of EC/GO term-set tables
#'   (any subset of `kinase_ec` / `kinase_go` / `phosphatase_ec` /
#'   `phosphatase_go`; missing tables fall back to the shipped defaults) under
#'   which the substrate / rigor / provenance columns are recomputed from raw
#'   per-gene evidence before the `mode` / `substrate` filters apply -- so the
#'   catalog can be re-typed under a customized term set without rebuilding.
#'   Supplying a term set with validation errors warns (it does not stop), unlike
#'   the default set which is held to a stricter contract. This path requires the
#'   `dplyr`, `purrr`, and `stringr` packages; the default path does not.
#' @return A [tibble][tibble::tibble] of human kinases, one row per gene.
#' @examples
#' k <- get_kinases()
#' nrow(k)
#' get_kinases(mode = "strict", substrate = "protein")
#' @export
get_kinases <- function(mode = c("comprehensive", "strict"),
                        substrate = c("any", "protein", "nonprotein"),
                        term_sets = NULL) {
  mode <- match.arg(mode)
  substrate <- match.arg(substrate)
  if (is.null(term_sets)) {
    .pe_filter(.pe_load("human_kinases"), mode, substrate)
  } else {
    .pe_reclassify("human_kinases", "kinase", term_sets, mode, substrate)
  }
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
                             substrate = c("any", "protein", "nonprotein"),
                             term_sets = NULL) {
  mode <- match.arg(mode)
  substrate <- match.arg(substrate)
  if (is.null(term_sets)) {
    .pe_filter(.pe_load("human_phosphatases"), mode, substrate)
  } else {
    .pe_reclassify("human_phosphatases", "phosphatase", term_sets, mode, substrate)
  }
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
#' The `term_sets=` reclassification override is class-specific (it recomputes one class's columns
#' from that class's per-gene evidence) and so is not offered here. To re-type under a custom term
#' set, call [get_kinases()] / [get_phosphatases()] with `term_sets=` and join the results.
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
#' \donttest{
#' set_term_set("kinase", "ec", get_term_set("kinase", "ec"), dir = tempdir())
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
    # GO term_ids verified obsolete in the pinned ontology (QuickGO isObsolete == true). A retired
    # term resolves to an empty gene set, so it must never sit in a term set: hard error here,
    # distinct from a valid-but-unannotated term. Extend only with confirmed-obsolete IDs.
    obsolete_go <- c("GO:0004437")
    for (t in intersect(tbl$term_id, obsolete_go))
      add("error", nm, t, "obsolete GO term_id (retired in the ontology; remove or replace with its successor)")
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
# Normalize a user `term_sets=` argument into the four-table named list the recompute path needs.
# Accepts either the full four-table list (kinase_ec / kinase_go / phosphatase_ec / phosphatase_go)
# or the value of get_term_set() (a single table); a single table is insufficient on its own, so any
# missing default table is filled from the shipped defaults. Returns a named list of data frames.
.pe_normalize_term_sets <- function(term_sets) {
  keys <- c("kinase_ec", "kinase_go", "phosphatase_ec", "phosphatase_go")
  defaults <- stats::setNames(lapply(keys, function(k) {
    parts <- strsplit(k, "_", fixed = TRUE)[[1]]
    get_term_set(parts[1], parts[2])
  }), keys)
  if (is.null(names(term_sets)) || !all(names(term_sets) %in% keys)) {
    stop("term_sets must be a named list with names among: ",
         paste(keys, collapse = ", "),
         " (e.g. list(kinase_ec = ..., kinase_go = ..., phosphatase_ec = ..., phosphatase_go = ...))")
  }
  # Replace whole tables (never merge column-wise): a supplied table wins, the rest stay default.
  for (k in intersect(names(term_sets), keys)) defaults[[k]] <- term_sets[[k]]
  defaults
}

#' @keywords internal
#' @noRd
# Re-type the shipped master for one class under a user-supplied term set, then
# apply mode/substrate. This is the ONLY path that touches the dplyr/purrr/
# stringr/readr Suggests; it guards on them up front. The default (no term_sets=)
# accessors stay dependency-free base R. The recompute mirrors the build gate
# exactly via the runtime twin in R/term-sets-runtime.R, so passing the default
# term set reproduces the shipped table.
.pe_reclassify <- function(dataset, class, term_sets, mode, substrate) {
  needed <- c("dplyr", "purrr", "stringr", "readr")
  missing_pkgs <- needed[!vapply(
    needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("term_sets= reclassification requires the dplyr, purrr, stringr ",
         "and readr packages; install: ",
         paste(missing_pkgs, collapse = ", "))
  }
  tables <- .pe_normalize_term_sets(term_sets)

  # User term sets degrade to warnings (only the default set hard-errors); summarize any error rows.
  issues <- validate_term_set(tables)
  err <- issues[issues$severity == "error", , drop = FALSE]
  if (nrow(err)) {
    warning("supplied term_sets has ", nrow(err), " validation error(s); proceeding anyway:\n",
            paste0("  - [", err$table, "] ", ifelse(is.na(err$term_id), "", paste0(err$term_id, ": ")),
                   err$message, collapse = "\n"), call. = FALSE)
  }

  # The Ensembl-keyed GMT is not shipped in the installed package, so GO membership is resolved at
  # the accession level: the per-gene `go_terms` in the sidecar (membership among the candidate GO
  # accessions in the build-time GMT) is intersected with the term set's protein / non-protein /
  # per-subtype GO accession lists. EC typing needs no gene-set file at all. The result is identical
  # to the build's id-set membership for any term set drawn from the shipped candidate accessions.
  go_tbl <- tibble::as_tibble(tables[[paste0(class, "_go")]])
  ec_tbl <- tibble::as_tibble(tables[[paste0(class, "_ec")]])
  go_acc <- .pe_resolve_go_accessions_class(go_tbl)
  resolved_class <- .pe_resolve_ec_class(ec_tbl)

  sidecar_path <- system.file("extdata", "substrate_evidence.csv", package = "PhosphoEnzymes")
  if (!nzchar(sidecar_path)) stop("substrate_evidence sidecar not found in extdata")
  sidecar <- readr::read_csv(sidecar_path, show_col_types = FALSE) |>
    dplyr::filter(regulator_class == class)

  # Build the per-gene evidence tibble exactly as classify.R / classify_phosphatases.R do: EC codes
  # from splitting the sidecar's pipe field; the four co-equal GO/Chen substrate signals from the
  # gene's annotated GO accessions intersected with the resolved accession lists.
  evidence <- sidecar |>
    dplyr::mutate(
      all_ec_codes = purrr::map(ec_codes, .pe_split_pipe_delimited),
      gene_go      = purrr::map(go_terms, .pe_split_pipe_delimited),
      go_protein    = purrr::map_lgl(gene_go, \(g) any(g %in% go_acc$go_protein_accessions)),
      go_nonprotein = purrr::map_lgl(gene_go, \(g) any(g %in% go_acc$go_nonprotein_accessions)),
      go_nonprotein_subtype = purrr::map_chr(gene_go, \(g) {
        hits <- names(go_acc$go_nonprotein_subtype_accessions)[
          purrr::map_lgl(go_acc$go_nonprotein_subtype_accessions, \(acc) any(g %in% acc))]
        paste(hits, collapse = "|")
      }),
      chen_nonprotein       = as.logical(chen_nonprotein),
      in_structural_catalog = as.logical(in_structural_catalog),
      supplementary_support = as.logical(supplementary_support),
      go_experimental_protein = as.logical(go_experimental_protein),
      catalytic_status      = as.character(catalytic_status)) |>
    dplyr::select(ensembl_gene_id, all_ec_codes, go_protein, go_nonprotein, go_nonprotein_subtype,
                  chen_nonprotein, in_structural_catalog, supplementary_support,
                  go_experimental_protein, catalytic_status)

  recomputed <- .pe_apply_term_sets(evidence, resolved_class)

  # Columns the recompute owns: replace these on the master, keep everything else (taxonomy,
  # identifiers, labels, membership flags) from the shipped table.
  recomputed_cols <- c(
    "acts_on_protein", "acts_on_nonprotein", "nonprotein_substrate_type",
    "dual_protein_nonprotein", "substrate_call", "substrate_evidence",
    "substrate_concordance", "substrate_decider", "ec_protein", "ec_nonprotein",
    "go_protein", "go_nonprotein", "n_evidence_dimensions", "evidence_tier",
    "curated_core", "is_catalytic_background", "is_protein_catalytic_background")

  master <- .pe_load(dataset)
  keep_attr <- attr(master, "term_set_md5")
  master_kept <- master[, setdiff(names(master), setdiff(recomputed_cols, "ensembl_gene_id")), drop = FALSE]
  recompute_join <- recomputed[, intersect(c("ensembl_gene_id", recomputed_cols), names(recomputed)), drop = FALSE]

  out <- dplyr::left_join(master_kept, recompute_join, by = dplyr::join_by(ensembl_gene_id))
  out <- out[, names(master), drop = FALSE]   # restore the shipped column order
  if (!is.null(keep_attr)) attr(out, "term_set_md5") <- keep_attr
  .pe_filter(out, mode, substrate)
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
