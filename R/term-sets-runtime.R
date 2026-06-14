# Runtime twin of the build-engine term-set functions.
#
# The accessors expose a `term_sets=` override that re-types the shipped
# catalog under a user-supplied EC/GO term set without rebuilding. To do that
# off an installed package -- without sourcing inst/scripts/build -- the pure
# resolution and gate functions are ported here as internal (non-exported)
# copies. They are faithful to the build-engine versions in
# inst/scripts/build/term_sets.R and inst/scripts/build/utils.R: the build gate
# and this runtime path MUST stay in sync, so reclassifying with the default
# term set reproduces the shipped tables exactly. Change the gate in one place
# and you must change it in the other.
#
# The build resolves the GO axis into Ensembl-id sets via the pinned GMT. That
# GMT is not shipped in the installed package, so the runtime resolves the GO
# axis at the accession level instead (see .pe_resolve_go_accessions_class) and
# intersects it with each gene's pre-computed annotated accessions. The result
# is identical for any term set drawn from the build-time candidate accessions,
# which the round-trip identity test confirms on the default set.
#
# These functions use dplyr / purrr / stringr, which are Suggests (not
# Imports). They are reached ONLY through the `term_sets=` recompute path, which
# guards on those packages up front; the default (no `term_sets=`) accessors
# never touch this file.

# ---- ported from inst/scripts/build/utils.R --------------------------------

# Split a pipe-delimited field into a trimmed character vector (drops blanks).
# Literal double-quotes are stripped first so a quote-wrapped multi-value field
# does not leave a stray quote stuck to the first and last tokens.
.pe_split_pipe_delimited <- function(field_value) {
  if (is.na(field_value) || field_value == "") {
    character(0)
  } else {
    cleaned <- stringr::str_remove_all(field_value, '"')
    parts <- stringr::str_trim(
      stringr::str_split(cleaned, stringr::fixed("|"))[[1]])
    parts[parts != ""]
  }
}

# ---- ported from inst/scripts/build/term_sets.R ----------------------------

# A gene's EC code reduced to its 3-field subclass ("2.7.10.1" -> "2.7.10"); NA
# if shorter.
.pe_ec_subclass_of <- function(code) {
  fields <- stringr::str_split(code, stringr::fixed("."))[[1]]
  if (length(fields) >= 3) paste(fields[1:3], collapse = ".") else NA_character_
}

# TRUE if any of a gene's EC codes matches any rule. `rules` is a tibble with
# `code` (a term_id such as "2.7.10.-" or "3.1.3.16") and `scope`. subclass
# rules match by 3-field prefix; exact rules match the full 4-field code.
.pe_matches_ec_rules <- function(gene_codes, rules) {
  if (length(gene_codes) == 0 || nrow(rules) == 0) return(FALSE)
  subclass_targets <- stringr::str_remove(
    rules$code[rules$scope == "subclass"], "\\.-$")
  exact_targets <- rules$code[rules$scope == "exact"]
  gene_subclasses <- unique(stats::na.omit(
    purrr::map_chr(gene_codes, .pe_ec_subclass_of)))
  any(gene_subclasses %in% subclass_targets) ||
    any(gene_codes %in% exact_targets)
}

# GMT-free EC matcher resolver for one class. .pe_apply_term_sets reads only the
# `ec_*` matchers from a resolved class (EC typing needs no gene-set file), so
# the recompute path resolves them straight from the term-set EC table. Mirrors
# the `ec_*` slots of the build resolver.
.pe_resolve_ec_class <- function(ec_tbl) {
  ec_rigor <- dplyr::filter(ec_tbl, role == "rigor+substrate")
  list(
    ec_rigor      = dplyr::select(ec_rigor, code = term_id, scope),
    ec_protein    = ec_rigor |>
      dplyr::filter(substrate == "protein") |>
      dplyr::select(code = term_id, scope),
    ec_nonprotein = ec_rigor |>
      dplyr::filter(substrate == "nonprotein") |>
      dplyr::select(code = term_id, scope, subtype = substrate_subtype))
}

# GMT-free GO resolver for one class. The build builds GO id-sets from the GMT;
# here we keep the protein / non-protein / per-subtype GO *accession* lists read
# straight off the term-set table, and the recompute path intersects them with
# each gene's annotated accessions. Mirrors, in result, the build's GO id-sets.
.pe_resolve_go_accessions_class <- function(go_tbl) {
  go_rigor <- dplyr::filter(go_tbl, role == "rigor+substrate")
  go_np    <- dplyr::filter(go_rigor, substrate == "nonprotein")
  subtypes <- unique(go_np$substrate_subtype[
    !is.na(go_np$substrate_subtype) & go_np$substrate_subtype != ""])
  list(
    go_protein_accessions =
      dplyr::filter(go_rigor, substrate == "protein")$term_id,
    go_nonprotein_accessions = go_np$term_id,
    go_nonprotein_subtype_accessions = purrr::set_names(
      purrr::map(subtypes,
        \(st) go_np$term_id[go_np$substrate_subtype == st]),
      subtypes))
}

# Per-gene EC axis flags against a resolved class. Returns the rigor flag (E),
# the protein and nonprotein flags, and the firing nonprotein substrate
# subtype(s).
.pe_ec_axis_flags <- function(gene_codes, resolved_class) {
  np_rows <- resolved_class$ec_nonprotein
  fired_np <- if (length(gene_codes) == 0 || nrow(np_rows) == 0) {
    character(0)
  } else {
    np_rows$subtype[purrr::map_lgl(seq_len(nrow(np_rows)),
      \(i) .pe_matches_ec_rules(gene_codes, np_rows[i, c("code", "scope")]))]
  }
  list(
    ec_rigor      = .pe_matches_ec_rules(gene_codes, resolved_class$ec_rigor),
    ec_protein    = .pe_matches_ec_rules(gene_codes, resolved_class$ec_protein),
    ec_nonprotein = .pe_matches_ec_rules(
      gene_codes, resolved_class$ec_nonprotein),
    nonprotein_subtypes = unique(fired_np[!is.na(fired_np) & fired_np != ""]))
}

# Compute every term-set-dependent column from per-gene evidence and a resolved
# class. This is the single gate shared (in twinned form) by the build and the
# accessor term_sets= override.
.pe_apply_term_sets <- function(evidence, resolved_class) {
  ec_flags <- purrr::map(
    evidence$all_ec_codes, \(codes) .pe_ec_axis_flags(codes, resolved_class))
  evidence |> dplyr::mutate(
    ec_rigor      = purrr::map_lgl(ec_flags, "ec_rigor"),
    ec_protein    = purrr::map_lgl(ec_flags, "ec_protein"),
    ec_nonprotein = purrr::map_lgl(ec_flags, "ec_nonprotein"),
    ec_nonprotein_subtype = purrr::map_chr(
      ec_flags, \(f) paste(f$nonprotein_subtypes, collapse = "|")),

    # RIGOR (substrate-blind): L + E, where E = any class EC (protein or
    # non-protein).
    n_evidence_dimensions =
      as.integer(in_structural_catalog) + as.integer(ec_rigor),
    curated_core = n_evidence_dimensions >= 1L,
    evidence_tier = dplyr::case_when(
      n_evidence_dimensions == 2L                          ~ "Gold",
      n_evidence_dimensions == 1L & supplementary_support  ~ "Silver",
      n_evidence_dimensions == 1L                          ~ "Bronze",
      .default                                             = "Provisional"),

    # SUBSTRATE: co-equal flags; no lineage default.
    acts_on_protein    = go_protein | ec_protein,
    acts_on_nonprotein = go_nonprotein | ec_nonprotein | chen_nonprotein,
    nonprotein_substrate_type = purrr::pmap_chr(
      list(go_nonprotein_subtype, ec_nonprotein_subtype, acts_on_nonprotein),
      \(go_st, ec_st, any_np) {
        subs <- unique(c(.pe_split_pipe_delimited(go_st),
                         .pe_split_pipe_delimited(ec_st)))
        subs <- subs[nzchar(subs)]
        if (!any_np) {
          ""
        } else if (length(subs)) {
          paste(subs, collapse = "|")
        } else {
          "other"
        }
      }),
    dual_protein_nonprotein = acts_on_protein & acts_on_nonprotein,
    substrate_call = dplyr::case_when(
      acts_on_protein & acts_on_nonprotein ~ "dual",
      acts_on_protein                      ~ "protein",
      acts_on_nonprotein                   ~ "nonprotein",
      .default                             = "untyped"),

    # PROVENANCE.
    substrate_evidence = purrr::pmap_chr(
      list(go_protein | go_nonprotein, chen_nonprotein,
           ec_protein | ec_nonprotein),
      \(go, chen, ec) paste(c("GO", "Chen", "EC")[c(go, chen, ec)],
                            collapse = "+")),
    substrate_decider = dplyr::case_when(
      (go_protein | go_nonprotein) & go_experimental_protein ~ "GO-exp",
      (go_protein | go_nonprotein) ~ "GO-elec",
      chen_nonprotein ~ "Chen-flag",
      ec_protein | ec_nonprotein ~ "EC",
      .default = "precedence-default"),
    substrate_concordance = purrr::pmap_chr(
      list(go_protein, go_nonprotein, ec_protein, ec_nonprotein,
           chen_nonprotein, substrate_call),
      \(gp, gn, ep, en, cn, call) {
        if (call == "untyped") return("untyped")
        n_sources <- sum(gp | gn, ep | en, cn)
        if (n_sources <= 1) "single" else "concordant"
      }),

    # ENRICHMENT BACKGROUNDS.
    is_catalytic_background = catalytic_status == "active" & curated_core,
    is_protein_catalytic_background =
      catalytic_status == "active" & curated_core & acts_on_protein)
}
