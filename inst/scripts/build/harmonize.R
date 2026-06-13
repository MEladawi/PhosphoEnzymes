# Map the engine's kinase table to the shipped `human_kinases` package schema.
# The engine keeps its rich, internally-named columns (used by the QC report); this is the
# presentation layer that renames to the cross-class vocabulary and selects the final, ordered
# column set.
#
# Substrate is carried as co-equal parallel columns -- acts_on_protein, acts_on_nonprotein, and
# the pipe-delimited nonprotein_substrate_type (empty = protein-only) -- so a downstream
# non-protein-pipeline filter never silently drops a dual enzyme (e.g. PIK3CA). The granular
# label stays in substrate_subtype. All of these are computed in the gate; this layer only
# renames and orders.

harmonize_kinases_to_package_schema <- function(kinases_table) {
  kinases_table |>
    transmute(
      ensembl_gene_id,
      symbol            = hgnc_symbol,
      acts_on_protein,
      acts_on_nonprotein,
      nonprotein_substrate_type,
      substrate_subtype = kinase_type,
      dual_protein_nonprotein,
      catalytic_status,
      n_evidence_dimensions,
      evidence_tier,
      curated_core,
      is_catalytic_background,
      in_structural_catalog,
      is_protein_kinase_ec,
      go_experimental,
      has_uniprot_kw    = is_uniprot_kw_kinase,
      supplementary_support,
      membership_basis,
      kinase_family,
      classification_reason,
      # taxonomy + breadth (rich master, not part of the thin summary)
      kinase_group, kinase_subfamily, derived_family, uniprot_protein_family,
      n_membership_sources, is_pseudogene,
      # identifiers + metadata
      hgnc_id, gene_name, entrez_id, uniprot_ids, prev_symbol, alias_symbol,
      enzyme_id_EC, ec_kinase_subclass,
      hgnc_gene_group, locus_type, chromosomal_location, mane_select_transcript, iuphar_id,
      hgnc_kinase_gene_group,
      # per-source provenance flags (Axis 1 rolls these up; kept for per-source filtering)
      is_pkinfam, is_manning, is_kinhub, is_go_kinase_activity, is_ec_kinase,
      is_idg_dark_kinase) |>
    tibble::as_tibble()
}

# Map the phosphatase engine table to the shipped `human_phosphatases` schema -- parallel to the
# kinase master (same substrate/evidence columns) but with phosphatase taxonomy (fold/family/
# subfamily) and per-source flags. The engine already emits the package columns; this only renames
# the symbol column and promotes it to second position.
harmonize_phosphatases_to_package_schema <- function(phosphatases_table) {
  phosphatases_table |>
    rename(symbol = hgnc_symbol) |>
    relocate(symbol, .after = ensembl_gene_id) |>
    tibble::as_tibble()
}

# Semicolon-joined names of the per-source flags that are TRUE for each row.
.evidence_sources_string <- function(df, flag_columns, labels) {
  pmap_chr(df[flag_columns], function(...) paste(labels[c(...)], collapse = ";"))
}

KINASE_SOURCE_FLAGS  <- c("is_pkinfam", "is_manning", "is_kinhub", "is_go_kinase_activity",
                          "is_ec_kinase", "has_uniprot_kw", "is_idg_dark_kinase")
KINASE_SOURCE_LABELS <- c("pkinfam", "Manning", "KinHub", "GO", "EC", "UniProtKW", "IDG")
PHOSPHATASE_SOURCE_FLAGS  <- c("is_chen", "is_hgnc_phosphatase_group", "is_go_phosphatase_activity",
                               "is_phosphatase_ec", "is_uniprot_kw_phosphatase")
PHOSPHATASE_SOURCE_LABELS <- c("Chen", "HGNC_groups", "GO", "EC", "UniProtKW")

# The thin, class-agnostic unified summary derived from the two masters (one source of truth per
# class). Shared columns only -- no family/group column (the taxonomy vocabularies are disjoint;
# join to a master to recover them). Ships all rows incl. Provisional, each carrying evidence_tier.
build_unified_summary <- function(kinases_pkg, phosphatases_pkg) {
  thin <- function(df, cls, flags, labels) {
    df |> transmute(
      ensembl_gene_id, symbol, regulator_class = cls,
      acts_on_protein, acts_on_nonprotein, nonprotein_substrate_type, dual_protein_nonprotein,
      catalytic_status, n_evidence_dimensions,
      evidence_sources = .evidence_sources_string(df, flags, labels),
      evidence_tier, curated_core, is_catalytic_background)
  }
  bind_rows(thin(kinases_pkg,      "kinase",      KINASE_SOURCE_FLAGS,      KINASE_SOURCE_LABELS),
            thin(phosphatases_pkg, "phosphatase", PHOSPHATASE_SOURCE_FLAGS, PHOSPHATASE_SOURCE_LABELS)) |>
    tibble::as_tibble()
}

# Per-ENSG Axis-1 deriving-source map, making the single-lineage reality machine-checkable.
build_membership_provenance <- function(kinases_pkg, phosphatases_pkg) {
  one <- function(df, cls) {
    df |> filter(!is.na(membership_basis)) |>
      transmute(ensembl_gene_id, regulator_class = cls,
                deriving_source = str_remove(membership_basis, "^[^:]+:"),
                membership_basis)
  }
  bind_rows(one(kinases_pkg, "kinase"), one(phosphatases_pkg, "phosphatase")) |>
    tibble::as_tibble()
}
