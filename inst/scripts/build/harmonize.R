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
