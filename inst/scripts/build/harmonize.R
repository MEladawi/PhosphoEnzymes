# Map the engine's kinase table to the shipped `human_kinases` package schema.
# The engine keeps its rich, internally-named columns (used by the QC report); this is
# the presentation layer that renames to the cross-class vocabulary, flattens the granular
# kinase_type into the controlled substrate_type, and selects the final, ordered column set.
#
# substrate_type is the 5-value controlled vocabulary {protein, lipid, nucleotide,
# carbohydrate, other}. acts_on_protein is TRUE exactly when substrate_type == "protein",
# because the gate types a gene "Protein kinase" iff it passes the protein gate.

harmonize_kinases_to_package_schema <- function(kinases_table) {
  kinases_table |>
    mutate(
      substrate_type = case_when(
        kinase_type == "Protein kinase"                                  ~ "protein",
        kinase_type == "Lipid kinase"                                    ~ "lipid",
        kinase_type %in% c("Nucleotide/nucleoside kinase",
                           "Nucleotide kinase")                          ~ "nucleotide",
        kinase_type == "Carbohydrate/sugar kinase"                       ~ "carbohydrate",
        .default                                                         = "other")) |>
    transmute(
      ensembl_gene_id,
      symbol            = hgnc_symbol,
      acts_on_protein   = protein_kinase,
      substrate_type,
      substrate_subtype = kinase_type,
      n_evidence_dimensions,
      evidence_tier,
      curated_core,
      in_structural_catalog,
      is_protein_kinase_ec,
      go_experimental,
      has_uniprot_kw    = is_uniprot_kw_kinase,
      supplementary_support,
      kinase_family,
      classification_reason,
      # taxonomy + bifunctional flag + breadth (rich master, not part of the thin summary)
      kinase_group, kinase_subfamily, derived_family, uniprot_protein_family,
      dual_protein_nonprotein, n_membership_sources, is_pseudogene,
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
