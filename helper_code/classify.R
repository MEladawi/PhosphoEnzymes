# helper_code/classify.R
# Per-gene classification by enzymatic function (the protein-kinase gate) and assembly of
# the master table with one binary membership column per source. Fully vectorized: a single
# dplyr pipeline over the universe (no row-by-row loop).
#
#   universe_ensembl_ids : sorted union of all membership legs
#   hgnc_bridge          : provides $gene_metadata (unused here; metadata comes via ec$ec_table)
#   go_sets              : load_go_functional_sets() output (protein gate + non-protein classes)
#   ec                   : load_ec_kinome() output; $ec_table already carries all gene metadata
#   membership           : named list of per-source Ensembl vectors
#   taxonomy             : build_kinase_taxonomy() output (group/family/subfamily maps)

classify_kinases <- function(universe_ensembl_ids, hgnc_bridge, go_sets, ec, membership, taxonomy) {
  # Non-catalytic exclusions for the HGNC kinase-gene-group annotation flag.
  noncatalytic_pattern <- regex(
    "anchoring|phosphatase|activator|non-catalytic|guanylate kinases|subunits|MOB|binding RTK",
    ignore_case = TRUE)

  ec$ec_table %>%
    filter(ensembl_gene_id %in% universe_ensembl_ids) %>%
    mutate(
      # --- one binary membership column per source ---
      is_pkinfam            = ensembl_gene_id %in% membership$pkinfam,
      is_manning            = ensembl_gene_id %in% membership$manning,
      is_kinhub             = ensembl_gene_id %in% membership$kinhub,
      is_go_kinase_activity = ensembl_gene_id %in% membership$go_umbrella,
      is_ec_kinase          = ensembl_gene_id %in% membership$ec,
      is_uniprot_kw_kinase  = ensembl_gene_id %in% membership$uniprot_keyword,
      is_idg_dark_kinase    = ensembl_gene_id %in% membership$idg_dark,
      n_membership_sources  = is_pkinfam + is_manning + is_kinhub + is_go_kinase_activity +
                              is_ec_kinase + is_uniprot_kw_kinase + is_idg_dark_kinase,
      curated_core          = is_pkinfam | is_manning | is_kinhub | is_ec_kinase |
                              is_uniprot_kw_kinase | is_idg_dark_kinase,

      # --- functional gate ---
      in_protein_kinase_go = ensembl_gene_id %in% go_sets$protein_kinase_activity,
      # First matching non-protein functional class (priority order), else NA.
      nonprotein_class = case_when(
        ensembl_gene_id %in% go_sets$lipid_kinase              ~ "Lipid kinase",
        ensembl_gene_id %in% go_sets$inositol_phosphate_kinase ~ "Inositol-phosphate kinase",
        ensembl_gene_id %in% go_sets$carbohydrate_kinase       ~ "Carbohydrate/sugar kinase",
        ensembl_gene_id %in% go_sets$nucleotide_kinase         ~ "Nucleotide/nucleoside kinase",
        ensembl_gene_id %in% go_sets$creatine_kinase           ~ "Creatine kinase",
        TRUE                                                   ~ NA_character_),
      protein_kinase_evidence = in_protein_kinase_go | is_pkinfam | is_manning | is_kinhub |
                                is_protein_kinase_ec,
      # A non-protein class wins only if no reliable protein signal (GO protein activity or a
      # protein-kinase EC subclass, never lipid); pkinfam/Manning/KinHub lump PI3/PI4 so are
      # not trusted to override here.
      nonprotein_wins = !is.na(nonprotein_class) & !in_protein_kinase_go & !is_protein_kinase_ec,
      protein_kinase  = !nonprotein_wins & protein_kinase_evidence,
      # EC-subclass fallback type for genes with neither a non-protein class nor protein evidence.
      ec_fallback_type = case_when(
        map_lgl(ec_subclasses, ~ "2.7.3" %in% .x) ~ "Creatine/phosphagen kinase",
        map_lgl(ec_subclasses, ~ "2.7.4" %in% .x) ~ "Nucleotide kinase",
        map_lgl(ec_subclasses, ~ "2.7.6" %in% .x) ~ "Diphosphokinase",
        map_lgl(ec_subclasses, ~ "2.7.2" %in% .x) ~ "Carboxyl-group kinase",
        map_lgl(ec_subclasses, ~ "2.7.1" %in% .x) ~ "Small-molecule kinase (EC 2.7.1)",
        map_lgl(ec_subclasses, ~ "2.7.9" %in% .x) ~ "Dikinase (EC 2.7.9)",
        TRUE                                      ~ "Other/unclassified kinase"),
      kinase_type = case_when(
        nonprotein_wins ~ nonprotein_class,
        protein_kinase  ~ "Protein kinase",
        TRUE            ~ ec_fallback_type),

      # --- derived annotation columns ---
      dual_protein_and_nonprotein = protein_kinase & !is.na(nonprotein_class),
      confidence = if_else(n_membership_sources >= 2 | is_ec_kinase, "high", "low (single-source)"),
      is_pseudogene = str_detect(coalesce(locus_type, ""), regex("pseudogene", ignore_case = TRUE)),
      ec_kinase_subclass = map_chr(matched_kinase_subclasses, ~ paste(.x, collapse = ", ")),
      # Manning taxonomy (named-vector maps keyed by Ensembl ID); NA where absent.
      kinase_group           = unname(taxonomy$group[ensembl_gene_id]),
      kinase_family          = unname(taxonomy$family[ensembl_gene_id]),
      kinase_subfamily       = unname(taxonomy$subfamily[ensembl_gene_id]),
      uniprot_protein_family = unname(taxonomy$uniprot_family_raw[ensembl_gene_id]),
      hgnc_kinase_gene_group = map_lgl(gene_group, function(group_field) {
        terms <- split_pipe_delimited(group_field)
        any(str_detect(terms, regex("kinase", ignore_case = TRUE)) & !str_detect(terms, noncatalytic_pattern))
      })) %>%
    # Final column set and order (drops all intermediate columns).
    transmute(
      ensembl_gene_id,
      hgnc_symbol = symbol, hgnc_id, gene_name = name,
      kinase_type, protein_kinase,
      kinase_group, kinase_family, kinase_subfamily, uniprot_protein_family,
      dual_protein_and_nonprotein, confidence, n_membership_sources, curated_core,
      is_pseudogene,
      entrez_id, uniprot_ids, prev_symbol, alias_symbol,
      enzyme_id_EC = enzyme_id, ec_kinase_subclass,
      hgnc_gene_group = gene_group, locus_type,
      chromosomal_location = location, mane_select_transcript = mane_select, iuphar_id = iuphar,
      hgnc_kinase_gene_group,
      is_pkinfam, is_manning, is_kinhub, is_go_kinase_activity, is_ec_kinase,
      is_uniprot_kw_kinase, is_idg_dark_kinase) %>%
    arrange(desc(protein_kinase), kinase_type, hgnc_symbol)
}
