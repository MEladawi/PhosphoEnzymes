# helper_code/classify.R
# Per-gene classification by enzymatic function (the protein-kinase gate) and assembly
# of the master table with one binary membership column per source.
#
#   universe_ensembl_ids : sorted union of all membership legs
#   hgnc_bridge          : provides $gene_metadata
#   go_sets              : load_go_functional_sets() output (protein gate + non-protein classes)
#   ec                   : load_ec_kinome() output (provides $ec_table)
#   membership           : named list of per-source Ensembl vectors
#   taxonomy             : build_kinase_taxonomy() output (group/family/subfamily maps)

classify_kinases <- function(universe_ensembl_ids, hgnc_bridge, go_sets, ec, membership, taxonomy) {
  gene_metadata <- hgnc_bridge$gene_metadata
  ec_table      <- ec$ec_table
  metadata_row_of_ensembl <- setNames(seq_len(nrow(gene_metadata)), gene_metadata$ensembl_gene_id)
  ec_row_of_ensembl       <- setNames(seq_len(nrow(ec_table)),      ec_table$ensembl_gene_id)

  # Non-catalytic exclusions for the HGNC kinase-gene-group annotation flag.
  noncatalytic_pattern <- regex(
    "anchoring|phosphatase|activator|non-catalytic|guanylate kinases|subunits|MOB|binding RTK",
    ignore_case = TRUE)

  classify_single_gene <- function(ensembl_gene_id) {
    gene_row <- gene_metadata[metadata_row_of_ensembl[[ensembl_gene_id]], ]
    ec_row   <- ec_table[ec_row_of_ensembl[[ensembl_gene_id]], ]
    ec_subclasses_present <- ec_row$ec_subclasses[[1]]

    in_pkinfam              <- ensembl_gene_id %in% membership$pkinfam
    in_manning              <- ensembl_gene_id %in% membership$manning
    in_kinhub               <- ensembl_gene_id %in% membership$kinhub
    in_go_kinase_activity   <- ensembl_gene_id %in% membership$go_umbrella
    in_ec_kinase            <- ensembl_gene_id %in% membership$ec
    in_uniprot_keyword      <- ensembl_gene_id %in% membership$uniprot_keyword
    in_idg_dark             <- ensembl_gene_id %in% membership$idg_dark
    in_protein_kinase_go    <- ensembl_gene_id %in% go_sets$protein_kinase_activity
    has_protein_kinase_ec   <- isTRUE(ec_row$is_protein_kinase_ec)

    # First matching non-protein functional class (priority order), else NA.
    nonprotein_class <- dplyr::case_when(
      ensembl_gene_id %in% go_sets$lipid_kinase              ~ "Lipid kinase",
      ensembl_gene_id %in% go_sets$inositol_phosphate_kinase ~ "Inositol-phosphate kinase",
      ensembl_gene_id %in% go_sets$carbohydrate_kinase       ~ "Carbohydrate/sugar kinase",
      ensembl_gene_id %in% go_sets$nucleotide_kinase         ~ "Nucleotide/nucleoside kinase",
      ensembl_gene_id %in% go_sets$creatine_kinase           ~ "Creatine kinase",
      TRUE                                                   ~ NA_character_)

    # Protein-kinase evidence comes only from protein-specific resources.
    has_protein_kinase_evidence <- in_protein_kinase_go || in_pkinfam || in_manning ||
                                   in_kinhub || has_protein_kinase_ec

    # Gate: a non-protein class wins only if the gene is NOT flagged a protein kinase by a
    # reliable signal -- GO protein-kinase activity OR a protein-kinase EC subclass (2.7.10-14,
    # which are never lipid). pkinfam/Manning/KinHub are deliberately NOT trusted here because
    # they lump the PI3/PI4 (lipid+protein) family together.
    if (!is.na(nonprotein_class) && !in_protein_kinase_go && !has_protein_kinase_ec) {
      kinase_type <- nonprotein_class; is_protein_kinase <- FALSE
    } else if (has_protein_kinase_evidence) {
      kinase_type <- "Protein kinase"; is_protein_kinase <- TRUE
    } else {                                                    # non-protein fallback typed from EC subclass
      kinase_type <- dplyr::case_when(
        "2.7.3" %in% ec_subclasses_present ~ "Creatine/phosphagen kinase",
        "2.7.4" %in% ec_subclasses_present ~ "Nucleotide kinase",
        "2.7.6" %in% ec_subclasses_present ~ "Diphosphokinase",
        "2.7.2" %in% ec_subclasses_present ~ "Carboxyl-group kinase",
        "2.7.1" %in% ec_subclasses_present ~ "Small-molecule kinase (EC 2.7.1)",
        "2.7.9" %in% ec_subclasses_present ~ "Dikinase (EC 2.7.9)",
        TRUE                               ~ "Other/unclassified kinase")
      is_protein_kinase <- FALSE
    }

    membership_source_count <- sum(in_pkinfam, in_manning, in_kinhub, in_go_kinase_activity,
                                   in_ec_kinase, in_uniprot_keyword, in_idg_dark)
    gene_group_terms <- split_pipe_delimited(gene_row$gene_group)

    tibble(
      ensembl_gene_id = ensembl_gene_id,
      hgnc_symbol = gene_row$symbol, hgnc_id = gene_row$hgnc_id, gene_name = gene_row$name,
      kinase_type = kinase_type, protein_kinase = is_protein_kinase,
      kinase_group     = lookup_in_named_vector(taxonomy$group, ensembl_gene_id),
      kinase_family    = lookup_in_named_vector(taxonomy$family, ensembl_gene_id),
      kinase_subfamily = lookup_in_named_vector(taxonomy$subfamily, ensembl_gene_id),
      uniprot_protein_family = lookup_in_named_vector(taxonomy$uniprot_family_raw, ensembl_gene_id),
      dual_protein_and_nonprotein = is_protein_kinase && !is.na(nonprotein_class),
      confidence = if (membership_source_count >= 2 || in_ec_kinase) "high" else "low (single-source)",
      n_membership_sources = membership_source_count,
      curated_core = (in_pkinfam || in_manning || in_kinhub || in_ec_kinase ||
                      in_uniprot_keyword || in_idg_dark),
      is_pseudogene = str_detect(gene_row$locus_type %||% "", regex("pseudogene", ignore_case = TRUE)),
      entrez_id = gene_row$entrez_id, uniprot_ids = gene_row$uniprot_ids,
      prev_symbol = gene_row$prev_symbol, alias_symbol = gene_row$alias_symbol,
      enzyme_id_EC = gene_row$enzyme_id,
      ec_kinase_subclass = paste(ec_row$matched_kinase_subclasses[[1]], collapse = ", "),
      hgnc_gene_group = gene_row$gene_group, locus_type = gene_row$locus_type,
      chromosomal_location = gene_row$location, mane_select_transcript = gene_row$mane_select,
      iuphar_id = gene_row$iuphar,
      hgnc_kinase_gene_group = any(str_detect(gene_group_terms, regex("kinase", ignore_case = TRUE)) &
                                   !str_detect(gene_group_terms, noncatalytic_pattern)),
      # ---- one binary column per source ----
      is_pkinfam = in_pkinfam, is_manning = in_manning, is_kinhub = in_kinhub,
      is_go_kinase_activity = in_go_kinase_activity, is_ec_kinase = in_ec_kinase,
      is_uniprot_kw_kinase = in_uniprot_keyword, is_idg_dark_kinase = in_idg_dark)
  }

  map_dfr(universe_ensembl_ids, classify_single_gene) %>%
    arrange(desc(protein_kinase), kinase_type, hgnc_symbol)
}
