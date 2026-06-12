# helper_code/classify.R
# Per-gene classification by enzymatic function (the protein-kinase gate) and assembly of
# the master table with one binary membership column per source. A single dplyr pipeline
# over the universe: source/GO membership is added by tidy joins (not %in%), then the gate
# and derived columns are computed column-wise.
#
#   universe_ensembl_ids : sorted union of all membership legs
#   hgnc_bridge          : provides $gene_metadata (metadata also rides along in ec$ec_table)
#   go_sets              : load_go_functional_sets() output (protein gate + non-protein classes)
#   ec                   : load_ec_kinome() output; $ec_table already carries all gene metadata
#   membership           : named list of per-source Ensembl vectors
#   taxonomy             : build_kinase_taxonomy() output (group/family/subfamily maps)

classify_kinases <- function(universe_ensembl_ids, hgnc_bridge, go_sets, ec, membership, taxonomy,
                             go_experimental_ids = character(0)) {
  noncatalytic_pattern <- regex(
    "anchoring|phosphatase|activator|non-catalytic|guanylate kinases|subunits|MOB|binding RTK",
    ignore_case = TRUE)

  # Add a TRUE/FALSE column named {{ flag }} marking rows whose ensembl_gene_id is in `ids`
  # (a tidy left join instead of a bare `%in%` membership test). `flag` is a bare column name.
  add_set_flag <- function(data, ids, flag) {
    data |>
      left_join(tibble(ensembl_gene_id = unique(ids), present_in_set = TRUE),
                by = join_by(ensembl_gene_id)) |>
      mutate("{{ flag }}" := coalesce(present_in_set, FALSE)) |>
      select(-present_in_set)
  }

  ec$ec_table |>
    semi_join(tibble(ensembl_gene_id = universe_ensembl_ids), by = join_by(ensembl_gene_id)) |>
    # one binary membership column per source
    add_set_flag(membership$pkinfam,         is_pkinfam) |>
    add_set_flag(membership$manning,         is_manning) |>
    add_set_flag(membership$kinhub,          is_kinhub) |>
    add_set_flag(membership$go_umbrella,     is_go_kinase_activity) |>
    # is_ec_kinase already exists as a column in ec$ec_table (== membership$ec), so reuse it.
    add_set_flag(membership$uniprot_keyword, is_uniprot_kw_kinase) |>
    add_set_flag(membership$idg_dark,        is_idg_dark_kinase) |>
    # GO functional sets used by the gate
    add_set_flag(go_sets$protein_kinase_activity,   in_protein_kinase_go) |>
    add_set_flag(go_sets$lipid_kinase,              in_lipid) |>
    add_set_flag(go_sets$inositol_phosphate_kinase, in_inositol) |>
    add_set_flag(go_sets$carbohydrate_kinase,       in_carbohydrate) |>
    add_set_flag(go_sets$nucleotide_kinase,         in_nucleotide) |>
    add_set_flag(go_sets$creatine_kinase,           in_creatine) |>
    # Provenance proxy: non-electronic (experimental/curated) GO kinase-activity support.
    add_set_flag(go_experimental_ids,               go_experimental) |>
    mutate(
      n_membership_sources = is_pkinfam + is_manning + is_kinhub + is_go_kinase_activity +
                             is_ec_kinase + is_uniprot_kw_kinase + is_idg_dark_kinase,

      # Independent evidence axes: the rigor metric counts independent KINDS of confirmation,
      # not independent databases. Exactly two evidence TYPES answer genuinely different questions:
      #   Axis 1 -- structural/evolutionary catalog: does the gene carry the kinase sequence
      #             family? pkinfam / Manning / KinHub all answer this from overlapping
      #             scholarship, so they roll up into ONE axis -- counting them separately would
      #             inflate correlated agreement (the database-fame artifact that also demotes
      #             GO / UniProt).
      #   Axis 2 -- biochemical: a protein-specific EC (2.7.10-14) proving protein-directed
      #             catalysis. Non-protein EC still types substrate below but never scores here.
      # So axes == 2 means structure AND biochemistry independently agree (what "Gold" certifies).
      in_structural_catalog       = is_pkinfam | is_manning | is_kinhub,
      n_independent_evidence_axes = as.integer(in_structural_catalog) +
                                    as.integer(is_protein_kinase_ec),
      # curated_core = standing in at least one independent axis (a structural catalog or a
      # protein-EC number). The comprehensive-only remainder (GO/UniProt/IDG-only) is
      # curated_core FALSE and tiers as Provisional; this is the strict-mode population.
      curated_core = n_independent_evidence_axes >= 1L,
      # Supplementary support: experimental GO support OR the reviewed UniProt kinase keyword.
      # Neither counts toward the axes; together they split Silver from Bronze among single-axis
      # genes. It cannot manufacture standing from zero axes -- a zero-axis gene stays Provisional.
      supplementary_support = go_experimental | is_uniprot_kw_kinase,

      # First matching non-protein functional class (priority order), else NA.
      nonprotein_class = case_when(
        in_lipid        ~ "Lipid kinase",
        in_inositol     ~ "Inositol-phosphate kinase",
        in_carbohydrate ~ "Carbohydrate/sugar kinase",
        in_nucleotide   ~ "Nucleotide/nucleoside kinase",
        in_creatine     ~ "Creatine kinase",
        .default        = NA_character_),
      protein_kinase_evidence = in_protein_kinase_go | is_pkinfam | is_manning | is_kinhub |
                                is_protein_kinase_ec,
      # A non-protein class wins only if no reliable protein signal (GO protein activity or a
      # protein-kinase EC subclass, never lipid); pkinfam/Manning/KinHub lump PI3/PI4 so are
      # not trusted to override here.
      nonprotein_wins = !is.na(nonprotein_class) & !in_protein_kinase_go & !is_protein_kinase_ec,
      protein_kinase  = !nonprotein_wins & protein_kinase_evidence,
      # EC-subclass fallback type (checking a code within each gene's EC subclass list).
      ec_fallback_type = case_when(
        map_lgl(ec_subclasses, \(x) "2.7.3" %in% x) ~ "Creatine/phosphagen kinase",
        map_lgl(ec_subclasses, \(x) "2.7.4" %in% x) ~ "Nucleotide kinase",
        map_lgl(ec_subclasses, \(x) "2.7.6" %in% x) ~ "Diphosphokinase",
        map_lgl(ec_subclasses, \(x) "2.7.2" %in% x) ~ "Carboxyl-group kinase",
        map_lgl(ec_subclasses, \(x) "2.7.1" %in% x) ~ "Small-molecule kinase (EC 2.7.1)",
        map_lgl(ec_subclasses, \(x) "2.7.9" %in% x) ~ "Dikinase (EC 2.7.9)",
        .default                                    = "Other/unclassified kinase"),
      kinase_type = case_when(
        nonprotein_wins ~ nonprotein_class,
        protein_kinase  ~ "Protein kinase",
        .default        = ec_fallback_type),

      dual_protein_and_nonprotein = protein_kinase & !is.na(nonprotein_class),
      # evidence_tier: a documented PRIORITIZATION HEURISTIC over the two axes plus
      # supplementary support -- NOT a probability, evidence count, or confidence score. Gold
      # requires BOTH axes (structure and biochemistry agree); Silver/Bronze split the one-axis
      # genes by supplementary_support (experimental GO or reviewed UniProt keyword); Provisional
      # is the comprehensive-only remainder. GO / keyword never reach Gold by design -- they share
      # literature provenance with the catalogs, so admitting them would reintroduce the coupling
      # the axis count exists to exclude.
      evidence_tier = case_when(
        n_independent_evidence_axes == 2L                         ~ "Gold",
        n_independent_evidence_axes == 1L & supplementary_support ~ "Silver",
        n_independent_evidence_axes == 1L                         ~ "Bronze",
        .default                                                  = "Provisional"),
      is_pseudogene = str_detect(coalesce(locus_type, ""), regex("pseudogene", ignore_case = TRUE)),
      ec_kinase_subclass = map_chr(matched_kinase_subclasses, \(x) paste(x, collapse = ", ")),
      # Manning taxonomy (named-vector maps keyed by Ensembl ID); NA where absent.
      kinase_group           = unname(taxonomy$group[ensembl_gene_id]),
      kinase_family          = unname(taxonomy$family[ensembl_gene_id]),
      kinase_subfamily       = unname(taxonomy$subfamily[ensembl_gene_id]),
      uniprot_protein_family = unname(taxonomy$uniprot_family_raw[ensembl_gene_id]),
      # Fallback family descriptor for genes with no Manning kinase_family (typically
      # UniProt/GO-only atypical kinases). NON-Manning: the UniProt parsed family tier
      # (built once by parse_uniprot_protein_family via build_kinase_taxonomy), else the GO
      # functional class. Blank where a Manning kinase_family exists.
      derived_family = if_else(
        !is.na(kinase_family) & kinase_family != "",
        NA_character_,
        coalesce(unname(taxonomy$uniprot_family_tier[ensembl_gene_id]), nonprotein_class)),
      hgnc_kinase_gene_group = map_lgl(gene_group, \(group_field) {
        terms <- split_pipe_delimited(group_field)
        any(str_detect(terms, regex("kinase", ignore_case = TRUE)) & !str_detect(terms, noncatalytic_pattern))
      })) |>
    # Final column set and order (drops all intermediate flag/gate columns).
    transmute(
      ensembl_gene_id,
      hgnc_symbol = symbol, hgnc_id, gene_name = name,
      kinase_type, protein_kinase,
      kinase_group, kinase_family, derived_family, kinase_subfamily, uniprot_protein_family,
      dual_protein_and_nonprotein, evidence_tier, n_independent_evidence_axes,
      go_experimental, supplementary_support,
      n_membership_sources, curated_core,
      is_pseudogene,
      entrez_id, uniprot_ids, prev_symbol, alias_symbol,
      enzyme_id_EC = enzyme_id, ec_kinase_subclass,
      hgnc_gene_group = gene_group, locus_type,
      chromosomal_location = location, mane_select_transcript = mane_select, iuphar_id = iuphar,
      hgnc_kinase_gene_group,
      is_pkinfam, is_manning, is_kinhub, is_go_kinase_activity, is_ec_kinase,
      is_uniprot_kw_kinase, is_idg_dark_kinase) |>
    arrange(desc(protein_kinase), kinase_type, hgnc_symbol)
}
