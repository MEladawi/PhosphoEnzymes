# Per-gene phosphatase classification (the protein-phosphatase gate) and assembly of the master
# table, mirroring classify_kinases so the two masters share a schema. The gate is function-first
# and symmetric with the kinase side: a non-protein substrate signal wins unless the gene carries
# protein-phosphatase activity GO (the substrate-accurate signal). Catalog membership (Chen / HGNC
# groups) and protein-EC are protein EVIDENCE but never substrate OVERRIDES -- phosphatase
# protein-EC (3.1.3.16/48) is promiscuously assigned (PTEN, MTM lipid phosphatases carry it), so it
# scores rigor (Axis 2) without authoring the substrate call. The override keys on IEA-inclusive
# protein GO, guarded by the PSPH reverse canary (symmetric with the kinase gate).

# Non-protein substrate classes by 4-digit EC (the codes that carry a clear substrate meaning).
PPASE_EC_LIPID      <- c("3.1.3.4", "3.1.3.36", "3.1.3.64", "3.1.3.66", "3.1.3.67",
                         "3.1.3.78", "3.1.3.86", "3.1.3.95")
PPASE_EC_NUCLEOTIDE <- c("3.1.3.5", "3.1.3.6", "3.1.3.7", "3.1.3.31", "3.1.3.91", "3.6.1.1")
PPASE_EC_CARB       <- c("3.1.3.9", "3.1.3.10", "3.1.3.11", "3.1.3.37", "3.1.3.46", "3.1.3.58")
# Chen families that are unambiguously non-protein, for genes lacking a specific non-protein EC.
CHEN_FAMILY_LIPID      <- c("Myotubularin", "Sac", "PTEN", "INPP4", "Lipin",
                            "Phosphatidic acid phosphatase", "IPP5")
# P5NPSP is a mixed HAD family (pyrimidine-5'-nucleotidase AND phosphoserine phosphatase, e.g.
# PSPH), so it is deliberately NOT mapped to nucleotide -- its members fall back to "other".
CHEN_FAMILY_NUCLEOTIDE <- c("5N", "NagD", "Deoxyribonucleotidases",
                            "cN-I nucleotidases", "cN-II nucleotidases")
CHEN_FAMILY_CARB       <- c("Glc-6-Pase")
# AP = the alkaline/acid phosphatase fold: broad-specificity small-molecule phosphatases (typed
# "other"). Ensures AP-fold genes lacking a specific EC (e.g. ALPG) still carry non-protein typing.
CHEN_FAMILY_OTHER      <- c("AP")

# First-matching non-protein substrate class for one gene: a specific non-protein EC, then a
# non-protein GO class, then a clearly-non-protein Chen family, then a generic "other" when some
# non-protein signal exists but no finer class maps (so the enum is never left blank when the gate
# tripped a non-protein signal -- the §4A consistency invariant). NA when there is no signal.
classify_nonprotein_phosphatase_class <- function(ec_codes, chen_family, chen_nonprotein,
                                                  go_lipid, go_nucleotide, go_carbohydrate, go_inositol) {
  if (any(ec_codes %in% PPASE_EC_LIPID))       return("lipid")
  if (any(ec_codes %in% PPASE_EC_NUCLEOTIDE))  return("nucleotide")
  if (any(ec_codes %in% PPASE_EC_CARB))        return("carbohydrate")
  if (isTRUE(go_lipid))        return("lipid")
  if (isTRUE(go_nucleotide))   return("nucleotide")
  if (isTRUE(go_carbohydrate)) return("carbohydrate")
  if (isTRUE(go_inositol))     return("other")
  if (!is.na(chen_family) && chen_family %in% CHEN_FAMILY_LIPID)      return("lipid")
  if (!is.na(chen_family) && chen_family %in% CHEN_FAMILY_NUCLEOTIDE) return("nucleotide")
  if (!is.na(chen_family) && chen_family %in% CHEN_FAMILY_CARB)       return("carbohydrate")
  if (!is.na(chen_family) && chen_family %in% CHEN_FAMILY_OTHER)      return("other")
  if (length(ec_codes) > 0 || isTRUE(chen_nonprotein)) return("other")
  NA_character_
}

classify_phosphatases <- function(universe_ensembl_ids, hgnc_bridge, go_phosphatase_sets, ec_phosphatase,
                                  membership, chen_facts, go_experimental_ids = character(0)) {
  add_set_flag <- function(data, ids, flag) {
    data |>
      left_join(tibble(ensembl_gene_id = unique(ids), present_in_set = TRUE),
                by = join_by(ensembl_gene_id)) |>
      mutate("{{ flag }}" := coalesce(present_in_set, FALSE)) |>
      select(-present_in_set)
  }

  hgnc_bridge$gene_metadata |>
    semi_join(tibble(ensembl_gene_id = universe_ensembl_ids), by = join_by(ensembl_gene_id)) |>
    # EC table (is_protein_phosphatase_ec + non-protein EC codes), joined in.
    left_join(ec_phosphatase$ec_table |>
                select(ensembl_gene_id, is_protein_phosphatase_ec, nonprotein_phosphatase_ec),
              by = join_by(ensembl_gene_id)) |>
    # Chen per-gene facts (taxonomy + substrate flags + catalytic status), joined in.
    left_join(chen_facts, by = join_by(ensembl_gene_id)) |>
    # one binary membership column per source
    add_set_flag(membership$chen,                 is_chen) |>
    add_set_flag(membership$hgnc_protein_group,    is_hgnc_phosphatase_group) |>
    add_set_flag(membership$go_umbrella,           is_go_phosphatase_activity) |>
    add_set_flag(membership$ec,                    is_phosphatase_ec) |>
    add_set_flag(membership$uniprot_keyword,       is_uniprot_kw_phosphatase) |>
    # GO sets used by the gate / typing
    add_set_flag(go_phosphatase_sets$protein_phosphatase_activity, in_protein_phosphatase_go) |>
    add_set_flag(go_phosphatase_sets$lipid_phosphatase,            in_go_lipid) |>
    add_set_flag(go_phosphatase_sets$nucleotide_phosphatase,       in_go_nucleotide) |>
    add_set_flag(go_phosphatase_sets$carbohydrate_phosphatase,     in_go_carbohydrate) |>
    add_set_flag(go_phosphatase_sets$inositol_phosphatase,         in_go_inositol) |>
    add_set_flag(go_experimental_ids,                              go_experimental) |>
    mutate(
      is_protein_phosphatase_ec = coalesce(is_protein_phosphatase_ec, FALSE),
      nonprotein_phosphatase_ec = map(nonprotein_phosphatase_ec, ~ if (is.null(.x)) character(0) else .x),
      chen_nonprotein_substrate = coalesce(chen_nonprotein_substrate, FALSE),
      chen_protein_substrate    = coalesce(chen_protein_substrate, FALSE),

      # Axis 1 (structural catalog): Chen phosphatome OR an HGNC protein-phosphatase gene group.
      in_structural_catalog = is_chen | is_hgnc_phosphatase_group,
      # Axis 2 (biochemical): a protein-specific 4-digit EC (3.1.3.16 / 3.1.3.48).
      n_evidence_dimensions = as.integer(in_structural_catalog) + as.integer(is_protein_phosphatase_ec),
      curated_core = n_evidence_dimensions >= 1L,
      has_uniprot_kw = is_uniprot_kw_phosphatase,
      supplementary_support = go_experimental | has_uniprot_kw,

      # Non-protein substrate class (co-equal signals: EC, GO, Chen). NA when no signal.
      nonprotein_class = pmap_chr(
        list(nonprotein_phosphatase_ec, chen_family, chen_nonprotein_substrate,
             in_go_lipid, in_go_nucleotide, in_go_carbohydrate, in_go_inositol),
        classify_nonprotein_phosphatase_class),
      nonprotein_substrate_type = coalesce(nonprotein_class, ""),
      acts_on_nonprotein = nzchar(nonprotein_substrate_type),

      # Protein-phosphatase evidence (assignment uses all functional evidence; scoring is separate).
      protein_phosphatase_evidence = in_protein_phosphatase_go | in_structural_catalog | is_protein_phosphatase_ec,
      # Lineage never overrides substrate: the override keys on protein-phosphatase activity GO only,
      # never on the promiscuous protein-EC or on catalog membership (the PTEN principle).
      nonprotein_wins = acts_on_nonprotein & !in_protein_phosphatase_go,
      protein_phosphatase = !nonprotein_wins & protein_phosphatase_evidence,
      acts_on_protein = protein_phosphatase,
      dual_protein_nonprotein = acts_on_protein & acts_on_nonprotein,

      # Granular substrate label (parallels the kinase substrate_subtype).
      phosphatase_type = case_when(
        nonprotein_wins & nonprotein_substrate_type == "lipid"        ~ "Lipid phosphatase",
        nonprotein_wins & nonprotein_substrate_type == "nucleotide"   ~ "Nucleotide phosphatase",
        nonprotein_wins & nonprotein_substrate_type == "carbohydrate" ~ "Carbohydrate/sugar phosphatase",
        nonprotein_wins                                               ~ "Other small-molecule phosphatase",
        protein_phosphatase                                          ~ "Protein phosphatase",
        .default                                                     = "Unclassified phosphatase"),

      catalytic_status = coalesce(catalytic_status, "active"),
      is_pseudophosphatase = coalesce(is_pseudophosphatase, FALSE),
      is_catalytic_background = catalytic_status == "active" & curated_core,

      membership_basis = case_when(
        is_chen                   ~ "reconstructed:Chen2017",
        is_hgnc_phosphatase_group ~ "reconstructed:HGNC_groups",
        .default                  = NA_character_),

      protein_evidence_label = case_when(
        in_structural_catalog & is_protein_phosphatase_ec & in_protein_phosphatase_go ~ "structural catalog + protein-EC + GO",
        in_structural_catalog & is_protein_phosphatase_ec                            ~ "structural catalog + protein-EC",
        in_structural_catalog & in_protein_phosphatase_go                            ~ "structural catalog + GO",
        is_protein_phosphatase_ec & in_protein_phosphatase_go                        ~ "protein-EC + GO",
        in_structural_catalog                                                        ~ "structural catalog",
        is_protein_phosphatase_ec                                                    ~ "protein-EC",
        in_protein_phosphatase_go                                                    ~ "GO protein-phosphatase activity",
        .default                                                                     = "lineage evidence"),
      classification_reason = case_when(
        nonprotein_wins     ~ str_c(phosphatase_type, ": non-protein substrate; no protein-phosphatase GO"),
        protein_phosphatase ~ str_c("protein phosphatase: ", protein_evidence_label),
        .default            = "non-catalytic / unclassified: in a phosphatase set but no protein or non-protein substrate evidence"),

      evidence_tier = case_when(
        n_evidence_dimensions == 2L                          ~ "Gold",
        n_evidence_dimensions == 1L & supplementary_support  ~ "Silver",
        n_evidence_dimensions == 1L                          ~ "Bronze",
        .default                                             = "Provisional"),

      n_membership_sources = is_chen + is_hgnc_phosphatase_group + is_go_phosphatase_activity +
                             is_phosphatase_ec + is_uniprot_kw_phosphatase,
      is_pseudogene = str_detect(coalesce(locus_type, ""), regex("pseudogene", ignore_case = TRUE))) |>
    transmute(
      ensembl_gene_id,
      hgnc_symbol = symbol, hgnc_id, gene_name = name,
      acts_on_protein, acts_on_nonprotein, nonprotein_substrate_type,
      substrate_subtype = phosphatase_type,
      dual_protein_nonprotein,
      catalytic_status, is_catalytic_background, is_pseudophosphatase,
      n_evidence_dimensions, evidence_tier, curated_core,
      in_structural_catalog, is_protein_phosphatase_ec,
      go_experimental, has_uniprot_kw, supplementary_support, membership_basis,
      classification_reason,
      phosphatase_fold = chen_fold, phosphatase_family = chen_family, phosphatase_subfamily = chen_subfamily,
      n_membership_sources, is_pseudogene,
      entrez_id, uniprot_ids, prev_symbol, alias_symbol,
      enzyme_id_EC = enzyme_id,
      hgnc_gene_group = gene_group, locus_type,
      chromosomal_location = location, mane_select_transcript = mane_select, iuphar_id = iuphar,
      is_chen, is_hgnc_phosphatase_group, is_go_phosphatase_activity, is_phosphatase_ec,
      is_uniprot_kw_phosphatase) |>
    arrange(desc(acts_on_protein), substrate_subtype, hgnc_symbol)
}
