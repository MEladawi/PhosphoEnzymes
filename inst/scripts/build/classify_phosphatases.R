# Per-gene phosphatase classification (the phosphatase gate) and assembly of the master table,
# mirroring classify_kinases so the two masters share a schema. RIGOR is scored substrate-blind --
# a gene stands in Axis 2 when it carries ANY class EC (protein or non-protein), the phosphatase EC
# allow-list / non-protein 3.1.3.x codes -- and SUBSTRATE is typed by the four co-equal flags (GO
# protein / GO non-protein / EC protein / EC non-protein, plus the Chen curated non-protein flag)
# with NO lineage default: a structural catalog tells you the gene is a phosphatase, never which
# substrate it dephosphorylates. The substrate flags + tiering both come from the single shared
# apply_term_sets() path, so the build gate and the accessor term_sets= override agree.
#
#   universe_ensembl_ids : sorted union of all membership legs
#   hgnc_bridge          : provides $gene_metadata (identifiers + HGNC fields)
#   go_phosphatase_sets  : load_go_phosphatase_sets() output -- phosphatase_activity_umbrella,
#                          protein_phosphatase_activity, nonprotein_all, nonprotein_by_subtype
#   ec_phosphatase       : load_ec_phosphatome() output; $ec_table carries all_ec_codes
#   membership           : named list of per-source Ensembl vectors
#   chen_facts           : load_chen_phosphatome()$facts_table -- taxonomy + substrate flags + status
#   resolved_phosphatase : resolve_term_sets()$phosphatase -- the resolved EC/GO matchers

# Curated regulatory roles for the catalytically inactive myotubularins. These pseudophosphatases
# stay in the catalytic master as untyped phosphatome members (Chen includes them); their biology is
# that of an adapter/activator of an ACTIVE myotubularin, recorded here as an annotation ON the gene
# rather than relocated to the regulatory-subunit companion (which holds only genes disjoint from the
# catalytic master). `regulates` is the pipe-joined target symbol(s); `regulatory_role` carries the
# role and its primary citation. Uncharacterised pseudophosphatases (MTMR10, MTMR11) are left blank.
PHOSPHATASE_REGULATORY_ROLES <- tibble::tribble(
  ~symbol,  ~regulates,             ~regulatory_role,
  "SBF1",   "MTMR2",                "activator of the active phosphatase MTMR2 (catalytically inactive partner; Kim et al. PNAS 2003, 10.1073/pnas.0431052100)",
  "SBF2",   "MTMR2",                "activator of MTMR2 (Robinson & Dixon, Berger et al. Hum Mol Genet 2006)",
  "MTMR12", "MTMR2",                "adapter subunit (3-PAP) of the MTMR2 lipid phosphatase (Nandurkar et al. PNAS 2003)",
  "MTMR9",  "MTMR6|MTMR7|MTMR8",    "activator/adapter of the MTMR6/7/8 lipid phosphatases (Zou et al. J Biol Chem 2009)")

classify_phosphatases <- function(universe_ensembl_ids, hgnc_bridge, go_phosphatase_sets, ec_phosphatase,
                                  membership, chen_facts, resolved_phosphatase,
                                  go_experimental_ids = character(0)) {

  # Add a TRUE/FALSE column named {{ flag }} marking rows whose ensembl_gene_id is in `ids`
  # (a tidy left join instead of a bare `%in%` membership test). `flag` is a bare column name.
  add_set_flag <- function(data, ids, flag) {
    data |>
      left_join(tibble(ensembl_gene_id = unique(ids), present_in_set = TRUE),
                by = join_by(ensembl_gene_id)) |>
      mutate("{{ flag }}" := coalesce(present_in_set, FALSE)) |>
      select(-present_in_set)
  }

  # Per-gene non-protein GO subtype: pipe-joined names of the nonprotein_by_subtype sublists whose
  # id-vector contains the gene (empty string if none). Drives nonprotein_substrate_type in
  # apply_term_sets() alongside the EC subtype.
  np_subtype_of <- function(id) {
    hits <- names(go_phosphatase_sets$nonprotein_by_subtype)[
      map_lgl(go_phosphatase_sets$nonprotein_by_subtype, \(ids) id %in% ids)]
    paste(hits, collapse = "|")
  }

  classified <- hgnc_bridge$gene_metadata |>
    semi_join(tibble(ensembl_gene_id = universe_ensembl_ids), by = join_by(ensembl_gene_id)) |>
    # EC table -- only the all_ec_codes list-column; apply_term_sets() recomputes EC typing from
    # the resolved term set, so the leg's precomputed protein/non-protein EC flags are not joined.
    left_join(ec_phosphatase$ec_table |> select(ensembl_gene_id, all_ec_codes),
              by = join_by(ensembl_gene_id)) |>
    # Chen per-gene facts (taxonomy + substrate flags + catalytic status), joined in.
    left_join(chen_facts, by = join_by(ensembl_gene_id)) |>
    # one binary membership column per source
    add_set_flag(membership$chen,              is_chen) |>
    add_set_flag(membership$hgnc_protein_group, is_hgnc_phosphatase_group) |>
    add_set_flag(membership$go_umbrella,       is_go_phosphatase_activity) |>
    add_set_flag(membership$ec,                is_phosphatase_ec) |>
    add_set_flag(membership$uniprot_keyword,   is_uniprot_kw_phosphatase) |>
    # GO functional sets used by the substrate flags.
    add_set_flag(go_phosphatase_sets$protein_phosphatase_activity, in_protein_phosphatase_go) |>
    add_set_flag(go_phosphatase_sets$nonprotein_all,               in_nonprotein_go) |>
    # Provenance proxy: non-electronic (experimental/curated) GO phosphatase-activity support.
    add_set_flag(go_experimental_ids,                              go_experimental) |>
    mutate(
      # An empty EC code list for genes the EC leg never saw, so map() over all_ec_codes is safe.
      all_ec_codes = map(all_ec_codes, \(codes) if (is.null(codes)) character(0) else codes),
      chen_nonprotein_substrate = coalesce(chen_nonprotein_substrate, FALSE),

      n_membership_sources = is_chen + is_hgnc_phosphatase_group + is_go_phosphatase_activity +
                             is_phosphatase_ec + is_uniprot_kw_phosphatase,

      # Axis 1 -- structural catalog: the Chen phosphatome OR an HGNC protein-phosphatase gene group.
      # The complementary Axis 2 (biochemical EC) is computed substrate-blind inside apply_term_sets().
      in_structural_catalog = is_chen | is_hgnc_phosphatase_group,
      # Supplementary support: experimental GO support OR the reviewed UniProt phosphatase keyword.
      # Neither is an axis; together they split Silver from Bronze among single-axis genes and cannot
      # manufacture standing from zero axes.
      has_uniprot_kw = is_uniprot_kw_phosphatase,
      supplementary_support = go_experimental | has_uniprot_kw,
      # Coarse proxy for experimental GO support of the substrate call. The substrate_decider field
      # it feeds is informational, not gating -- it records HOW the call was made, never changes it.
      go_experimental_protein = go_experimental,
      # catalytic_status from the Chen facts (active / pseudo / uncertain); everything else defaults
      # to active. A SOFT signal, never a veto.
      catalytic_status = coalesce(catalytic_status, "active"),

      # Per-gene evidence inputs for apply_term_sets(): the four co-equal substrate signals plus the
      # Chen curated non-protein flag.
      go_protein            = in_protein_phosphatase_go,
      go_nonprotein         = in_nonprotein_go,
      go_nonprotein_subtype = map_chr(ensembl_gene_id, np_subtype_of),
      chen_nonprotein       = chen_nonprotein_substrate)

  # The shared term-set path: EC axis flags (rigor substrate-blind), the four substrate flags, the
  # substrate call + provenance, the rigor tier, and the enrichment backgrounds. Binds its columns
  # back onto the classified table.
  classified <- apply_term_sets(classified, resolved_phosphatase)

  classified |>
    mutate(
      # BACK-COMPAT: the rest of the pipeline reads these names.
      is_protein_phosphatase_ec = ec_protein,
      # phosphatase_type: protein/dual -> protein; otherwise the firing non-protein substrate subtype,
      # in priority order; else a generic small-molecule label; else unclassified.
      phosphatase_type = case_when(
        substrate_call %in% c("protein", "dual")              ~ "Protein phosphatase",
        str_detect(nonprotein_substrate_type, "lipid")        ~ "Lipid phosphatase",
        str_detect(nonprotein_substrate_type, "nucleotide")   ~ "Nucleotide phosphatase",
        str_detect(nonprotein_substrate_type, "carbohydrate") ~ "Carbohydrate/sugar phosphatase",
        acts_on_nonprotein                                    ~ "Other small-molecule phosphatase",
        .default                                              = "Unclassified phosphatase"),
      classification_reason = case_when(
        substrate_call %in% c("protein", "dual") ~ str_c("protein phosphatase: substrate evidence ", substrate_evidence),
        acts_on_nonprotein                       ~ str_c(phosphatase_type, ": non-protein substrate (",
                                                         nonprotein_substrate_type, ")"),
        .default                                 = "in a phosphatase set but no protein or non-protein substrate evidence"),

      is_pseudophosphatase = coalesce(is_pseudophosphatase, FALSE),
      is_pseudogene = str_detect(coalesce(locus_type, ""), regex("pseudogene", ignore_case = TRUE)),
      # Chen taxonomy carried through to the master.
      phosphatase_fold      = chen_fold,
      phosphatase_family    = chen_family,
      phosphatase_subfamily = chen_subfamily,
      # membership_basis: the deriving source of the Axis-1 (structural-catalog) call, anchored on the
      # reconstructed Chen facts first, then the HGNC gene groups; NA when in no structural catalog.
      membership_basis = case_when(
        is_chen                   ~ "reconstructed:Chen2017",
        is_hgnc_phosphatase_group ~ "reconstructed:HGNC_groups",
        .default                  = NA_character_)) |>
    # Annotate the inactive-myotubularin adapters with their regulatory target + role.
    left_join(PHOSPHATASE_REGULATORY_ROLES, by = join_by(symbol)) |>
    mutate(regulates       = coalesce(regulates, NA_character_),
           regulatory_role = coalesce(regulatory_role, NA_character_)) |>
    # Final column set and order (drops all intermediate flag/gate columns).
    transmute(
      ensembl_gene_id,
      hgnc_symbol = symbol, hgnc_id, gene_name = name,
      acts_on_protein, acts_on_nonprotein, nonprotein_substrate_type,
      substrate_subtype = phosphatase_type,
      substrate_call, substrate_evidence, substrate_concordance, substrate_decider,
      ec_protein, ec_nonprotein, go_protein, go_nonprotein,
      dual_protein_nonprotein,
      catalytic_status, is_catalytic_background, is_protein_catalytic_background, is_pseudophosphatase,
      regulates, regulatory_role,
      n_evidence_dimensions, evidence_tier, curated_core,
      in_structural_catalog, is_protein_phosphatase_ec,
      go_experimental, has_uniprot_kw, supplementary_support, membership_basis,
      classification_reason,
      phosphatase_fold, phosphatase_family, phosphatase_subfamily,
      n_membership_sources, is_pseudogene,
      entrez_id, uniprot_ids, prev_symbol, alias_symbol,
      enzyme_id_EC = enzyme_id,
      hgnc_gene_group = gene_group, locus_type,
      chromosomal_location = location, mane_select_transcript = mane_select, iuphar_id = iuphar,
      is_chen, is_hgnc_phosphatase_group, is_go_phosphatase_activity, is_phosphatase_ec,
      is_uniprot_kw_phosphatase) |>
    arrange(desc(acts_on_protein), substrate_subtype, hgnc_symbol)
}
