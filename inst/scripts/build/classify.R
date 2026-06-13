# helper_code/classify.R
# Per-gene classification by enzymatic function (the kinase gate) and assembly of the master
# table with one binary membership column per source. RIGOR is scored substrate-blind -- a gene
# stands in Axis 2 when it carries ANY class EC (protein or non-protein), the kinase EC subclasses
# -- and SUBSTRATE is typed by the four co-equal flags (GO protein / GO non-protein / EC protein /
# EC non-protein) with NO lineage default: a sequence-family catalog tells you the gene is a kinase,
# never which substrate it phosphorylates. The substrate flags + tiering both come from the single
# shared apply_term_sets() path, so the build gate and the accessor term_sets= override agree.
#
#   universe_ensembl_ids : sorted union of all membership legs
#   hgnc_bridge          : provides $gene_metadata (metadata also rides along in ec$ec_table)
#   go_sets              : load_go_functional_sets() output -- kinase_activity_umbrella,
#                          protein_kinase_activity, nonprotein_all, nonprotein_by_subtype
#   ec                   : load_ec_kinome() output; $ec_table carries gene metadata + all_ec_codes
#   membership           : named list of per-source Ensembl vectors
#   taxonomy             : build_kinase_taxonomy() output (group/family/subfamily maps)
#   resolved_kinase      : resolve_term_sets()$kinase -- the resolved EC/GO matchers

classify_kinases <- function(universe_ensembl_ids, hgnc_bridge, go_sets, ec, membership, taxonomy,
                             resolved_kinase,
                             go_experimental_ids = character(0),
                             pseudokinase_ensembl_ids = character(0)) {
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

  # Per-gene non-protein GO subtype: pipe-joined names of the nonprotein_by_subtype sublists whose
  # id-vector contains the gene (empty string if none). Drives nonprotein_substrate_type in
  # apply_term_sets() alongside the EC subtype.
  np_subtype_of <- function(id) {
    hits <- names(go_sets$nonprotein_by_subtype)[
      map_lgl(go_sets$nonprotein_by_subtype, \(ids) id %in% ids)]
    paste(hits, collapse = "|")
  }

  classified <- ec$ec_table |>
    semi_join(tibble(ensembl_gene_id = universe_ensembl_ids), by = join_by(ensembl_gene_id)) |>
    # one binary membership column per source
    add_set_flag(membership$pkinfam,         is_pkinfam) |>
    add_set_flag(membership$manning,         is_manning) |>
    add_set_flag(membership$kinhub,          is_kinhub) |>
    add_set_flag(membership$go_umbrella,     is_go_kinase_activity) |>
    # is_ec_kinase already exists as a column in ec$ec_table (== membership$ec), so reuse it.
    add_set_flag(membership$uniprot_keyword, is_uniprot_kw_kinase) |>
    add_set_flag(membership$idg_dark,        is_idg_dark_kinase) |>
    # GO functional sets used by the substrate flags.
    add_set_flag(go_sets$protein_kinase_activity, in_protein_kinase_go) |>
    add_set_flag(go_sets$nonprotein_all,          in_nonprotein_go) |>
    # Provenance proxy: non-electronic (experimental/curated) GO kinase-activity support.
    add_set_flag(go_experimental_ids,             go_experimental) |>
    mutate(
      n_membership_sources = is_pkinfam + is_manning + is_kinhub + is_go_kinase_activity +
                             is_ec_kinase + is_uniprot_kw_kinase + is_idg_dark_kinase,

      # Axis 1 -- structural/evolutionary sequence-family catalog: does the gene carry the kinase
      # sequence family? pkinfam / Manning / KinHub answer this from overlapping scholarship, so they
      # roll up into ONE dimension (counting them separately would inflate correlated agreement). The
      # complementary Axis 2 (biochemical EC) is computed substrate-blind inside apply_term_sets().
      in_structural_catalog = is_pkinfam | is_manning | is_kinhub,
      # Supplementary support: experimental GO support OR the reviewed UniProt kinase keyword. Neither
      # is an axis; together they split Silver from Bronze among single-axis genes and cannot
      # manufacture standing from zero axes.
      supplementary_support = go_experimental | is_uniprot_kw_kinase,
      # Coarse proxy for experimental GO support of the substrate call. The substrate_decider field
      # it feeds is informational, not gating -- it records HOW the call was made, never changes it.
      go_experimental_protein = go_experimental,
      # catalytic_status from the curated pseudokinase set (lineage TRUE, catalytically dead);
      # everything else defaults to active. A SOFT signal, never a veto.
      catalytic_status = if_else(ensembl_gene_id %in% pseudokinase_ensembl_ids, "pseudo", "active"),

      # Per-gene evidence inputs for apply_term_sets(): the four co-equal substrate signals.
      go_protein            = in_protein_kinase_go,
      go_nonprotein         = in_nonprotein_go,
      go_nonprotein_subtype = map_chr(ensembl_gene_id, np_subtype_of),
      chen_nonprotein       = FALSE)   # kinases have no Chen flag

  # The shared term-set path: EC axis flags (rigor substrate-blind), the four substrate flags, the
  # substrate call + provenance, the rigor tier, and the enrichment backgrounds. Binds its columns
  # back onto the classified table.
  classified <- apply_term_sets(classified, resolved_kinase)

  classified |>
    mutate(
      # BACK-COMPAT: the rest of the pipeline reads these names.
      protein_kinase        = acts_on_protein,
      is_protein_kinase_ec  = ec_protein,
      # kinase_type: protein/dual -> protein; otherwise the firing non-protein substrate subtype,
      # in priority order; else a generic small-molecule label; else unclassified.
      kinase_type = case_when(
        substrate_call %in% c("protein", "dual")              ~ "Protein kinase",
        str_detect(nonprotein_substrate_type, "lipid")        ~ "Lipid kinase",
        str_detect(nonprotein_substrate_type, "carbohydrate") ~ "Carbohydrate/sugar kinase",
        str_detect(nonprotein_substrate_type, "nucleotide")   ~ "Nucleotide/nucleoside kinase",
        str_detect(nonprotein_substrate_type, "metabolite")   ~ "Metabolite kinase",
        acts_on_nonprotein                                    ~ "Other small-molecule kinase",
        .default                                              = "Unclassified kinase"),
      classification_reason = case_when(
        substrate_call %in% c("protein", "dual") ~ str_c("protein kinase: substrate evidence ", substrate_evidence),
        acts_on_nonprotein                       ~ str_c(kinase_type, ": non-protein substrate (",
                                                         nonprotein_substrate_type, ")"),
        .default                                 = "in a kinase set but no protein or non-protein substrate evidence"),

      is_pseudogene = str_detect(coalesce(locus_type, ""), regex("pseudogene", ignore_case = TRUE)),
      ec_kinase_subclass = map_chr(matched_kinase_subclasses, \(x) paste(x, collapse = ", ")),
      # Manning taxonomy (named-vector maps keyed by Ensembl ID); NA where absent.
      kinase_group           = unname(taxonomy$group[ensembl_gene_id]),
      kinase_family          = unname(taxonomy$family[ensembl_gene_id]),
      kinase_subfamily       = unname(taxonomy$subfamily[ensembl_gene_id]),
      uniprot_protein_family = unname(taxonomy$uniprot_family_raw[ensembl_gene_id]),
      # Fallback family descriptor for genes with no Manning kinase_family (typically UniProt/GO-only
      # atypical kinases): the UniProt parsed family tier, else the firing non-protein substrate type.
      derived_family = if_else(
        !is.na(kinase_family) & kinase_family != "",
        NA_character_,
        coalesce(unname(taxonomy$uniprot_family_tier[ensembl_gene_id]),
                 na_if(nonprotein_substrate_type, ""))),
      # membership_basis: the deriving source of the Axis-1 (structural-catalog) call, anchored on the
      # cleanly-licensed leg first (pkinfam CC-BY), then the reconstructed Manning facts, then KinHub
      # as a cross-check; NA when the gene is in no structural catalog.
      membership_basis = case_when(
        is_pkinfam ~ "reconstructed:pkinfam",
        is_manning ~ "reconstructed:kinase.com",
        is_kinhub  ~ "crosscheck:KinHub",
        .default   = NA_character_),
      hgnc_kinase_gene_group = map_lgl(gene_group, \(group_field) {
        terms <- split_pipe_delimited(group_field)
        any(str_detect(terms, regex("kinase", ignore_case = TRUE)) & !str_detect(terms, noncatalytic_pattern))
      })) |>
    # Final column set and order (drops all intermediate flag/gate columns).
    transmute(
      ensembl_gene_id,
      hgnc_symbol = symbol, hgnc_id, gene_name = name,
      kinase_type, protein_kinase,
      acts_on_protein, acts_on_nonprotein, nonprotein_substrate_type,
      substrate_call, substrate_evidence, substrate_concordance, substrate_decider,
      ec_protein, ec_nonprotein, go_protein, go_nonprotein,
      catalytic_status, is_catalytic_background, is_protein_catalytic_background,
      kinase_group, kinase_family, derived_family, kinase_subfamily, uniprot_protein_family,
      dual_protein_nonprotein, evidence_tier, n_evidence_dimensions,
      go_experimental, supplementary_support,
      in_structural_catalog, is_protein_kinase_ec, classification_reason, membership_basis,
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
