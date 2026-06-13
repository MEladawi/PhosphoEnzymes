# Phosphatase source legs. Mirror the kinase legs: each parses its source, resolves to a base
# Ensembl gene ID through the shared HGNC bridge, splits mapped/unmapped, and returns the
# membership set plus any taxonomy/per-gene tables the gate needs. Axis 1 (the structural
# catalog) is anchored on the Chen 2017 phosphatome (reconstructed facts), supplemented by the
# HGNC protein-phosphatase gene groups (CC0). Axis 2 is the surgical protein-phosphatase EC
# allow-list (3.1.3.16 / 3.1.3.48), drawn from HGNC enzyme_id and the UniProt EC column.

# --- Chen 2017 phosphatome (reconstructed facts) -----------------------------------------------
# The authoritative human phosphatome: membership in the 189-gene protein set, the
# fold/family/subfamily taxonomy, per-gene curated substrate flags (protein / non-protein), and
# catalytic-activity status (active / pseudophosphatase / uncertain). Resolved by Entrez ID then
# symbol.
load_chen_phosphatome <- function(chen_facts_path, hgnc_bridge) {
  # Resolve by Entrez ID then symbol (as the other legs do). A handful of Chen rows carry informal
  # names (Laforin, CIN) or pre-rename symbols (LPPR1, PPAPDC1A, TMEM55B) that resolve via neither
  # and flow into the unmapped record rather than being dropped silently; none are sanity genes.
  chen <- read_tsv(chen_facts_path, show_col_types = FALSE,
                   col_types = cols(entrez_id = col_character())) |>
    mutate(ensembl_gene_id = pmap_chr(list(entrez_id, symbol),
                                      ~ hgnc_bridge$resolve_to_ensembl(entrez_id = ..1,
                                                                       source_symbols = ..2)))
  split_result <- split_mapped_and_unmapped(chen, "Chen2017", id_column = "entrez_id")
  mapped <- split_result$mapped |> distinct(ensembl_gene_id, .keep_all = TRUE)
  list(
    # Axis-1 membership = the 189-gene curated protein phosphatome.
    ensembl_ids = mapped |> filter(included_in_phosphatome) |> pull(ensembl_gene_id) |> unique(),
    # All Chen rows (incl. the non-protein members) so the gate can read substrate flags + taxonomy.
    facts_table = mapped |>
      transmute(ensembl_gene_id,
                chen_fold = fold, chen_family = family, chen_subfamily = subfamily,
                chen_in_phosphatome = included_in_phosphatome,
                chen_protein_substrate    = protein_substrates    == "Yes",
                chen_nonprotein_substrate = nonprotein_substrates == "Yes",
                # active / pseudo / uncertain, mirroring the kinase catalytic_status enum.
                catalytic_status = case_when(
                  is_pseudophosphatase                      ~ "pseudo",
                  catalytic_activity == "Yes"               ~ "active",
                  .default                                  = "uncertain"),
                is_pseudophosphatase),
    unmapped = split_result$unmapped)
}

# --- HGNC phosphatase gene groups (CC0) --------------------------------------------------------
# HGNC organises phosphatases into curated gene groups. The protein-phosphatase catalytic groups
# supplement the Chen structural catalog (Axis-1 lineage); the regulatory/scaffold groups are the
# basis of the (deferred) regulatory companion and are returned separately so they never enter the
# catalytic master; the non-protein groups corroborate substrate typing. Operates directly on the
# HGNC metadata (already Ensembl-keyed; no resolve step).
HGNC_PROTEIN_PHOSPHATASE_GROUP <- regex(paste(
  "protein tyrosine phosphatase", "dual specificity phosphatase", "MAP kinase phosphatase",
  "protein phosphatase catalytic", "Mg2\\+/Mn2\\+ dependent", "CDC14", "CDC25",
  "CTD family phosphatase", "PTEN protein phosphatase", "EYA ", "Slingshot",
  "Serine/threonine phosphatase", "Cys-based", "HAD Asp-based protein", "LAR protein",
  "Protein tyrosine phosphatase 4A", sep = "|"), ignore_case = TRUE)
HGNC_REGULATORY_GROUP <- regex(paste(
  "regulatory subunit", "modulatory subunit", "scaffold subunit", "targeting", "actin regulator",
  sep = "|"), ignore_case = TRUE)

load_hgnc_phosphatase_groups <- function(gene_metadata) {
  group_terms <- function(group_field) split_pipe_delimited(group_field)
  has_group <- function(group_field, pattern)
    any(str_detect(group_terms(group_field), pattern))
  tagged <- gene_metadata |>
    mutate(
      in_protein_phosphatase_group = map_lgl(gene_group, ~ has_group(.x, HGNC_PROTEIN_PHOSPHATASE_GROUP)),
      in_regulatory_group          = map_lgl(gene_group, ~ has_group(.x, HGNC_REGULATORY_GROUP)),
      # any HGNC group mentioning "phosphatase" (broad umbrella incl. non-protein groups)
      in_any_phosphatase_group     = map_lgl(gene_group, ~ has_group(.x, regex("phosphatase", ignore_case = TRUE))))
  list(
    protein_phosphatase_ids = tagged |> filter(in_protein_phosphatase_group) |> pull(ensembl_gene_id),
    regulatory_ids          = tagged |> filter(in_regulatory_group, !in_protein_phosphatase_group) |> pull(ensembl_gene_id),
    any_phosphatase_ids     = tagged |> filter(in_any_phosphatase_group) |> pull(ensembl_gene_id))
}

# --- UniProt reviewed human, keyword KW-0904 ("Protein phosphatase") ---------------------------
# Supplementary membership/keyword signal, the actively-curated UniProt family taxonomy, and -- the
# reason this leg matters for Axis 2 -- the UniProt EC column, which assigns 3.1.3.48 to the EYA
# aspartate-based phosphatases that HGNC enzyme_id leaves blank. Mapped by accession, then symbol.
load_uniprot_keyword_phosphatome <- function(uniprot_tsv_path, hgnc_bridge) {
  resolved <- read_tsv(uniprot_tsv_path, col_types = cols(.default = col_character())) |>
    transmute(uniprot_accession = Entry, symbol = `Gene Names (primary)`,
              ec = `EC number`, protein_families = `Protein families`) |>
    mutate(ensembl_gene_id = pmap_chr(list(uniprot_accession, symbol),
                                      ~ hgnc_bridge$resolve_to_ensembl(uniprot_accessions = ..1, source_symbols = ..2)))
  split_result <- split_mapped_and_unmapped(resolved, "UniProt_KW0904", id_column = "uniprot_accession")
  mapped <- split_result$mapped |> distinct(ensembl_gene_id, .keep_all = TRUE)
  list(
    ensembl_ids = unique(mapped$ensembl_gene_id),
    # ensembl -> UniProt EC string (pipe/space-separated 4-digit codes), merged into the EC leg.
    ec_table = mapped |> filter(!is.na(ec) & ec != "") |> select(ensembl_gene_id, uniprot_ec = ec),
    taxonomy_table = mapped |> transmute(ensembl_gene_id, uniprot_protein_family = protein_families),
    unmapped = split_result$unmapped)
}

# --- Phosphatase EC leg (4-digit allow-list) ---------------------------------------------------
# Unlike kinases (whole 2.7.10-13 subclasses are protein), protein phosphatases are a few 4-digit
# needles inside the large, mostly non-protein EC 3.1.3 class. Axis 2 = 3.1.3.16 (Ser/Thr) or
# 3.1.3.48 (Tyr). Non-protein 4-digit 3.1.3.x (and 3.6.1.x nucleotide pyrophosphatases) type
# substrate but never score Axis 2. EC is read from HGNC enzyme_id, supplemented by the UniProt EC
# column (which carries EYA's 3.1.3.48).
PROTEIN_PHOSPHATASE_EC <- c("3.1.3.16", "3.1.3.48")
load_ec_phosphatome <- function(gene_metadata, uniprot_ec_table) {
  ec_long <- gene_metadata |>
    select(ensembl_gene_id, enzyme_id) |>
    left_join(uniprot_ec_table, by = join_by(ensembl_gene_id)) |>
    mutate(all_ec = map2(enzyme_id, uniprot_ec, function(hgnc_ec, up_ec) {
      codes <- c(split_pipe_delimited(hgnc_ec),
                 if (!is.na(up_ec)) str_trim(str_split(up_ec, "[;|]")[[1]]) else character(0))
      unique(codes[nzchar(codes) & codes != "-"])
    })) |>
    mutate(
      is_protein_phosphatase_ec = map_lgl(all_ec, ~ any(.x %in% PROTEIN_PHOSPHATASE_EC)),
      # non-protein 3.1.3.x / 3.6.1.x codes (fully specified 4-digit only -- a wildcard like
      # "3.1.3.-" names no substrate and must not trigger a non-protein call) minus the two
      # protein needles, for substrate typing.
      nonprotein_phosphatase_ec = map(all_ec, ~ setdiff(
        .x[str_detect(.x, "^3\\.1\\.3\\.[0-9]|^3\\.6\\.1\\.[0-9]")], PROTEIN_PHOSPHATASE_EC)))
  list(
    ensembl_ids = ec_long |> filter(map_lgl(all_ec, ~ any(str_detect(.x, "^3\\.1\\.3\\.")))) |> pull(ensembl_gene_id),
    ec_table = ec_long)
}
