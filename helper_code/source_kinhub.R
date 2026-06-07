# helper_code/source_kinhub.R
# KinHub (kinhub.org) human kinase list, parsed from the live HTML table. Carries the
# Manning group/family/subfamily classification with current HGNC names. Mapped by
# UniProt accession, then HGNC symbol.

load_kinhub_kinome <- function(kinhub_html_path, hgnc_bridge) {
  kinhub_table <- read_html(kinhub_html_path) |> html_element("table") |> html_table()
  names(kinhub_table) <- str_squish(names(kinhub_table))           # normalise header whitespace/nbsp
  resolved <- kinhub_table |>
    transmute(symbol = `HGNC Name`, uniprot_accession = UniprotID, manning_name = `Manning Name`,
              group = Group, family = Family, subfamily = SubFamily) |>
    mutate(across(everything(), ~ na_if(str_trim(.x), "")),
           ensembl_gene_id = pmap_chr(list(uniprot_accession, symbol),
                                      ~ hgnc_bridge$resolve_to_ensembl(uniprot_accessions = ..1, source_symbols = ..2)))
  split_result <- split_mapped_and_unmapped(resolved, "KinHub", id_column = "uniprot_accession")
  list(ensembl_ids    = unique(split_result$mapped$ensembl_gene_id),
       taxonomy_table = split_result$mapped,   # ensembl_gene_id + group/family/subfamily
       unmapped       = split_result$unmapped)
}

# Taxonomy keyed by Ensembl, combined across sources by per-field priority.
# UniProt (actively curated, current) is primary for group and subfamily; the Manning
# intermediate "family" tier is taken from KinHub/kinase.com (UniProt has no equivalent
# tier). Any field is filled by the first source that has a value for that gene.
build_kinase_taxonomy <- function(uniprot_taxonomy_table, kinhub_taxonomy_table, manning_taxonomy_table) {
  first_row_per_gene <- function(taxonomy_table) taxonomy_table |> distinct(ensembl_gene_id, .keep_all = TRUE)
  uniprot <- first_row_per_gene(uniprot_taxonomy_table)
  kinhub  <- first_row_per_gene(kinhub_taxonomy_table)
  manning <- first_row_per_gene(manning_taxonomy_table)

  # Named (ensembl -> value) vector for one column, with blanks dropped.
  column_map <- function(taxonomy_table, column_name) {
    if (!column_name %in% names(taxonomy_table)) {
      setNames(character(0), character(0))
    } else {
      values <- setNames(taxonomy_table[[column_name]], taxonomy_table$ensembl_gene_id)
      values[!is.na(values) & values != ""]
    }
  }
  # First non-blank value per gene across the maps, in the order given (highest priority first).
  coalesce_maps <- function(...) {
    combined <- character(0)
    for (one_map in list(...)) combined <- c(combined, one_map[setdiff(names(one_map), names(combined))])
    combined
  }

  list(
    group     = coalesce_maps(column_map(uniprot, "group"),     column_map(kinhub, "group"),     column_map(manning, "group")),
    # kinase_family stays in the Manning short-label vocabulary (Akt, CDK, ...); the verbose
    # UniProt family string is NOT mixed in here -- it is kept separately in uniprot_protein_family.
    family    = coalesce_maps(column_map(kinhub,  "family"),    column_map(manning, "family")),
    subfamily = coalesce_maps(column_map(uniprot, "subfamily"), column_map(kinhub, "subfamily"), column_map(manning, "subfamily")),
    uniprot_family_raw = column_map(uniprot, "uniprot_protein_family"))
}
