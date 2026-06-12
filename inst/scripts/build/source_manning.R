# helper_code/source_manning.R
# Manning human kinome (kinase.com, 2002). Provides the group/family/subfamily taxonomy.
# Pseudogene entries are excluded. Mapped by Entrez GeneID, then symbol/Manning name.

load_manning_kinome <- function(manning_xls_path, hgnc_bridge) {
  manning_table <- read_excel(manning_xls_path, sheet = 1, col_types = "text") |>
    filter(is.na(.data[["Pseudogene?"]]) | .data[["Pseudogene?"]] != "Y")
  resolved <- manning_table |>
    transmute(manning_name = Name, group = Group, family = Family, subfamily = Subfamily,
              entrez_id = Entrez_GeneID, symbol = Entrez_Symbol) |>
    mutate(ensembl_gene_id = pmap_chr(list(entrez_id, symbol, manning_name),
                                      ~ hgnc_bridge$resolve_to_ensembl(entrez_id = ..1,
                                                                       source_symbols = c(..2, ..3))))
  split_result <- split_mapped_and_unmapped(resolved, "Manning", id_column = "entrez_id")
  list(ensembl_ids   = unique(split_result$mapped$ensembl_gene_id),
       taxonomy_table = split_result$mapped,   # ensembl_gene_id + group/family/subfamily
       unmapped       = split_result$unmapped)
}
