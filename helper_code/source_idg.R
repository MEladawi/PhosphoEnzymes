# helper_code/source_idg.R
# IDG understudied ("dark") kinome: a curated list of HGNC symbols (classic-Mac CR
# line endings). Mapped by symbol.

load_idg_dark_kinome <- function(idg_csv_path, hgnc_bridge) {
  symbols <- str_trim(unlist(str_split(read_lines(idg_csv_path), "\r")))
  symbols <- symbols[symbols != "" & symbols != "Kinase_HGNC_ID"]
  resolved <- tibble(symbol = symbols) |>
    mutate(ensembl_gene_id = map_chr(symbol, ~ hgnc_bridge$resolve_to_ensembl(source_symbols = .x)))
  split_result <- split_mapped_and_unmapped(resolved, "IDG_dark")
  list(ensembl_ids = unique(split_result$mapped$ensembl_gene_id), unmapped = split_result$unmapped)
}
