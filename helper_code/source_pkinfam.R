# helper_code/source_pkinfam.R
# UniProt pkinfam: the curated human protein kinome. Returns the Ensembl membership
# set and any rows that could not be mapped.

load_pkinfam_kinome <- function(pkinfam_path, hgnc_bridge) {
  file_lines <- read_lines(pkinfam_path, locale = locale(encoding = "latin1"))   # pkinfam is latin-1
  data_start_index <- grep("Swiss-Prot entries for protein kinases", file_lines)[1]
  if (is.na(data_start_index))
    stop("pkinfam: section marker 'Swiss-Prot entries for protein kinases' not found in ",
         pkinfam_path, " (the source layout may have changed).")
  # A data line looks like: "AKT1   AKT1_HUMAN (P31749 )  AKT1_MOUSE (P31750)".
  entry_pattern <- "^(\\S+)\\s+\\S+_HUMAN\\s*\\(\\s*([A-Z0-9]+)\\s*\\)"          # symbol + human accession

  data_lines <- if (data_start_index < length(file_lines))
    file_lines[(data_start_index + 1):length(file_lines)] else character(0)
  matched <- str_match(data_lines, entry_pattern)   # non-entry lines yield NA rows, dropped below

  resolved <- tibble(symbol = matched[, 2], uniprot_accession = matched[, 3]) %>%
    filter(!is.na(symbol)) %>%
    mutate(ensembl_gene_id = pmap_chr(list(uniprot_accession, symbol),
                                      ~ hgnc_bridge$resolve_to_ensembl(uniprot_accessions = ..1, source_symbols = ..2)))
  split_result <- split_mapped_and_unmapped(resolved, "pkinfam", id_column = "uniprot_accession")
  list(ensembl_ids = unique(split_result$mapped$ensembl_gene_id), unmapped = split_result$unmapped)
}
