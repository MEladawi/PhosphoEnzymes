# helper_code/source_ec.R
# HGNC EC 2.7 phosphotransferase sweep. Reduces each EC number on a gene to its subclass
# (2.7.X) and flags kinase subclasses. 2.7.5 / 2.7.7 (nucleotidyltransferases) / 2.7.8 are
# NOT kinases. Operates directly on the HGNC metadata (no resolve step needed).

load_ec_kinome <- function(gene_metadata) {
  KINASE_EC_SUBCLASSES         <- c("2.7.1","2.7.2","2.7.3","2.7.4","2.7.6","2.7.9",
                                    "2.7.10","2.7.11","2.7.12","2.7.13","2.7.14")
  PROTEIN_KINASE_EC_SUBCLASSES <- c("2.7.10","2.7.11","2.7.12","2.7.13","2.7.14")

  ec_subclasses_of <- function(enzyme_id_field) {
    unique(na.omit(map_chr(split_pipe_delimited(enzyme_id_field), function(ec_number) {
      dotted_fields <- str_split(ec_number, fixed("."))[[1]]
      if (length(dotted_fields) >= 3) paste(dotted_fields[1:3], collapse = ".") else NA_character_
    })))
  }

  ec_table <- gene_metadata |> mutate(
    ec_subclasses          = map(enzyme_id, ec_subclasses_of),
    matched_kinase_subclasses = map(ec_subclasses, ~ intersect(.x, KINASE_EC_SUBCLASSES)),
    is_ec_kinase           = map_lgl(matched_kinase_subclasses, ~ length(.x) > 0),
    is_protein_kinase_ec   = map_lgl(ec_subclasses, ~ any(.x %in% PROTEIN_KINASE_EC_SUBCLASSES)))

  list(ensembl_ids = ec_table |> filter(is_ec_kinase) |> pull(ensembl_gene_id),
       ec_table    = ec_table)
}
