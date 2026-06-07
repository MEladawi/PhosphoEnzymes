# helper_code/source_uniprot_kw.R
# UniProtKB reviewed human proteins carrying keyword KW-0418 ("Kinase"). Broad,
# all-kinase-types functional annotation. Mapped by UniProt accession, then symbol.
# Also yields the actively-curated UniProt protein-family classification (group /
# family / subfamily), used as the primary, up-to-date taxonomy source.

# Manning kinase groups. The leading token of a UniProt protein-kinase family string is
# the Manning group ONLY for these codes; anything else (e.g. "Ser/Thr" from a prefix-less
# "Ser/Thr protein kinase family", or "Hexokinase"/"Adenylate" from a non-protein-kinase
# family) is not a group and must not populate kinase_group.
MANNING_GROUP_CODES <- c("AGC","CAMK","CK1","CMGC","NEK","RGC","STE","TK","TKL","Other","Atypical")

# Parse a UniProt "Protein families" string, e.g.
#   "Protein kinase superfamily, AGC Ser/Thr protein kinase family, RAC subfamily"
# into group ("AGC"), family ("AGC Ser/Thr protein kinase family"), subfamily ("RAC").
parse_uniprot_protein_family <- function(families_string) {
  if (is.na(families_string) || families_string == "") {
    list(group = NA_character_, family = NA_character_, subfamily = NA_character_)
  } else {
    parts <- str_trim(str_split(families_string, ",")[[1]])
    parts <- parts[parts != ""]
    is_superfamily <- str_detect(parts, regex("superfamily$", ignore_case = TRUE))
    is_subfamily   <- str_detect(parts, regex("subfamily$",   ignore_case = TRUE))
    is_family      <- str_detect(parts, regex("family$",      ignore_case = TRUE)) & !is_superfamily & !is_subfamily

    family    <- if (any(is_family))    parts[is_family][1] else NA_character_
    subfamily <- if (any(is_subfamily)) str_trim(str_remove(parts[is_subfamily][1],
                                         regex("\\s*subfamily$", ignore_case = TRUE))) else NA_character_
    group <- if (!is.na(family)) str_split(family, "\\s+")[[1]][1] else NA_character_  # leading token of the family
    if (!is.na(group) && group == "Tyr") group <- "TK"                                # UniProt label -> Manning label
    if (!is.na(group) && !(group %in% MANNING_GROUP_CODES)) group <- NA_character_     # only real Manning groups
    list(group = group, family = family, subfamily = subfamily)
  }
}

load_uniprot_keyword_kinome <- function(uniprot_tsv_path, hgnc_bridge) {
  resolved <- read_tsv(uniprot_tsv_path, col_types = cols(.default = col_character())) |>
    transmute(uniprot_accession = Entry, symbol = `Gene Names (primary)`,
              protein_families = `Protein families`) |>
    mutate(ensembl_gene_id = pmap_chr(list(uniprot_accession, symbol),
                                      ~ hgnc_bridge$resolve_to_ensembl(uniprot_accessions = ..1, source_symbols = ..2)))
  split_result <- split_mapped_and_unmapped(resolved, "UniProt_KW0418", id_column = "uniprot_accession")

  parsed <- map(split_result$mapped$protein_families, parse_uniprot_protein_family)
  taxonomy_table <- split_result$mapped |>
    transmute(ensembl_gene_id,
              group     = map_chr(parsed, "group"),
              family    = map_chr(parsed, "family"),
              subfamily = map_chr(parsed, "subfamily"),
              uniprot_protein_family = protein_families)

  list(ensembl_ids    = unique(split_result$mapped$ensembl_gene_id),
       taxonomy_table = taxonomy_table,
       unmapped       = split_result$unmapped)
}
