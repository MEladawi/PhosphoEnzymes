# helper_code/utils.R
# Shared utilities used across the kinase build. Assumes the tidyverse packages are
# already attached by the calling script (build_kinases.R).

# Coalesce to a usable scalar: treat NULL / empty / NA / "" as missing. Length>1 values are
# returned as-is (the NA/"" emptiness test only applies to a single scalar).
`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0) {
    fallback
  } else if (length(value) == 1 && (is.na(value) || value == "")) {
    fallback
  } else {
    value
  }
}

# Split a pipe-delimited HGNC field into a trimmed character vector (drops blanks). Literal
# double-quotes are stripped first: some HGNC multi-value fields arrive quote-wrapped
# ("3.1.3.16|3.1.3.48|3.1.3.67"), which would otherwise leave a stray quote stuck to the first
# and last tokens and break exact matching (e.g. an EC code that no longer equals "3.1.3.16").
split_pipe_delimited <- function(field_value) {
  if (is.na(field_value) || field_value == "") {
    character(0)
  } else {
    parts <- str_trim(str_split(str_remove_all(field_value, '"'), fixed("|"))[[1]])
    parts[parts != ""]
  }
}

# Look up a value in a named vector; returns NA when the key is absent or blank.
# Uses single-bracket indexing (one keyed lookup) rather than an `%in% names()` scan.
lookup_in_named_vector <- function(named_vector, key) {
  if (is.na(key) || key == "") {
    NA_character_
  } else {
    value <- unname(named_vector[key])
    if (is.na(value)) NA_character_ else value
  }
}

# Split a resolved source table (must have `symbol` + `ensembl_gene_id` columns) into the
# rows that mapped and a tidy record of the rows that failed to map (carrying `id_column`
# as the source identifier). Nothing is dropped silently; callers collect $unmapped.
split_mapped_and_unmapped <- function(resolved_table, source_name, id_column = NULL) {
  failed_to_map <- is.na(resolved_table$ensembl_gene_id)
  source_ids <- if (is.null(id_column)) rep(NA_character_, sum(failed_to_map))
                else resolved_table[[id_column]][failed_to_map]
  list(
    mapped   = resolved_table[!failed_to_map, , drop = FALSE],
    unmapped = tibble(source = source_name,
                      symbol = resolved_table$symbol[failed_to_map],
                      id     = source_ids))
}

# Parse an Ensembl-keyed GMT into a named list of Ensembl-id vectors, keyed by the GO
# accession in each set-name's last "%"-field ("PROTEIN KINASE ACTIVITY%GOMF%GO:0004672").
# Rows with fewer than three tab-fields (no members) are dropped. Shared by the GO loader
# and the term-set resolver.
read_gmt_accession_map <- function(gmt_path) {
  rows <- str_split(read_lines(gmt_path), fixed("\t"))
  rows <- rows[lengths(rows) >= 3]
  accession_of <- map_chr(rows, \(row) {
    fields <- str_split(row[1], fixed("%"))[[1]]
    fields[length(fields)]
  })
  members <- map(rows, \(row) row[-(1:2)])
  set_names(members, accession_of)
}
