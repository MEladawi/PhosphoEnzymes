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

# Split a pipe-delimited HGNC field into a trimmed character vector (drops blanks).
split_pipe_delimited <- function(field_value) {
  if (is.na(field_value) || field_value == "") return(character(0))
  parts <- str_trim(str_split(field_value, fixed("|"))[[1]])
  parts[parts != ""]
}

# Look up a value in a named vector; returns NA when the key is absent or blank.
# Uses single-bracket indexing (one keyed lookup) rather than an `%in% names()` scan.
lookup_in_named_vector <- function(named_vector, key) {
  if (is.na(key) || key == "") return(NA_character_)
  value <- unname(named_vector[key])
  if (is.na(value)) NA_character_ else value
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
