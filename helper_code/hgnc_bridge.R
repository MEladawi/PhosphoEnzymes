# helper_code/hgnc_bridge.R
# The HGNC complete set is the sole identifier authority. build_hgnc_bridge() reads it
# and returns the lookup maps, per-gene metadata, and a resolve_to_ensembl() closure
# that maps a source gene (by Entrez, UniProt, or symbol) to a base Ensembl gene ID.

build_hgnc_bridge <- function(hgnc_complete_set_path) {
  hgnc_complete_set <- read_tsv(hgnc_complete_set_path,
                                col_types = cols(.default = col_character()),
                                na = c("", "NA"), quote = "")
  # Fail at this boundary with a clear message if the upstream schema changed, rather than
  # erroring deep inside a transmute (the file is auto-refreshed, so columns can move/rename).
  required_columns <- c("status", "ensembl_gene_id", "hgnc_id", "symbol", "name",
                        "prev_symbol", "alias_symbol", "entrez_id", "uniprot_ids", "locus_type",
                        "gene_group", "enzyme_id", "iuphar", "mane_select", "location")
  missing_columns <- setdiff(required_columns, names(hgnc_complete_set))
  if (length(missing_columns))
    stop("HGNC complete set is missing expected column(s): ", paste(missing_columns, collapse = ", "),
         ". The source schema may have changed; check ", hgnc_complete_set_path, ".")
  approved_genes <- hgnc_complete_set %>% filter(status == "Approved")
  # Only rows that carry an Ensembl gene ID can be keyed.
  genes_with_ensembl <- approved_genes %>% filter(!is.na(ensembl_gene_id), ensembl_gene_id != "")

  # KEY -> ensembl maps, keeping the first occurrence of each key (file order).
  entrez_id_to_ensembl <- genes_with_ensembl %>%
    filter(!is.na(entrez_id), entrez_id != "") %>%
    distinct(entrez_id, .keep_all = TRUE) %>%
    { setNames(.$ensembl_gene_id, .$entrez_id) }

  uniprot_accession_to_ensembl <- genes_with_ensembl %>%
    transmute(ensembl_gene_id, uniprot_accession = map(uniprot_ids, split_pipe_delimited)) %>%
    unnest(uniprot_accession) %>%
    distinct(uniprot_accession, .keep_all = TRUE) %>%
    { setNames(.$ensembl_gene_id, .$uniprot_accession) }

  current_symbol_to_ensembl <- genes_with_ensembl %>%
    distinct(symbol, .keep_all = TRUE) %>%
    { setNames(.$ensembl_gene_id, .$symbol) }
  alias_and_previous_symbol_to_ensembl <- genes_with_ensembl %>%
    transmute(ensembl_gene_id,
              alias = map(alias_symbol, split_pipe_delimited),
              previous = map(prev_symbol, split_pipe_delimited)) %>%
    pivot_longer(c(alias, previous), values_to = "historical_symbol") %>%
    unnest(historical_symbol) %>%
    filter(!historical_symbol %in% names(current_symbol_to_ensembl)) %>%
    distinct(historical_symbol, .keep_all = TRUE) %>%
    { setNames(.$ensembl_gene_id, .$historical_symbol) }
  symbol_to_ensembl <- c(current_symbol_to_ensembl, alias_and_previous_symbol_to_ensembl)
  # Case-insensitive symbol index (last resort; resolves e.g. "SGK494" -> "SgK494"/RSKR).
  symbol_to_ensembl_caseinsensitive <- symbol_to_ensembl
  names(symbol_to_ensembl_caseinsensitive) <- toupper(names(symbol_to_ensembl))
  symbol_to_ensembl_caseinsensitive <-
    symbol_to_ensembl_caseinsensitive[!duplicated(names(symbol_to_ensembl_caseinsensitive))]

  # Locus type per Ensembl ID, used to reject implausible hits.
  ensembl_to_locus_type <- setNames(genes_with_ensembl$locus_type, genes_with_ensembl$ensembl_gene_id)
  ensembl_to_locus_type <- ensembl_to_locus_type[!duplicated(names(ensembl_to_locus_type))]
  # A kinase is a protein product (or pseudogene), never an RNA gene. A stale source ID
  # can point at a locus since reassigned to an ncRNA (e.g. a 2002 Entrez ID now belonging
  # to an antisense RNA), so candidates on an "RNA, ..." locus are skipped.
  candidate_is_rna_gene <- function(candidate_ensembl) {
    str_detect(lookup_in_named_vector(ensembl_to_locus_type, candidate_ensembl) %||% "",
               regex("RNA", ignore_case = TRUE))
  }

  # Every name a gene has ever carried (current + alias + previous), upper-cased. HGNC's
  # alias_symbol/prev_symbol ARE the symbol history, so this lets us compare a source's
  # (possibly outdated) symbol against every name the gene has had.
  symbol_history_long <- genes_with_ensembl %>%
    transmute(ensembl_gene_id,
              historical_symbol = pmap(list(symbol, alias_symbol, prev_symbol),
                function(current, aliases, previous)
                  toupper(c(current, split_pipe_delimited(aliases), split_pipe_delimited(previous))))) %>%
    unnest(historical_symbol) %>% distinct(ensembl_gene_id, historical_symbol)
  ensembl_to_symbol_history <- split(symbol_history_long$historical_symbol,
                                     symbol_history_long$ensembl_gene_id)
  # Does a candidate gene carry any of the source symbols (now or historically)?
  # No source symbol -> cannot disprove, treated as agreement.
  candidate_symbol_matches <- function(candidate_ensembl, source_symbols) {
    if (length(source_symbols) == 0) return(TRUE)
    known_symbols <- ensembl_to_symbol_history[[candidate_ensembl]]
    !is.null(known_symbols) && any(toupper(source_symbols) %in% known_symbols)
  }

  # Resolve by Entrez, then UniProt, then symbol/alias/previous (exact, then case-insensitive).
  # Among the hits, prefer one that is a plausible locus (not an RNA gene) AND whose symbol
  # history agrees with the source symbol -- this rejects stale Entrez/UniProt IDs that now
  # point at a different gene. Fall back to any non-RNA hit if none agrees.
  resolve_to_ensembl <- function(entrez_id = NA_character_,
                                 uniprot_accessions = character(0),
                                 source_symbols = character(0)) {
    source_symbols <- source_symbols[!is.na(source_symbols) & source_symbols != ""]
    candidate_ensembl_ids <- c(
      lookup_in_named_vector(entrez_id_to_ensembl, entrez_id %||% NA_character_),
      vapply(uniprot_accessions, function(accession) lookup_in_named_vector(uniprot_accession_to_ensembl, accession), character(1)),
      vapply(source_symbols,     function(symbol)    lookup_in_named_vector(symbol_to_ensembl, symbol), character(1)),
      vapply(source_symbols,     function(symbol)    lookup_in_named_vector(symbol_to_ensembl_caseinsensitive, toupper(symbol)), character(1)))
    candidate_ensembl_ids <- candidate_ensembl_ids[!is.na(candidate_ensembl_ids)]

    for (candidate in candidate_ensembl_ids)
      if (!candidate_is_rna_gene(candidate) && candidate_symbol_matches(candidate, source_symbols)) return(candidate)
    for (candidate in candidate_ensembl_ids)
      if (!candidate_is_rna_gene(candidate)) return(candidate)
    NA_character_
  }

  # Per-Ensembl metadata (first HGNC row per Ensembl ID).
  gene_metadata <- genes_with_ensembl %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    transmute(ensembl_gene_id, hgnc_id, symbol, name, prev_symbol, alias_symbol,
              entrez_id, uniprot_ids, locus_type, gene_group, enzyme_id,
              iuphar, mane_select, location)

  message(sprintf("  HGNC approved genes: %d | with Ensembl ID: %d",
                  nrow(approved_genes), nrow(genes_with_ensembl)))
  list(entrez_id_to_ensembl = entrez_id_to_ensembl,
       gene_metadata        = gene_metadata,
       resolve_to_ensembl   = resolve_to_ensembl)
}
