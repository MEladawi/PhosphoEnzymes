# The four EC/GO term-set tables, externalized as cited data. Loads them, resolves the
# declarative rows into EC/GO matchers the gate uses, validates them against the pinned
# release, and applies them to per-gene evidence. One code path serves both the build
# gate and the accessor term_sets= override.

TERM_SET_FILES <- c(
  kinase_ec      = "kinase_ec_terms.csv",
  kinase_go      = "kinase_go_terms.csv",
  phosphatase_ec = "phosphatase_ec_terms.csv",
  phosphatase_go = "phosphatase_go_terms.csv")

# Read the four CSVs from a directory; record each file's md5.
load_term_sets <- function(extdata_dir) {
  paths <- set_names(file.path(extdata_dir, TERM_SET_FILES), names(TERM_SET_FILES))
  tables <- map(paths, \(p) read_csv(p, col_types = cols(.default = col_character())))
  list(tables = tables,
       md5    = set_names(unname(tools::md5sum(paths)), names(paths)))
}

# A gene's EC code reduced to its 3-field subclass ("2.7.10.1" -> "2.7.10"); NA if shorter.
ec_subclass_of <- function(code) {
  fields <- str_split(code, fixed("."))[[1]]
  if (length(fields) >= 3) paste(fields[1:3], collapse = ".") else NA_character_
}

# TRUE if any of a gene's EC codes matches any rule. `rules` is a tibble with `code` (a term_id
# such as "2.7.10.-" or "3.1.3.16") and `scope`. subclass rules match by 3-field prefix; exact
# rules match the full 4-field code.
matches_ec_rules <- function(gene_codes, rules) {
  if (length(gene_codes) == 0 || nrow(rules) == 0) return(FALSE)
  subclass_targets <- str_remove(rules$code[rules$scope == "subclass"], "\\.-$")
  exact_targets    <- rules$code[rules$scope == "exact"]
  gene_subclasses  <- unique(na.omit(map_chr(gene_codes, ec_subclass_of)))
  any(gene_subclasses %in% subclass_targets) || any(gene_codes %in% exact_targets)
}

# Resolve the loaded tables against a GMT into per-class matchers + GO id sets.
resolve_term_sets <- function(term_sets, go_gmt_path) {
  gmt <- read_gmt_accession_map(go_gmt_path)
  ids_for <- function(accessions) {
    present <- accessions[accessions %in% names(gmt)]
    unique(list_c(gmt[present]))
  }
  resolve_class <- function(ec_tbl, go_tbl) {
    ec_rigor <- filter(ec_tbl, role == "rigor+substrate")
    go_rigor <- filter(go_tbl, role == "rigor+substrate")
    go_np    <- filter(go_rigor, substrate == "nonprotein")
    subtypes <- unique(go_np$substrate_subtype[!is.na(go_np$substrate_subtype) & go_np$substrate_subtype != ""])
    list(
      ec_rigor       = select(ec_rigor, code = term_id, scope),
      ec_protein     = ec_rigor |> filter(substrate == "protein")    |> select(code = term_id, scope),
      ec_nonprotein  = ec_rigor |> filter(substrate == "nonprotein") |> select(code = term_id, scope, subtype = substrate_subtype),
      go_protein_ids    = ids_for(filter(go_rigor, substrate == "protein")$term_id),
      go_nonprotein_ids = ids_for(go_np$term_id),
      go_nonprotein_subtype_ids = set_names(
        map(subtypes, \(st) ids_for(go_np$term_id[go_np$substrate_subtype == st])), subtypes),
      go_umbrella_ids = ids_for(filter(go_tbl, role == "rigor_umbrella")$term_id),
      go_accessions   = go_tbl$term_id,
      ec_codes        = ec_tbl$term_id)
  }
  list(
    kinase      = resolve_class(term_sets$tables$kinase_ec,      term_sets$tables$kinase_go),
    phosphatase = resolve_class(term_sets$tables$phosphatase_ec, term_sets$tables$phosphatase_go),
    gmt_accessions = names(gmt),
    md5            = term_sets$md5)
}
