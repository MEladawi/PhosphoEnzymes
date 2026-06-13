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

# Structural + provenance validation. Returns a tibble(severity, table, term_id, message).
# `severity == "error"` rows are fatal for the DEFAULT set (the caller stops); for a user set
# they downgrade to warnings at the call site. When `resolved` is supplied (built against a
# GMT), term-id resolution and the reverse canaries are also checked.
validate_term_set <- function(term_sets, resolved = NULL,
                              reverse_canaries = TERM_SET_REVERSE_CANARIES) {
  issues <- list()
  add <- function(severity, table, term_id, message)
    issues[[length(issues) + 1]] <<- tibble(severity = severity, table = table,
                                            term_id = term_id, message = message)

  for (nm in names(term_sets$tables)) {
    tbl <- term_sets$tables[[nm]]
    need <- c("term_id", "class", "substrate", "substrate_subtype", "role", "scope", "citation", "note")
    miss <- setdiff(need, names(tbl))
    if (length(miss)) add("error", nm, NA, paste("missing columns:", paste(miss, collapse = ",")))
    no_cite <- tbl$term_id[is.na(tbl$citation) | tbl$citation == ""]
    for (t in no_cite) add("error", nm, t, "missing citation")
    rig <- filter(tbl, role == "rigor+substrate")
    bad_np <- rig$term_id[rig$substrate == "nonprotein" &
                            (is.na(rig$substrate_subtype) | rig$substrate_subtype == "")]
    for (t in bad_np) add("error", nm, t, "nonprotein rigor row lacks substrate_subtype")
    umb <- filter(tbl, role == "rigor_umbrella")
    bad_umb <- umb$term_id[umb$substrate != "na"]
    for (t in bad_umb) add("error", nm, t, "umbrella row must have substrate == na")
  }
  for (nm in names(term_sets$tables)) {
    tbl <- filter(term_sets$tables[[nm]], role == "rigor+substrate")
    dup <- tbl |> summarise(n_sub = n_distinct(substrate), .by = term_id) |> filter(n_sub > 1)
    for (t in dup$term_id) add("error", nm, t, "substrate overlap: term_id tagged both protein and nonprotein")
  }
  if (!is.null(resolved)) {
    for (cls in c("kinase", "phosphatase")) {
      go_tbl <- term_sets$tables[[paste0(cls, "_go")]]
      unresolved <- setdiff(go_tbl$term_id, resolved$gmt_accessions)
      for (t in unresolved) add("error", paste0(cls, "_go"), t, "GO term_id not present in the pinned GMT release")
    }
    for (canary in reverse_canaries) {
      hits <- resolved[[canary$class]]$go_protein_ids
      if (canary$ensembl %in% hits)
        add("error", paste0(canary$class, "_go"), canary$symbol,
            paste0("reverse canary: ", canary$symbol, " must have go_protein == FALSE"))
    }
  }
  if (length(issues)) list_rbind(issues) else
    tibble(severity = character(), table = character(), term_id = character(), message = character())
}

# Pure non-protein enzymes that must never carry a protein-activity GO leaf (the guard against
# EC-derived electronic GO leaking a non-protein enzyme into the protein set). Ensembl IDs are
# pinned; update if an Ensembl release retires one.
TERM_SET_REVERSE_CANARIES <- list(
  list(class = "kinase",      symbol = "PI4KA", ensembl = "ENSG00000241973"),
  list(class = "phosphatase", symbol = "PSPH",  ensembl = "ENSG00000146733"))

# Per-gene EC axis flags against a resolved class. Returns the rigor flag (E), the protein and
# nonprotein flags, and the firing nonprotein substrate subtype(s).
ec_axis_flags <- function(gene_codes, resolved_class) {
  np_rows <- resolved_class$ec_nonprotein
  fired_np <- if (length(gene_codes) == 0 || nrow(np_rows) == 0) character(0) else {
    np_rows$subtype[map_lgl(seq_len(nrow(np_rows)),
      \(i) matches_ec_rules(gene_codes, np_rows[i, c("code", "scope")]))]
  }
  list(
    ec_rigor      = matches_ec_rules(gene_codes, resolved_class$ec_rigor),
    ec_protein    = matches_ec_rules(gene_codes, resolved_class$ec_protein),
    ec_nonprotein = matches_ec_rules(gene_codes, resolved_class$ec_nonprotein),
    nonprotein_subtypes = unique(fired_np[!is.na(fired_np) & fired_np != ""]))
}
