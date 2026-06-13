# GO Molecular Function gene sets for kinases and phosphatases. Both loaders now receive a
# resolved class object from resolve_term_sets() (term_sets.R) and delegate accession
# selection entirely to the curated CSV term sets. The GMT is ancestor-propagated, so each
# curated parent accession already carries its full subtree of gene members.

# Build the kinase GO sets from the resolved kinase term set (the GMT is ancestor-propagated,
# so each accession already carries its subtree). The protein set is asserted propagated.
load_go_functional_sets <- function(go_gmt_path, resolved_kinase) {
  go_sets <- list(
    kinase_activity_umbrella = resolved_kinase$go_umbrella_ids,
    protein_kinase_activity  = resolved_kinase$go_protein_ids,
    nonprotein_all           = resolved_kinase$go_nonprotein_ids,
    nonprotein_by_subtype    = resolved_kinase$go_nonprotein_subtype_ids)
  assert_protein_kinase_set_is_propagated(go_sets$protein_kinase_activity)
  go_sets
}

# Genes carrying NON-ELECTRONIC GO activity support, read from the no-IEA GMT variant (the
# upstream distribution with IEA / Inferred-from-Electronic-Annotation evidence stripped). This
# is the file-level proxy for experimental GO support: membership in a given accession's no-IEA
# set means the gene has at least one experimental/curated annotation to that term, not merely an
# electronic one. `go_accessions` selects which terms to union (default: the kinase umbrella
# GO:0016301 + protein-kinase GO:0004672, the supplementary support signal); pass a single
# protein-activity accession (GO:0004672 kinase, GO:0004721 phosphatase) to get the experimental
# protein set the substrate override keys on. Deliberately does NOT run the propagation assertion
# -- the no-IEA protein set is legitimately smaller than the IEA-inclusive one. Returns
# character(0) if the no-IEA file is absent, so the experimental signal is then FALSE everywhere.
load_go_experimental_ids <- function(no_iea_gmt_path,
                                     go_accessions = c("GO:0016301", "GO:0004672")) {
  if (!file.exists(no_iea_gmt_path)) {
    message("    [GO] no-IEA variant absent; experimental GO support = FALSE for all genes")
    return(character(0))
  }
  gmt_rows <- str_split(read_lines(no_iea_gmt_path), fixed("\t"))
  gmt_rows <- gmt_rows[lengths(gmt_rows) >= 3]
  go_accession_of_row <- map_chr(gmt_rows, function(row) {
    name_fields <- str_split(row[1], fixed("%"))[[1]]
    name_fields[length(name_fields)]
  })
  selected_rows <- go_accession_of_row %in% go_accessions
  gmt_rows[selected_rows] |> map(~ .x[-(1:2)]) |> list_c() |> unique()
}

# The whole functional gate rests on GO:0004672 (protein kinase activity) being
# ANCESTOR-PROPAGATED: the set must contain every gene annotated to any descendant term
# (Tyr / Ser-Thr / His ... kinase activity), not just genes annotated directly to the parent.
# A direct-only GMT would list only a handful of genes here and silently break the gate, so we
# fail the build loudly instead. Canaries are kinases annotated only to CHILD terms; their
# presence in GO:0004672 proves propagation. (Ensembl IDs are pinned with comments; if a future
# Ensembl release retires one, update it here.)
GO_PROTEIN_KINASE_CANARIES <- c(
  EGFR  = "ENSG00000146648",   # receptor tyrosine kinase  (child: GO:0004714)
  PRKCA = "ENSG00000154229",   # Ser/Thr protein kinase    (child: GO:0004674)
  NME1  = "ENSG00000239672")   # protein-histidine kinase  (child: GO:0004673)
assert_protein_kinase_set_is_propagated <- function(protein_kinase_members) {
  member_count <- length(protein_kinase_members)
  if (member_count < 300)
    stop("GO:0004672 (protein kinase activity) has only ", member_count, " members; expected ",
         "~500+. The GO MF GMT appears to carry DIRECT annotations only (not ancestor-",
         "propagated), which silently breaks the protein-kinase functional gate. Use an ",
         "ancestor-propagated GO MF GMT (the Bader Lab EM_Genesets distribution is propagated).")
  missing_canaries <- GO_PROTEIN_KINASE_CANARIES[!(GO_PROTEIN_KINASE_CANARIES %in% protein_kinase_members)]
  if (length(missing_canaries))
    stop("GO:0004672 is missing child-annotated canary protein kinase(s): ",
         paste(names(missing_canaries), collapse = ", "),
         ". The GMT may not be ancestor-propagated, or these Ensembl IDs were retired upstream.")
  message(sprintf("  [GO] protein kinase activity (GO:0004672): %d members (ancestor-propagated)",
                  member_count))
}

# Phosphatase GO sets from the resolved phosphatase term set. Protein set asserted propagated
# (incl. its reverse canary).
load_go_phosphatase_sets <- function(go_gmt_path, resolved_phosphatase) {
  go_sets <- list(
    phosphatase_activity_umbrella = resolved_phosphatase$go_umbrella_ids,
    protein_phosphatase_activity  = resolved_phosphatase$go_protein_ids,
    nonprotein_all                = resolved_phosphatase$go_nonprotein_ids,
    nonprotein_by_subtype         = resolved_phosphatase$go_nonprotein_subtype_ids)
  assert_protein_phosphatase_set_is_propagated(go_sets$protein_phosphatase_activity)
  go_sets
}

# GO:0004721 must be ancestor-propagated (children: Ser/Thr GO:0004722, Tyr GO:0004725, dual
# GO:0008138). Forward canaries are child-only protein phosphatases; the reverse canary (PSPH, a
# small-molecule phosphatase) must be ABSENT, guarding against a script that grants the protein
# child to every umbrella gene.
GO_PROTEIN_PHOSPHATASE_CANARIES <- c(
  PPP1CA = "ENSG00000172531",   # Ser/Thr protein phosphatase (child: GO:0004722)
  PTPN1  = "ENSG00000196396",   # protein tyrosine phosphatase (child: GO:0004725)
  DUSP1  = "ENSG00000120129")   # dual-specificity phosphatase (child: GO:0008138)
GO_PHOSPHATASE_REVERSE_CANARY <- c(PSPH = "ENSG00000146733")  # small-molecule; must NOT be present
assert_protein_phosphatase_set_is_propagated <- function(protein_phosphatase_members) {
  member_count <- length(protein_phosphatase_members)
  if (member_count < 100)
    stop("GO:0004721 (phosphoprotein phosphatase activity) has only ", member_count, " members; ",
         "expected ~150+. The GO MF GMT appears not to be ancestor-propagated, which silently ",
         "breaks the protein-phosphatase gate. Use the ancestor-propagated Bader Lab GMT.")
  missing_canaries <- GO_PROTEIN_PHOSPHATASE_CANARIES[!(GO_PROTEIN_PHOSPHATASE_CANARIES %in% protein_phosphatase_members)]
  if (length(missing_canaries))
    stop("GO:0004721 is missing child-annotated canary protein phosphatase(s): ",
         paste(names(missing_canaries), collapse = ", "), ". The GMT may not be propagated.")
  leaked <- GO_PHOSPHATASE_REVERSE_CANARY[GO_PHOSPHATASE_REVERSE_CANARY %in% protein_phosphatase_members]
  if (length(leaked))
    stop("Reverse canary failure: small-molecule phosphatase(s) ", paste(names(leaked), collapse = ", "),
         " appear in GO:0004721 (protein phosphatase activity). The GMT propagation is inverted or ",
         "contaminated -- a non-protein enzyme must not carry the protein-phosphatase child term.")
  message(sprintf("  [GO] protein phosphatase activity (GO:0004721): %d members (ancestor-propagated)",
                  member_count))
}
