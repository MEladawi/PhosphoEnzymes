# helper_code/go_functional_sets.R
# GO Molecular Function gene sets (Bader Lab EM_Genesets, Ensembl-native GMT). Members
# are already Ensembl gene IDs, so no Entrez bridge is needed for this leg. Sets are
# selected by stable GO accession (the 3rd "%"-field of each set name), which is robust
# to term-name wording. load_go_functional_sets() returns the functional umbrella (a
# membership leg), the protein-kinase discriminator, and the non-protein classes used
# to type non-protein kinases.

load_go_functional_sets <- function(go_gmt_path) {
  gmt_rows <- str_split(read_lines(go_gmt_path), fixed("\t"))
  gmt_rows <- gmt_rows[lengths(gmt_rows) >= 3]
  # Set-name field looks like "PROTEIN KINASE ACTIVITY%GOMF%GO:0004672"; take the GO id.
  go_accession_of_row <- map_chr(gmt_rows, function(row) {
    name_fields <- str_split(row[1], fixed("%"))[[1]]
    name_fields[length(name_fields)]
  })
  ensembl_ids_by_go_accession <- map(gmt_rows, ~ .x[-(1:2)])
  names(ensembl_ids_by_go_accession) <- go_accession_of_row

  # Union of the Ensembl members across one or more GO accessions; absent ones are skipped.
  ensembl_ids_for_go_accessions <- function(go_accessions) {
    present <- go_accessions[go_accessions %in% names(ensembl_ids_by_go_accession)]
    missing <- setdiff(go_accessions, present)
    if (length(missing))
      message("    [GO] not in this release, skipped: ", paste(missing, collapse = ", "))
    ensembl_ids_by_go_accession[present] |> list_c() |> unique()
  }

  go_sets <- list(
    kinase_activity_umbrella = ensembl_ids_for_go_accessions("GO:0016301"),  # KINASE ACTIVITY
    protein_kinase_activity  = ensembl_ids_for_go_accessions("GO:0004672"),  # PROTEIN KINASE ACTIVITY
    lipid_kinase = ensembl_ids_for_go_accessions(c(
      "GO:0001727",   # lipid kinase activity
      "GO:0052742",   # phosphatidylinositol kinase activity
      "GO:0016303",   # 1-phosphatidylinositol-3-kinase activity
      "GO:0046934",   # 1-phosphatidylinositol-4,5-bisphosphate 3-kinase activity
      "GO:0035005",   # 1-phosphatidylinositol-4-phosphate 3-kinase activity
      "GO:0016308",   # 1-phosphatidylinositol-4-phosphate 5-kinase activity
      "GO:0004143")), # ATP-dependent diacylglycerol kinase activity
    inositol_phosphate_kinase = ensembl_ids_for_go_accessions(c(
      "GO:0000828",   # inositol hexakisphosphate kinase activity
      "GO:0051766",   # inositol trisphosphate kinase activity
      "GO:0180030",   # inositol phosphate kinase activity
      "GO:0000827")), # inositol-1,3,4,5,6-pentakisphosphate kinase activity
    carbohydrate_kinase = ensembl_ids_for_go_accessions(c(
      "GO:0019200",   # carbohydrate kinase activity
      "GO:0004396",   # hexokinase activity
      "GO:0004340",   # glucokinase activity
      "GO:0008443")), # phosphofructokinase activity
    nucleotide_kinase = ensembl_ids_for_go_accessions(c(
      "GO:0019206",   # nucleoside kinase activity
      "GO:0050145",   # nucleoside monophosphate kinase activity
      "GO:0004550",   # nucleoside diphosphate kinase activity
      "GO:0019205",   # nucleobase-containing compound kinase activity
      "GO:0019136",   # deoxynucleoside kinase activity
      "GO:0047507",   # deoxynucleoside phosphate kinase activity, ATP as phosphate donor
      "GO:0036431",   # dCMP kinase activity
      "GO:0004385",   # GMP kinase activity
      "GO:0004017")), # AMP kinase activity
    creatine_kinase = ensembl_ids_for_go_accessions("GO:0004111"))  # creatine kinase activity

  assert_protein_kinase_set_is_propagated(go_sets$protein_kinase_activity)
  go_sets
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
