# helper_code/outputs.R
# Writing the master table (CSV + styled workbook), the Ensembl/symbol lists, the
# unmapped report, and the source-version manifest; plus a quality-control report with
# hard sanity-gene checks.

# Per-source gene counts, in display order.
count_genes_per_source <- function(kinases_table) {
  tibble(
    source = c("pkinfam","Manning (kinase.com)","KinHub","GO kinase activity",
               "HGNC EC 2.7","UniProt KW-0418","IDG dark kinome"),
    n = c(sum(kinases_table$is_pkinfam), sum(kinases_table$is_manning), sum(kinases_table$is_kinhub),
          sum(kinases_table$is_go_kinase_activity), sum(kinases_table$is_ec_kinase),
          sum(kinases_table$is_uniprot_kw_kinase), sum(kinases_table$is_idg_dark_kinase)))
}

write_outputs <- function(kinases_table, unmapped_table, output_dir, source_manifest, build_date) {
  write_csv(kinases_table, file.path(output_dir, "human_kinases_master.csv"))
  write_lines(kinases_table$ensembl_gene_id, file.path(output_dir, "kinases_ensembl_all.txt"))
  write_lines(kinases_table$ensembl_gene_id[kinases_table$protein_kinase],
              file.path(output_dir, "kinases_ensembl_protein.txt"))
  write_lines(kinases_table$ensembl_gene_id[kinases_table$confidence == "high"],
              file.path(output_dir, "kinases_ensembl_highconf.txt"))
  write_lines(kinases_table$hgnc_symbol, file.path(output_dir, "kinases_symbols_all.txt"))
  write_lines(kinases_table$hgnc_symbol[kinases_table$protein_kinase],
              file.path(output_dir, "kinases_symbols_protein.txt"))
  write_csv(distinct(unmapped_table), file.path(output_dir, "kinases_unmapped.csv"))
  write_tsv(source_manifest, file.path(output_dir, "source_versions.tsv"))

  type_counts   <- kinases_table %>% count(kinase_type, sort = TRUE)
  source_counts <- count_genes_per_source(kinases_table)

  workbook <- createWorkbook()
  addWorksheet(workbook, "Kinases")
  writeData(workbook, "Kinases", kinases_table, withFilter = TRUE)
  freezePane(workbook, "Kinases", firstActiveRow = 2)
  addStyle(workbook, "Kinases",
           createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                       fontColour = "white", fgFill = "#1F3864", halign = "left", border = "Bottom"),
           rows = 1, cols = seq_len(ncol(kinases_table)), gridExpand = TRUE)
  addStyle(workbook, "Kinases", createStyle(fontName = "Arial", fontSize = 10),
           rows = 2:(nrow(kinases_table) + 1), cols = seq_len(ncol(kinases_table)),
           gridExpand = TRUE, stack = TRUE)
  setColWidths(workbook, "Kinases", cols = seq_len(ncol(kinases_table)), widths = "auto")

  addWorksheet(workbook, "README")
  readme_lines <- c(
    paste("Human kinase reference table  |  build date:", build_date),
    "Key column: ensembl_gene_id (base, unversioned).",
    "",
    "SOURCES (version / fetched):",
    sprintf("  %-42s %s  [fetched %s]",
            source_manifest$source, source_manifest$version, source_manifest$fetched),
    "",
    "METHOD:",
    "  Membership = union of seven legs, each mapped to Ensembl through HGNC",
    "  (the GO leg is already Ensembl-keyed). Genes are typed by FUNCTION: a gene",
    "  with a non-protein GO kinase activity is typed non-protein ONLY IF it is not",
    "  in GO:0004672 (protein kinase activity). This gate classifies e.g.",
    "  PI4KA/SPHK1/DGKA as lipid kinases but PRKCA/ATM/MTOR/PIK3CA as protein",
    "  kinases. Non-protein priority: lipid > inositol-phosphate > carbohydrate >",
    "  nucleotide > creatine, then an EC-subclass fallback. EC kinase subclasses:",
    "  2.7.1-4,6,9-14 (2.7.5/2.7.7/2.7.8 excluded).",
    "",
    "COLUMN NOTES:",
    "  is_* ................ binary membership in each source (filter per source).",
    "  kinase_group ........ UniProt -> KinHub -> kinase.com (UniProt token kept only",
    "                        if it is a real Manning group; non-protein kinases get none).",
    "  kinase_family ....... Manning short tier (Akt, CDK ...): KinHub -> kinase.com.",
    "  kinase_subfamily .... UniProt -> KinHub -> kinase.com.",
    "  uniprot_protein_family : raw UniProt 'Protein families' string.",
    "  protein_kinase ...... TRUE for protein kinases; FALSE for small-molecule",
    "                        kinases (lipid / sugar / nucleotide / etc.).",
    "  curated_core ........ in a curated/functional leg (not GO-umbrella only);",
    "                        GO-only singletons hold most false positives.",
    "  confidence .......... high if >=2 sources or any EC kinase subclass.",
    "  dual_protein_and_nonprotein : protein kinase that also has a non-protein",
    "                        kinase function (e.g. PI3K family, NME).",
    "",
    "COUNTS BY kinase_type:",
    sprintf("  %-34s %d", type_counts$kinase_type, type_counts$n),
    "",
    "GENES PER SOURCE:",
    sprintf("  %-22s %d", source_counts$source, source_counts$n))
  writeData(workbook, "README", tibble(README = readme_lines))
  setColWidths(workbook, "README", cols = 1, widths = 95)
  saveWorkbook(workbook, file.path(output_dir, "human_kinases_master.xlsx"), overwrite = TRUE)
}

# Check the sanity genes (the cases the functional gate must get right) and, when
# verbose, print the QC summary. Returns TRUE if every sanity gene passed.
qc_report <- function(kinases_table, unmapped_table, verbose = TRUE) {
  emit <- function(...) if (verbose) cat(...)

  emit("\n==================== QC REPORT ====================\n")
  emit(sprintf("Total genes ............ %d\n", nrow(kinases_table)))
  emit(sprintf("Protein kinases ........ %d\n", sum(kinases_table$protein_kinase)))
  emit(sprintf("Non-protein kinases .... %d\n", sum(!kinases_table$protein_kinase)))
  emit(sprintf("High confidence ........ %d\n", sum(kinases_table$confidence == "high")))
  emit(sprintf("Curated core ........... %d\n", sum(kinases_table$curated_core)))
  emit(sprintf("Dual (protein+nonprot).. %d\n", sum(kinases_table$dual_protein_and_nonprotein)))
  emit(sprintf("Pseudogenes ............ %d\n", sum(kinases_table$is_pseudogene)))
  emit(sprintf("Unmapped (reported) .... %d\n", nrow(distinct(unmapped_table))))
  if (verbose) {
    cat("\nGenes per source:\n");      print(count_genes_per_source(kinases_table))
    cat("\nCounts by kinase_type:\n"); print(count(kinases_table, kinase_type, sort = TRUE), n = Inf)
  }

  sanity_genes <- tribble(
    ~symbol,  ~expected_type,                 ~expected_protein_kinase,
    "PI4KA",  "Lipid kinase",                 FALSE, "SPHK1",  "Lipid kinase", FALSE,
    "DGKA",   "Lipid kinase",                 FALSE, "HK1",    "Carbohydrate/sugar kinase", FALSE,
    "AK1",    "Nucleotide/nucleoside kinase", FALSE, "PIK3CA", "Protein kinase", TRUE,
    "NME1",   "Protein kinase",               TRUE,  "ATM",    "Protein kinase", TRUE,
    "MTOR",   "Protein kinase",               TRUE,  "PRKCA",  "Protein kinase", TRUE,
    "TRIB1",  "Protein kinase",               TRUE,  "CASK",   "Protein kinase", TRUE,
    "POMK",   "Protein kinase",               TRUE,  "FAM20C", "Protein kinase", TRUE)
  emit("\nSanity genes:\n"); all_passed <- TRUE
  for (i in seq_len(nrow(sanity_genes))) {
    gene_row <- kinases_table[toupper(kinases_table$hgnc_symbol) == toupper(sanity_genes$symbol[i]), ]
    passed <- nrow(gene_row) == 1 &&
              gene_row$kinase_type[1] == sanity_genes$expected_type[i] &&
              gene_row$protein_kinase[1] == sanity_genes$expected_protein_kinase[i]
    all_passed <- all_passed && passed
    emit(sprintf("  [%s] %-8s %-30s %s\n", if (passed) "PASS" else "FAIL", sanity_genes$symbol[i],
                 if (nrow(gene_row)) gene_row$kinase_type[1] else "ABSENT",
                 if (nrow(gene_row)) gene_row$protein_kinase[1] else ""))
  }
  emit(sprintf("\nSanity check: %s\n", if (all_passed) "ALL PASS" else "FAILURES ABOVE"))
  all_passed
}
