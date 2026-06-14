# Recovering class taxonomy from the unified summary.
#
# get_phosphoenzymes() keeps only the columns shared by both classes, so
# class-specific taxonomy (kinase group/family, phosphatase fold/family) stays
# on the masters. Recover it with a CLASS-ROUTED JOIN: split the summary by
# regulator_class and join each part back to the matching master on
# ensembl_gene_id.

library(PhosphoEnzymes)

pe <- get_phosphoenzymes()

# Kinase rows -> back to the kinase master for kinase_family.
kin_tax <- merge(
  pe[pe$regulator_class == "kinase", c("ensembl_gene_id", "evidence_tier")],
  get_kinases()[, c("ensembl_gene_id", "symbol", "kinase_family")],
  by = "ensembl_gene_id"
)
cat("Kinase summary rows joined to kinase_family:\n")
print(head(kin_tax[, c("symbol", "evidence_tier", "kinase_family")], 5))

# Phosphatase rows -> back to the phosphatase master for phosphatase_fold.
phos <- pe[pe$regulator_class == "phosphatase",
           c("ensembl_gene_id", "substrate_call")]
phos_tax <- merge(
  phos,
  get_phosphatases()[, c("ensembl_gene_id", "symbol", "phosphatase_fold")],
  by = "ensembl_gene_id"
)
cat("\nPhosphatase summary rows joined to phosphatase_fold:\n")
print(head(phos_tax[, c("symbol", "substrate_call", "phosphatase_fold")], 5))
