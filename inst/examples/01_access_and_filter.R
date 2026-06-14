# Accessing the tables and the two orthogonal filters.
#
# PhosphoEnzymes ships three tables, one row per gene, keyed on base Ensembl ID:
#   get_kinases()        - the human kinase master
#   get_phosphatases()   - the human protein-phosphatase master
#   get_phosphoenzymes() - a thin class-agnostic summary spanning both
#
# get_kinases() / get_phosphatases() take two filters that never couple:
#   mode      - rigor:   "comprehensive" (default) or "strict" (>= 1 evidence dim)
#   substrate - acts-on: "any" (default), "protein", or "nonprotein"

library(PhosphoEnzymes)

k  <- get_kinases()
p  <- get_phosphatases()
pe <- get_phosphoenzymes()

cat("Table sizes (rows):\n")
print(c(kinases = nrow(k), phosphatases = nrow(p), phosphoenzymes = nrow(pe)))

# The two knobs are independent. "strict" filters rigor, not substrate: a
# strictly-curated lipid kinase is kept until you also ask for protein.
cat("\nKinase funnel across the two knobs:\n")
print(c(
  comprehensive = nrow(get_kinases()),
  strict = nrow(get_kinases(mode = "strict")),
  strict_protein = nrow(get_kinases(mode = "strict", substrate = "protein")),
  strict_nonprot = nrow(get_kinases(mode = "strict", substrate = "nonprotein"))
))

# The unified summary uses regulator_class to distinguish the two classes.
cat("\nUnified summary by class:\n")
print(table(pe$regulator_class))
