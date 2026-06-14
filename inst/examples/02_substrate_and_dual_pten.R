# Substrate typing and the dual enzyme PTEN.
#
# Substrate is carried as two CO-EQUAL booleans, not one label:
#   acts_on_protein, acts_on_nonprotein
# summarised by substrate_call (protein / nonprotein / dual / untyped), with
# nonprotein_substrate_type naming the chemical class. A bifunctional enzyme is
# TRUE for both flags and stays "dual" rather than being collapsed to one label.

library(PhosphoEnzymes)

k <- get_kinases()

# The two flags are not mutually exclusive: the off-diagonal "dual" cell is real.
cat("Kinase substrate flags (acts_on_protein x acts_on_nonprotein):\n")
print(table(
  acts_on_protein    = k$acts_on_protein,
  acts_on_nonprotein = k$acts_on_nonprotein
))

# PTEN is the canonical dual case: a Gold-tier phosphatase acting on BOTH a
# protein substrate and the lipid PIP3.
cat("\nPTEN, straight from the shipped table:\n")
pten <- subset(
  get_phosphatases(),
  symbol == "PTEN",
  select = c(symbol, evidence_tier, acts_on_protein, acts_on_nonprotein,
             substrate_call, nonprotein_substrate_type)
)
print(as.data.frame(pten))

# Every substrate call is auditable: which evidence fired, which source decided,
# and whether the available evidence agreed.
cat("\nPTEN substrate provenance:\n")
print(as.data.frame(subset(
  get_phosphatases(),
  symbol == "PTEN",
  select = c(substrate_evidence, substrate_decider, substrate_concordance)
)))
