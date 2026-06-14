# Rigor and substrate are independent axes.
#
# evidence_tier answers "how well is this gene's catalytic identity set?" and is
# SUBSTRATE-BLIND. The substrate columns answer a separate question. The
# cross-tabulation below makes the independence concrete: tier (rows) does not
# determine substrate (columns) -- some Gold-tier enzymes act on lipids.

library(PhosphoEnzymes)

p <- get_phosphatases()

cat("Phosphatase evidence_tier x substrate_call:\n")
print(with(p, table(evidence_tier = evidence_tier, substrate = substrate_call)))

# Concretely: Gold-tier (both rigor dimensions present) genes that are NOT
# protein-only -- their high rigor says nothing about what they act on.
cat("\nGold-tier phosphatases that act on a non-protein substrate:\n")
gold_np <- subset(
  p,
  evidence_tier == "Gold" & acts_on_nonprotein,
  select = c(symbol, substrate_call, nonprotein_substrate_type)
)
print(head(as.data.frame(gold_np), 10))

# The same point from the filter side: a strict (high-rigor) run still contains
# non-protein enzymes until you additionally filter on substrate.
cat("\nStrict kinases by substrate (rigor filter keeps non-protein enzymes):\n")
print(c(
  strict_any     = nrow(get_kinases(mode = "strict")),
  strict_nonprot = nrow(get_kinases(mode = "strict", substrate = "nonprotein"))
))
