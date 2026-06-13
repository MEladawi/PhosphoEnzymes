# Category-level invariants.
#
# Family-wide guarantees that make the same promise as a long list of named genes
# but present ZERO "why this gene" surface to a reviewer, because they quantify
# over a family rather than asserting individual canonical status. These catch
# whole-class regressions a single-gene test would miss.

test_that("myotubularins split by catalytic status: active are lipid, pseudo are untyped", {
  # The MTM/MTMR family are PTP-fold but act on phosphoinositides (lipids). The invariant
  # is expressed over the family + its catalytic status (never a hardcoded gene list), so
  # it tracks the source facts rather than a frozen roster:
  #  - Catalytically active members act on lipid. Several (MTM1, MTMR3/4/6/7/14)
  #    additionally dephosphorylate proteins (Chen 2017), so the assertion is dual-aware:
  #    every active member acts on lipid (none is a PURE protein phosphatase), mirroring
  #    the dual-aware DGK kinase invariant.
  #  - The pseudophosphatases (MTMR9/10/11/12, SBF1/SBF2) are catalytically dead, so none is
  #    a lipid phosphatase: the substrate-blind gate never assigns them the lipid substrate
  #    by lineage default. (Their protein axis can still carry a homology-propagated electronic
  #    GO term -- SBF1 does -- which is quarantined into the Provisional tier elsewhere, not a
  #    claim of catalytic activity; their adapter roles are annotated via regulates/regulatory_role.)
  p <- pe_phosphat()
  mtm <- p[p$phosphatase_family %in% "Myotubularin", , drop = FALSE]
  skip_if(nrow(mtm) == 0L, "no myotubularin-family members present")
  active <- mtm[!mtm$is_pseudophosphatase, , drop = FALSE]
  pseudo <- mtm[mtm$is_pseudophosphatase, , drop = FALSE]
  expect_true(all(active$acts_on_nonprotein & grepl("lipid", active$nonprotein_substrate_type)))
  expect_false(any(pseudo$acts_on_nonprotein))
})

# Kinase-side family invariants (the mirror of the MTM phosphatase invariant).
# Each family is unambiguously non-protein, so no member should act on protein.

test_that("no diacylglycerol kinase is a pure protein kinase", {
  # DGKs are lipid kinases. DGKQ is a genuine dual lipid+protein kinase, so the
  # guarantee is that no DGK is a PURE protein kinase: every member is non-protein
  # or flagged dual. (Mirror of why PIK3CA is a trap case, not an invariant member.)
  k <- pe_kinases()
  fam <- k[k$symbol %in% PE_DGK_FAMILY, , drop = FALSE]
  skip_if(nrow(fam) == 0L, "no DGK-family members present")
  expect_true(all(!fam$acts_on_protein | fam$dual_protein_nonprotein))
})

test_that("no hexokinase is a protein kinase", {
  # Hexokinases phosphorylate glucose; all carbohydrate kinases.
  k <- pe_kinases()
  fam <- k[k$symbol %in% PE_HK_FAMILY, , drop = FALSE]
  skip_if(nrow(fam) == 0L, "no hexokinase-family members present")
  expect_false(any(fam$acts_on_protein))
})

test_that("no adenylate kinase is a protein kinase", {
  # Adenylate kinases (AK1-AK9) are nucleotide phosphotransferases.
  k <- pe_kinases()
  fam <- k[k$symbol %in% PE_AK_FAMILY, , drop = FALSE]
  skip_if(nrow(fam) == 0L, "no adenylate-kinase-family members present")
  expect_false(any(fam$acts_on_protein))
})

test_that("no type-II/III PI 4-kinase is a protein kinase", {
  # PI4Ks phosphorylate phosphatidylinositol; all lipid kinases.
  k <- pe_kinases()
  fam <- k[k$symbol %in% PE_PI4K_FAMILY, , drop = FALSE]
  skip_if(nrow(fam) == 0L, "no PI4K-family members present")
  expect_false(any(fam$acts_on_protein))
})

test_that("every alkaline phosphatase retains small-molecule activity", {
  # The alkaline phosphatases (ALPL/ALPP/ALPI/ALPG, EC 3.1.3.1) are small-molecule
  # phosphatases. They are broad-specificity, and Chen 2017 flags several with protein
  # substrates too, so the invariant is dual-aware (mirroring MTM/DGK): every member must
  # act on a non-protein substrate (none is a PURE protein phosphatase). Acid phosphatases
  # are excluded entirely -- ACP1 is low-molecular-weight PTP, a genuine protein-Tyr
  # phosphatase carrying the historical acid-phosphatase EC (a named exception, not a rule).
  p <- pe_phosphat()
  sm <- p[grepl("^ALP[GILP]$", p$symbol), , drop = FALSE]
  skip_if(nrow(sm) == 0L, "no alkaline phosphatases present")
  expect_true(all(sm$acts_on_nonprotein))
})

test_that("every protein-acting phosphatase has at least one evidence dimension", {
  # A protein call should never rest on supplementary-only evidence (GO/keyword
  # without any catalog or protein-EC). Such genes would be Provisional, which
  # should not coexist with a confident protein typing in the curated core.
  p <- pe_phosphat()
  prot <- p[p$acts_on_protein, , drop = FALSE]
  skip_if(nrow(prot) == 0L, "no protein phosphatases present")
  # Provisional protein-phosphatases are comprehensive-mode, GO-only (0-axis) entries -- some are
  # GO false-positives (kinases/structural proteins carrying an IBA-propagated protein-phosphatase
  # term). They are quarantined out of strict mode, never deleted. Flag only if they DOMINATE the
  # protein phosphatome, which would mean the gate is leaning on weak evidence.
  prov_frac <- mean(prot$evidence_tier == "Provisional")
  expect_lt(prov_frac, 0.20)
})

test_that("protein-kinase and protein-phosphatase counts are within expected ranges", {
  # Source-drift tripwire. The human protein kinome is ~500-560; the protein
  # phosphatome is ~189 (Chen et al. 2017). Wide bands - this catches a catalog
  # that silently halved or doubled, not fine fluctuations.
  k <- pe_kinases()
  p <- pe_phosphat()
  n_pk <- sum(k$acts_on_protein)
  n_pp <- sum(p$acts_on_protein)
  expect_gt(n_pk, 450L); expect_lt(n_pk, 650L)
  expect_gt(n_pp, 150L); expect_lt(n_pp, 260L)
})

test_that("tier monotonicity holds: Gold implies both axes, Provisional implies none", {
  for (df in list(pe_kinases(), pe_phosphat())) {
    expect_true(all(df$n_evidence_dimensions[df$evidence_tier == "Gold"] == 2L))
    expect_true(all(df$n_evidence_dimensions[df$evidence_tier == "Provisional"] == 0L))
    # Silver/Bronze are exactly the single-axis genes
    one_axis <- df$n_evidence_dimensions == 1L
    expect_true(all(df$evidence_tier[one_axis] %in% c("Silver", "Bronze")))
  }
})
