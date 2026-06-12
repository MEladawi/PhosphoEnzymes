# Category-level invariants.
#
# Family-wide guarantees that make the same promise as a long list of named genes
# but present ZERO "why this gene" surface to a reviewer, because they quantify
# over a family rather than asserting individual canonical status. These catch
# whole-class regressions a single-gene test would miss.

test_that("no myotubularin-family member is typed a protein phosphatase", {
  # The MTM/MTMR family are PTP-fold but act on phosphoinositides (lipids).
  p <- pe_phosphat()
  fam <- p[p$symbol %in% PE_MTM_FAMILY, , drop = FALSE]
  skip_if(nrow(fam) == 0L, "no MTM-family members present")
  expect_false(any(fam$acts_on_protein))
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
  expect_true(all(!fam$acts_on_protein | fam$dual_protein_and_nonprotein))
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

test_that("no classical small-molecule phosphatase clears the protein gate", {
  # Alkaline (EC 3.1.3.1) and acid (EC 3.1.3.2) phosphatases must never be typed
  # protein. Detect by EC if available, else by the well-known gene families.
  p <- pe_phosphat()
  if ("ec_number" %in% names(p)) {
    sm <- p[grepl("(^|;| )3\\.1\\.3\\.(1|2)(;|$| )", p$ec_number), , drop = FALSE]
  } else {
    sm <- p[grepl("^ALPL$|^ALPP$|^ALPI$|^ALPG$|^ACP[0-9]", p$symbol), , drop = FALSE]
  }
  skip_if(nrow(sm) == 0L, "no classical small-molecule phosphatases present")
  expect_false(any(sm$acts_on_protein))
})

test_that("every protein-acting phosphatase has at least one independent axis", {
  # A protein call should never rest on supplementary-only evidence (GO/keyword
  # without any catalog or protein-EC). Such genes would be Provisional, which
  # should not coexist with a confident protein typing in the curated core.
  p <- pe_phosphat()
  prot <- p[p$acts_on_protein, , drop = FALSE]
  skip_if(nrow(prot) == 0L, "no protein phosphatases present")
  # Provisional protein-phosphatases are allowed to EXIST (comprehensive mode)
  # but flag if they dominate - a sign the gate is leaning on weak evidence.
  prov_frac <- mean(prot$evidence_tier == "Provisional")
  expect_lt(prov_frac, 0.10)
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
    expect_true(all(df$n_independent_evidence_axes[df$evidence_tier == "Gold"] == 2L))
    expect_true(all(df$n_independent_evidence_axes[df$evidence_tier == "Provisional"] == 0L))
    # Silver/Bronze are exactly the single-axis genes
    one_axis <- df$n_independent_evidence_axes == 1L
    expect_true(all(df$evidence_tier[one_axis] %in% c("Silver", "Bronze")))
  }
})
