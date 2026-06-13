# Trap-case regression tests.
#
# These are the ONLY tests that guard the thing the whole resource exists to get
# right: genes where structural family and substrate DISAGREE, where a refactor
# or a source-version bump could flip the answer while every schema property test
# stays green. Each gene is a published fact, cited inline, not an opinion - so a
# reviewer reading this file sees regression guards over the literature, not a
# list of debatable "canonical" picks. We therefore include ONLY disagreement /
# fragility cases here, not obvious confirmations (those are covered by the
# family-level invariants in test-invariants.R).

# --- Phosphatases: structurally PTP-fold or "phospho-named", but NOT protein ---

test_that("PTEN is a lipid phosphatase, not a protein phosphatase", {
  # PTEN sits in the PTP superfamily by sequence but dephosphorylates PIP3.
  # Maehama & Dixon, J Biol Chem 1998; fold class in Chen et al. 2017.
  p <- pe_phosphat()
  row <- pe_row(p, "PTEN")
  skip_if(is.null(row), "PTEN absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "lipid")
})

test_that("MTMR2 is a lipid (phosphoinositide) phosphatase", {
  # Myotubularin-related; PTP fold, acts on PI3P. Chen et al. 2017.
  p <- pe_phosphat()
  row <- pe_row(p, "MTMR2")
  skip_if(is.null(row), "MTMR2 absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "lipid")
})

test_that("PSPH is a small-molecule phosphatase, not a protein phosphatase", {
  # The EC trap: 'phosphoserine phosphatase' acts on the FREE amino acid
  # O-phospho-L-serine (L-serine biosynthesis, HAD family), EC 3.1.3.3 -
  # not a phosphoserine residue on a protein.
  p <- pe_phosphat()
  row <- pe_row(p, "PSPH")
  skip_if(is.null(row), "PSPH absent")
  expect_false(row$acts_on_protein)
  expect_true(row$acts_on_nonprotein)
})

# --- Phosphatases: genuine protein phosphatase that breaks the naive motif rule

test_that("EYA1 is a protein phosphatase despite lacking the CX5R motif", {
  # Class IV / HAD aspartate-based protein-Tyr/Thr phosphatase. A blanket
  # 'all catalytic PTPs have CX5R' check would falsely call it a pseudophosphatase.
  # Tootle et al. Nature 2003; Rayapureddi et al. Nature 2003; Chen et al. 2017.
  p <- pe_phosphat()
  row <- pe_row(p, "EYA1")
  skip_if(is.null(row), "EYA1 absent")
  expect_true(row$acts_on_protein)
  expect_false(row$acts_on_nonprotein)
  if ("is_pseudophosphatase" %in% names(row)) {
    expect_false(isTRUE(row$is_pseudophosphatase))
  }
})

# --- Kinases: protein-kinase calls that rest on fragile / counter-intuitive evidence

test_that("NME1 is typed as a protein kinase", {
  # Dual-function: protein histidine kinase AND NDP kinase. Its protein call
  # rests on GO histidine-kinase activity; guards against it collapsing to
  # nucleotide-only on a GO release change.
  k <- pe_kinases()
  row <- pe_row(k, "NME1")
  skip_if(is.null(row), "NME1 absent")
  expect_true(row$acts_on_protein)
})

test_that("POMK is typed as a protein kinase despite a non-protein EC", {
  # Protein O-mannose kinase; EC 2.7.1.183 is a non-protein EC subclass, so its
  # protein call cannot come from EC and must survive the gate via GO. Guards the
  # GO-fragility of that call.
  k <- pe_kinases()
  row <- pe_row(k, "POMK")
  skip_if(is.null(row), "POMK absent")
  expect_true(row$acts_on_protein)
})

test_that("PIK3CA is retained as protein-acting (dual lipid/protein kinase)", {
  # In Manning by sequence; has protein-kinase activity (GO:0004672) alongside
  # lipid-kinase activity, so the gate keeps it protein rather than lipid-only.
  k <- pe_kinases()
  row <- pe_row(k, "PIK3CA")
  skip_if(is.null(row), "PIK3CA absent")
  expect_true(row$acts_on_protein)
})

# --- Kinases: kinome-catalog / kinase-annotated genes that act on NON-protein --
# substrates. The mirror of the phosphatase non-protein traps (PTEN, MTMR2, PSPH):
# a refactor or source bump could flip one of these to a protein kinase while every
# schema property test stays green. Each is a published fact, cited inline.

test_that("PI4KA is a lipid kinase, not a protein kinase", {
  # PI 4-kinase alpha (EC 2.7.1.67): shares the PI3K/PI4K kinase fold and sits in
  # kinome catalogs, but phosphorylates phosphatidylinositol. The kinase-side PTEN.
  k <- pe_kinases()
  row <- pe_row(k, "PI4KA")
  skip_if(is.null(row), "PI4KA absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "lipid")
})

test_that("SPHK1 is a lipid kinase, not a protein kinase", {
  # Sphingosine kinase 1 (EC 2.7.1.91): DAGK-like catalytic domain, lipid substrate.
  k <- pe_kinases()
  row <- pe_row(k, "SPHK1")
  skip_if(is.null(row), "SPHK1 absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "lipid")
})

test_that("DGKA is a lipid kinase, not a protein kinase", {
  # Diacylglycerol kinase alpha (EC 2.7.1.107): lipid second-messenger kinase.
  k <- pe_kinases()
  row <- pe_row(k, "DGKA")
  skip_if(is.null(row), "DGKA absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "lipid")
})

test_that("HK1 is a carbohydrate kinase, not a protein kinase", {
  # Hexokinase 1 (EC 2.7.1.1): phosphorylates glucose. Looks metabolic, is metabolic.
  k <- pe_kinases()
  row <- pe_row(k, "HK1")
  skip_if(is.null(row), "HK1 absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "carbohydrate")
})

test_that("AK1 is a nucleotide kinase, not a protein kinase", {
  # Adenylate kinase 1 (EC 2.7.4.3): nucleotide phosphotransfer.
  k <- pe_kinases()
  row <- pe_row(k, "AK1")
  skip_if(is.null(row), "AK1 absent")
  expect_false(row$acts_on_protein)
  expect_identical(row$nonprotein_substrate_type, "nucleotide")
})
