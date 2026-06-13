# Substrate-typing guards for the phosphoinositide-phosphatase class.
#
# The lipid subtype for this class is reached through one broad GO parent
# (GO:0052866, phosphatidylinositol phosphate phosphatase activity) whose
# ancestor-propagated subtree is the whole class. A broad subtree parent is only
# safe while it stays bounded, so these tests pin its edges: the class types as
# lipid, the one soluble-inositol exception does not, genuine duals are not
# flattened to lipid-only, and the dead members carry no substrate.

# The active phosphoinositide phosphatases the GO:0052866 parent must reach.
PE_PHOSPHOINOSITIDE_CLASS <- c(
  "MTM1", "MTMR1", "MTMR2", "MTMR3", "MTMR4", "MTMR6", "MTMR7", "MTMR8", "MTMR14",
  "INPP4A", "INPP4B", "INPP5B", "INPP5E", "INPP5F", "INPP5J", "INPP5K", "INPPL1",
  "OCRL", "SYNJ1", "SYNJ2", "SACM1L", "FIG4")

test_that("the phosphoinositide-phosphatase class is lipid-typed via the GO parent", {
  p <- pe_phosphat()
  cls <- p[p$symbol %in% PE_PHOSPHOINOSITIDE_CLASS, , drop = FALSE]
  skip_if(nrow(cls) == 0L, "no phosphoinositide-class members present")
  # Each present member acts on a non-protein lipid substrate (duals included).
  expect_true(all(cls$acts_on_nonprotein))
  expect_true(all(grepl("lipid", cls$nonprotein_substrate_type)))
  # None is silently left at the substrate-less "other" fallback.
  expect_false(any(cls$nonprotein_substrate_type == "other"))
})

test_that("INPP5A boundary canary: a soluble inositol 5-phosphatase is NOT lipid-typed", {
  # INPP5A hydrolyses soluble Ins(1,4,5)P3, not a membrane phosphoinositide, so it sits
  # outside GO:0052866 and must never be pulled into the lipid class by the broad parent.
  # This is the specific gene whose mis-inclusion would signal the parent's boundary moved.
  p <- pe_phosphat()
  row <- pe_row(p, "INPP5A")
  skip_if(is.null(row), "INPP5A absent")
  expect_false(grepl("lipid", row$nonprotein_substrate_type))
})

test_that("the broad lipid parent does not flatten genuine duals to lipid-only", {
  # A class member carrying protein evidence (GO or EC protein) must still be typed dual --
  # adding the lipid axis must never suppress an existing protein call.
  p <- pe_phosphat()
  cls <- p[p$symbol %in% PE_PHOSPHOINOSITIDE_CLASS, , drop = FALSE]
  skip_if(nrow(cls) == 0L, "no phosphoinositide-class members present")
  with_protein_evidence <- cls[cls$go_protein | cls$ec_protein, , drop = FALSE]
  if (nrow(with_protein_evidence)) {
    expect_true(all(with_protein_evidence$acts_on_protein))
    expect_true(all(with_protein_evidence$dual_protein_nonprotein))
  }
})

test_that("inactive myotubularins are not lipid-typed and carry their regulatory annotation", {
  p <- pe_phosphat()
  mtm_pseudo <- p[p$phosphatase_family %in% "Myotubularin" & p$is_pseudophosphatase, , drop = FALSE]
  skip_if(nrow(mtm_pseudo) == 0L, "no myotubularin pseudophosphatases present")
  # Catalytically dead: none is a lipid phosphatase (no non-protein substrate assigned by
  # lineage). A homology-propagated electronic protein-GO term may still flag the protein
  # axis (SBF1) -- that is quarantined into Provisional elsewhere, not asserted here.
  expect_false(any(mtm_pseudo$acts_on_nonprotein))
  # The documented adapters carry a regulatory target rather than a bare untyped row.
  for (sym in c("SBF1", "MTMR12", "MTMR9")) {
    row <- pe_row(p, sym)
    if (!is.null(row)) {
      expect_false(is.na(row$regulates))
      expect_true(nzchar(row$regulatory_role))
    }
  }
})
