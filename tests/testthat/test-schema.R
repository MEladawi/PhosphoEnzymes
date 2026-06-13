# Generic schema / property tests.
# These verify the MACHINERY (columns, types, domains) - not biology.
# They will go green on a structurally valid table even if the gate is wrong;
# the biological correctness is guarded by test-trap-cases.R and test-invariants.R.

ALLOWED_SUBSTRATE <- c("protein", "lipid", "nucleotide", "carbohydrate", "other")
ALLOWED_TIER      <- c("Gold", "Silver", "Bronze", "Provisional")

REQUIRED_COLS <- c(
  "ensembl_gene_id", "symbol", "acts_on_protein", "substrate_type",
  "n_evidence_dimensions", "evidence_tier", "curated_core"
)

check_table_schema <- function(df) {
  expect_true(all(REQUIRED_COLS %in% names(df)))

  # keys
  expect_true(all(nzchar(df$ensembl_gene_id)))
  expect_false(anyNA(df$ensembl_gene_id))
  expect_false(any(duplicated(df$ensembl_gene_id)))
  expect_true(all(grepl("^ENSG[0-9]{11}$", df$ensembl_gene_id)))

  # column types / domains
  expect_type(df$acts_on_protein, "logical")
  expect_false(anyNA(df$acts_on_protein))
  expect_type(df$curated_core, "logical")

  expect_true(all(df$substrate_type %in% ALLOWED_SUBSTRATE))
  expect_true(all(df$evidence_tier %in% ALLOWED_TIER))
  expect_true(all(df$n_evidence_dimensions %in% 0:2))

  # cross-field consistency
  # acts_on_protein implies substrate_type == "protein", and vice versa
  expect_equal(df$acts_on_protein, df$substrate_type == "protein")
  # Provisional <=> zero axes <=> not curated_core
  expect_equal(df$evidence_tier == "Provisional",
               df$n_evidence_dimensions == 0L)
  expect_equal(df$curated_core, df$n_evidence_dimensions >= 1L)
  # Gold requires both axes
  expect_true(all(df$n_evidence_dimensions[df$evidence_tier == "Gold"] == 2L))
}

test_that("kinase table satisfies the schema contract", {
  check_table_schema(pe_kinases())
})

test_that("phosphatase table satisfies the schema contract", {
  check_table_schema(pe_phosphat())
})

test_that("unified summary has the agreed thin schema", {
  pe <- pe_unified()
  expect_true(all(c("ensembl_gene_id", "symbol", "regulator_class",
                    "acts_on_protein", "substrate_type",
                    "n_evidence_dimensions", "evidence_tier",
                    "curated_core") %in% names(pe)))
  expect_true(all(pe$regulator_class %in% c("kinase", "phosphatase")))
  # deliberately NO structural family column in the thin summary
  expect_false(any(c("kinase_family", "phosphatase_fold", "family") %in% names(pe)))
})

test_that("strict mode is a subset that drops exactly the Provisional rows", {
  full   <- pe_kinases("comprehensive")
  strict <- pe_kinases("strict")
  expect_true(nrow(strict) <= nrow(full))
  expect_true(all(strict$ensembl_gene_id %in% full$ensembl_gene_id))
  expect_false(any(strict$evidence_tier == "Provisional"))
  expect_setequal(strict$ensembl_gene_id,
                  full$ensembl_gene_id[full$curated_core])
})
