#' PhosphoEnzymes: a reproducible reference of human kinases and phosphatases
#'
#' A function-first reference resource for the human phosphorylation machinery,
#' keyed on base Ensembl gene IDs and typed by substrate class. It ships the
#' kinase master ([get_kinases()]), the protein-phosphatase master
#' ([get_phosphatases()]), and a thin class-agnostic summary spanning both
#' ([get_phosphoenzymes()]); a regulatory-subunit companion is forthcoming. See
#' the package vignette for the evidence model, the `evidence_tier` definition,
#' and usage.
#'
#' @section Evidence model (brief):
#' `n_evidence_dimensions` (0-2) counts distinct evidence *classes*:
#' (1) a structural/evolutionary sequence-family catalog and (2) a
#' protein-specific 4-digit EC. It is the rigor metric. `evidence_tier`
#' (Gold/Silver/Bronze/Provisional) is a practical prioritization heuristic
#' that additionally uses supplementary GO + reviewed UniProt keyword support;
#' it is not a probability or confidence score.
#'
#' @keywords internal
"_PACKAGE"

# Data-masked column names referenced inside the dplyr verbs of the term-set
# runtime twin (R/term-sets-runtime.R). Declared so R CMD check does not flag
# them as undefined globals; they are columns of the term-set tables and the
# per-gene evidence tibble.
utils::globalVariables(c(
  "role", "substrate", "substrate_subtype", "term_id", "scope",
  "all_ec_codes", "in_structural_catalog", "supplementary_support",
  "go_protein", "go_nonprotein", "go_nonprotein_subtype", "chen_nonprotein",
  "go_experimental_protein", "catalytic_status",
  "ec_rigor", "ec_protein", "ec_nonprotein", "ec_nonprotein_subtype",
  "n_evidence_dimensions", "curated_core", "acts_on_protein",
  "acts_on_nonprotein", "substrate_call",
  # referenced inside the dplyr verbs of .pe_reclassify (R/accessors.R)
  "regulator_class", "ec_codes", "ensembl_gene_id", "go_terms", "gene_go"))
