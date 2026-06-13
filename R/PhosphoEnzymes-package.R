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
