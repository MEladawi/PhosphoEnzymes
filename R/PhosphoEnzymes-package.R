#' PhosphoEnzymes: a reproducible reference of human kinases and phosphatases
#'
#' A function-first reference resource for the human phosphorylation machinery.
#' Two class-specific master tables (kinases, phosphatases) are keyed on base
#' Ensembl gene IDs and typed by substrate class; a derived unified summary
#' ([get_phosphoenzymes()]) provides a thin cross-class view. See the package
#' vignette for the evidence model, the `evidence_tier` definition, and usage.
#'
#' @section Evidence model (brief):
#' `n_independent_evidence_axes` (0-2) counts independent evidence *types*:
#' (1) a structural/evolutionary sequence-family catalog and (2) a
#' protein-specific 4-digit EC. It is the rigor metric. `evidence_tier`
#' (Gold/Silver/Bronze/Provisional) is a practical prioritization heuristic
#' that additionally uses supplementary GO + reviewed UniProt keyword support;
#' it is not a probability or confidence score.
#'
#' @keywords internal
"_PACKAGE"
