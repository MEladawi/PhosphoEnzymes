#' Access the human kinase reference table
#'
#' @param mode `"comprehensive"` (the full membership union, default) or
#'   `"strict"` (the canonical set: genes with at least one independent
#'   evidence axis, i.e. an expert structural catalog and/or a protein-specific
#'   EC; equivalently `curated_core == TRUE` / not Provisional).
#' @return A [tibble][tibble::tibble] of human kinases, one row per gene.
#' @examples
#' k <- get_kinases()
#' nrow(k)
#' get_kinases(mode = "strict")
#' @export
get_kinases <- function(mode = c("comprehensive", "strict")) {
  mode <- match.arg(mode)
  .pe_get_class("human_kinases", mode)
}

#' Access the human phosphatase reference table
#'
#' @inheritParams get_kinases
#' @return A [tibble][tibble::tibble] of human phosphatases, one row per
#'   catalytic gene. Regulatory/scaffold subunits are not included here; see
#'   [get_phosphatase_regulators()].
#' @examples
#' \dontrun{
#' p <- get_phosphatases()
#' table(p$substrate_type)
#' }
#' @export
get_phosphatases <- function(mode = c("comprehensive", "strict")) {
  mode <- match.arg(mode)
  .pe_get_class("human_phosphatases", mode)
}

#' Access the unified phospho-enzyme summary
#'
#' The thin, derived cross-class summary (kinases + phosphatases) with the
#' shared, class-agnostic columns only. "All protein-acting phospho-enzymes"
#' is `subset(get_phosphoenzymes(), acts_on_protein)`.
#'
#' @return A [tibble][tibble::tibble] combining both classes.
#' @examples
#' \dontrun{
#' pe <- get_phosphoenzymes()
#' with(pe, table(regulator_class, acts_on_protein))
#' }
#' @export
get_phosphoenzymes <- function() {
  .pe_load("human_phosphoenzymes")
}

#' Access the phosphatase regulatory/scaffold subunit companion table
#'
#' Non-catalytic subunits (e.g. `PPP1R*`, `PPP2R*`) mapped to their catalytic
#' complexes. These are deliberately excluded from [get_phosphatases()] to keep
#' the phosphatase reference catalytic-only (the analog of excluding cyclins
#' from a kinome).
#'
#' @return A [tibble][tibble::tibble] of regulatory/scaffold subunits.
#' @export
get_phosphatase_regulators <- function() {
  .pe_load("phosphatase_regulatory_subunits")
}

# ---- internal helpers -------------------------------------------------------

#' @keywords internal
#' @noRd
.pe_get_class <- function(dataset, mode) {
  df <- .pe_load(dataset)
  if (mode == "strict") {
    # Strict = canonical set: at least one independent evidence axis.
    # Equivalent to curated_core == TRUE / evidence_tier != "Provisional".
    df <- df[df$n_independent_evidence_axes >= 1L, , drop = FALSE]
  }
  df
}

#' @keywords internal
#' @noRd
.pe_load <- function(dataset) {
  e <- new.env(parent = emptyenv())
  # An absent dataset makes utils::data() warn and create nothing; suppress that
  # warning so the get() below raises the single, clean "not built yet" error that
  # the test helpers turn into a skip().
  suppressWarnings(utils::data(list = dataset, package = "PhosphoEnzymes", envir = e))
  tibble::as_tibble(get(dataset, envir = e))
}
