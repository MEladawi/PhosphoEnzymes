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
  utils::data(list = dataset, package = "PhosphoEnzymes", envir = e)
  tibble::as_tibble(get(dataset, envir = e))
}
