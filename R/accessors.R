#' Access the human kinase reference table
#'
#' @param mode `"comprehensive"` (the full membership union, default) or
#'   `"strict"` (the canonical set: genes with at least one evidence dimension,
#'   i.e. an expert structural catalog and/or a protein-specific EC; equivalently
#'   `curated_core == TRUE` / not Provisional).
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

#' Access the human protein-phosphatase reference table
#'
#' @inheritParams get_kinases
#' @return A [tibble][tibble::tibble] of human phosphatases, one row per gene.
#' @examples
#' p <- get_phosphatases()
#' nrow(p)
#' get_phosphatases(mode = "strict")
#' @export
get_phosphatases <- function(mode = c("comprehensive", "strict")) {
  mode <- match.arg(mode)
  .pe_get_class("human_phosphatases", mode)
}

#' Access the unified human phospho-enzyme summary
#'
#' A thin, class-agnostic table spanning both kinases and phosphatases (one row per gene,
#' `regulator_class` distinguishing them), carrying only the shared substrate/evidence columns.
#' Join to [get_kinases()] / [get_phosphatases()] by `ensembl_gene_id` to recover class-specific
#' taxonomy. Ships all rows including Provisional; filter on `evidence_tier` / `curated_core`.
#'
#' @return A [tibble][tibble::tibble] spanning both enzyme classes.
#' @examples
#' pe <- get_phosphoenzymes()
#' table(pe$regulator_class)
#' @export
get_phosphoenzymes <- function() {
  .pe_load("human_phosphoenzymes")
}

# ---- internal helpers -------------------------------------------------------

#' @keywords internal
#' @noRd
.pe_get_class <- function(dataset, mode) {
  df <- .pe_load(dataset)
  if (mode == "strict") {
    # Strict = canonical set: at least one evidence dimension.
    # Equivalent to curated_core == TRUE / evidence_tier != "Provisional".
    df <- df[df$n_evidence_dimensions >= 1L, , drop = FALSE]
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
