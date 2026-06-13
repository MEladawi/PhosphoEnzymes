test_that("apply_term_sets computes substrate-blind rigor, axis flags, and provenance", {
  source(build_file("utils.R"),     local = TRUE)
  source(build_file("term_sets.R"), local = TRUE)
  ts  <- load_term_sets(extdata_dir())
  tmp <- tempfile(fileext = ".gmt"); writeLines("X%GOMF%GO:0004672\tna\tENSGP", tmp)
  res <- resolve_term_sets(ts, tmp)
  ev <- tibble::tibble(
    ensembl_gene_id = c("g_pi4ka","g_lhpp_like","g_untyped","g_dual"),
    all_ec_codes    = list("2.7.1.67", "3.6.1.1", character(0), c("3.1.3.16","3.1.3.67")),
    go_protein            = c(FALSE, FALSE, FALSE, TRUE),
    go_nonprotein         = c(TRUE,  FALSE, FALSE, FALSE),
    go_nonprotein_subtype = c("lipid","","",""),
    chen_nonprotein       = c(FALSE, TRUE,  FALSE, FALSE),
    in_structural_catalog = c(TRUE,  TRUE,  TRUE,  TRUE),
    supplementary_support = c(TRUE,  TRUE,  FALSE, TRUE),
    go_experimental_protein = c(FALSE, FALSE, FALSE, TRUE),
    catalytic_status      = "active")
  outk <- apply_term_sets(ev[1, ], res$kinase)
  expect_true(outk$ec_rigor); expect_equal(outk$n_evidence_dimensions, 2L)
  expect_false(outk$acts_on_protein); expect_true(outk$acts_on_nonprotein)
  expect_identical(outk$nonprotein_substrate_type, "lipid")
  outp <- apply_term_sets(ev[2:4, ], res$phosphatase)
  expect_equal(outp$evidence_tier[1], "Gold")
  expect_identical(outp$substrate_call[2], "untyped")
  expect_true(outp$acts_on_protein[3] && outp$acts_on_nonprotein[3])
  expect_identical(outp$substrate_call[3], "dual")
})
