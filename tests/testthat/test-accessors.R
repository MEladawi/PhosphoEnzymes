test_that("two knobs filter rigor and substrate independently; strict is not a protein filter", {
  k_all    <- PhosphoEnzymes::get_kinases()
  k_strict <- PhosphoEnzymes::get_kinases(mode = "strict")
  k_prot   <- PhosphoEnzymes::get_kinases(substrate = "protein")
  expect_true(all(k_strict$curated_core))
  expect_true(all(k_prot$acts_on_protein))
  expect_true(nrow(k_strict) <= nrow(k_all))
  # strict must NOT imply protein: there is at least one strict, non-protein gene
  expect_true(any(k_strict$acts_on_nonprotein & !k_strict$acts_on_protein))
  # the dissertation path
  kpp <- PhosphoEnzymes::get_kinases(mode = "strict", substrate = "protein")
  expect_true(all(kpp$curated_core & kpp$acts_on_protein))
})

test_that("get_phosphatases / get_phosphoenzymes honor substrate", {
  p <- PhosphoEnzymes::get_phosphatases(substrate = "nonprotein")
  expect_true(all(p$acts_on_nonprotein))
  pe <- PhosphoEnzymes::get_phosphoenzymes(substrate = "protein")
  expect_true(all(pe$acts_on_protein))
})

test_that("provenance returns the stamped term-set md5s", {
  pv <- PhosphoEnzymes::provenance(PhosphoEnzymes::get_kinases())
  expect_true(all(c("kinase_ec","kinase_go","phosphatase_ec","phosphatase_go") %in% names(pv$term_set_md5)))
})

test_that("get_term_set returns the schema; validate_term_set passes default and catches a fault", {
  ts <- PhosphoEnzymes::get_term_set("kinase", "ec")
  expect_true(all(c("term_id","class","substrate","role","scope") %in% names(ts)))
  issues <- PhosphoEnzymes::validate_term_set()
  expect_false(any(issues$severity == "error"))
  bad <- PhosphoEnzymes::get_term_set("kinase", "ec")
  dup <- bad[bad$term_id == "2.7.10.-", , drop = FALSE]; dup$substrate <- "nonprotein"
  bad2 <- rbind(bad, dup)
  iss <- PhosphoEnzymes::validate_term_set(list(kinase_ec = bad2))
  expect_true(any(iss$severity == "error" & grepl("overlap", iss$message)))
})
