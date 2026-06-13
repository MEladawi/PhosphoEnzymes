# PhosphoEnzymes 0.99.0

* Initial skeleton.
* Evidence model: two evidence dimensions by evidence class (structural/evolutionary
  catalog; protein-specific EC), `n_evidence_dimensions` in 0-2 for both
  classes.
* `evidence_tier` (Gold/Silver/Bronze/Provisional): Gold requires both axes;
  Silver = one axis + supplementary support (experimental GO OR reviewed UniProt
  keyword); Bronze = one axis, no supplementary; Provisional = zero axes.
* Two class-specific master tables + a derived unified summary; catalytic-only
  cores with a regulatory-subunit companion.
* Substrate typing carried as co-equal columns: `acts_on_protein`,
  `acts_on_nonprotein`, and the pipe-delimited `nonprotein_substrate_type`
  (empty = protein-only), so dual enzymes (e.g. PIK3CA) are never dropped by a
  single-substrate filter. `catalytic_status` (active/pseudo/uncertain) and the
  derived `is_catalytic_background` distinguish pseudoenzymes; `membership_basis`
  records the deriving source of each structural-catalog call.
