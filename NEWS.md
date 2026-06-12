# PhosphoEnzymes 0.99.0

* Initial skeleton.
* Evidence model: two independent axes by evidence type (structural/evolutionary
  catalog; protein-specific EC), `n_independent_evidence_axes` in 0-2 for both
  classes.
* `evidence_tier` (Gold/Silver/Bronze/Provisional): Gold requires both axes;
  Silver = one axis + supplementary support (experimental GO OR reviewed UniProt
  keyword); Bronze = one axis, no supplementary; Provisional = zero axes.
* Two class-specific master tables + a derived unified summary; catalytic-only
  cores with a regulatory-subunit companion.
