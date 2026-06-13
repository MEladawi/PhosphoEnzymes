# PhosphoEnzymes 0.99.1

* Phosphoinositide phosphatases are now lipid-typed through the
  `GO:0052866` (phosphatidylinositol phosphate phosphatase activity) parent,
  whose ancestor-propagated subtree covers the whole class (myotubularins,
  `INPP4A`/`INPP4B`, the lipid `INPP5` members, `OCRL`, `SYNJ1`/`SYNJ2`,
  `SACM1L`, `FIG4`). Previously several of these carried the subtype-less
  `other` label. `INPP5A` (a soluble Ins(1,4,5)P3 5-phosphatase) is correctly
  excluded and stays non-lipid.
* The obsolete term `GO:0004437` (retired in the ontology, zero annotations)
  was removed from the phosphatase GO term set. `validate_term_set()` now
  hard-errors on a pinned denylist of confirmed-obsolete GO ids, and the build
  re-verifies every GO id's obsolescence against the live ontology on a source
  refresh so the denylist cannot silently go stale (soft-depends on jsonlite).
* The myotubularin family invariant is refined to catalytically active members;
  the family pseudophosphatases (`MTMR9`/`MTMR10`/`MTMR11`/`MTMR12`,
  `SBF1`/`SBF2`) are not lipid-typed.
* New `regulates` / `regulatory_role` columns on the phosphatase master record
  the adapter/activator relationships of the inactive myotubularins in place,
  rather than relocating them out of the catalytic phosphatome.

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
