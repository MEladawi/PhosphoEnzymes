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
* The obsolete term `GO:0035004` (retired as an unnecessary grouping class for
  "phosphatidylinositol 3-kinase activity") was removed from the kinase GO term
  set and added to the obsolete-GO denylist. PI3K lipid coverage is unchanged:
  it is retained by `GO:0001727` (lipid kinase activity, subtree) and
  `GO:0016303` (1-PI 3-kinase activity, exact), so no gene's membership or typing
  moved.
* The myotubularin family invariant is refined to catalytically active members;
  the family pseudophosphatases (`MTMR9`/`MTMR10`/`MTMR11`/`MTMR12`,
  `SBF1`/`SBF2`) are not lipid-typed.
* New `regulates` / `regulatory_role` columns on the phosphatase master record
  the adapter/activator relationships of the inactive myotubularins in place,
  rather than relocating them out of the catalytic phosphatome.
* `get_kinases()` and `get_phosphatases()` gain a two-knob filter API: `mode`
  (`"comprehensive"` / `"strict"`) selects on rigor (`curated_core`) and
  `substrate` (`"any"` / `"protein"` / `"nonprotein"`) selects on what the
  enzyme acts on. The knobs are orthogonal, so a strictly-curated lipid kinase
  survives `mode = "strict"` unless `substrate = "protein"` is also requested.
* The EC/GO rules that drive substrate typing are now externalized as data:
  four declarative term-set CSVs (kinase/phosphatase x EC/GO) ship in the
  package, readable with `get_term_set()` and lintable with
  `validate_term_set()`.
* `get_kinases()` / `get_phosphatases()` accept a `term_sets =` override that
  re-types the shipped catalog under a user-supplied EC/GO term set without a
  rebuild. The recompute is gene-set-file-free -- it reads each gene's recorded
  evidence from a shipped sidecar -- and the default term set reproduces the
  shipped table exactly (round-trip identity). This path is the only one that
  uses `dplyr` / `purrr` / `stringr` / `readr`, which remain `Suggests`.
* New substrate-provenance columns (`substrate_evidence`, `substrate_decider`,
  `substrate_concordance`) make every substrate call auditable, and a per-gene
  `substrate_evidence.csv` sidecar ships the term-set-independent evidence that
  powers the `term_sets =` recompute.
* Two enrichment-background columns (`is_catalytic_background`,
  `is_protein_catalytic_background`) define the catalytic and protein-catalytic
  universes for downstream over-representation testing.
* The build manifest now stamps the md5 fingerprint of each shipped term-set
  CSV, and the consolidated QC step re-verifies the on-disk term sets against
  those recorded fingerprints, so a term set cannot drift from its provenance
  unnoticed.

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
