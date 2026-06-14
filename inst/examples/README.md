# PhosphoEnzymes — worked examples

Runnable scripts demonstrating the package surface. Each is self-contained and
prints its output; run any of them after `library(PhosphoEnzymes)`.

From an installed package you can locate and source one with:

```r
library(PhosphoEnzymes)
source(system.file("examples", "01_access_and_filter.R", package = "PhosphoEnzymes"))
```

or just open the file and step through it.

| Script | Shows |
|---|---|
| `01_access_and_filter.R` | the three accessors and the two orthogonal `mode` / `substrate` knobs |
| `02_substrate_and_dual_pten.R` | the co-equal substrate columns and the dual enzyme PTEN, end to end with its provenance |
| `03_rigor_vs_substrate.R` | that rigor (`evidence_tier`) and substrate are independent axes |
| `04_custom_term_sets.R` | reading, linting, and overriding the cited EC/GO term sets (`term_sets=`), including the round-trip identity |
| `05_unified_class_routed_join.R` | recovering class-specific taxonomy from the unified summary via a class-routed join |

The default accessor path uses only base R. The `term_sets=` recompute in
`04_custom_term_sets.R` is the one path that needs the `Suggests` packages
(`dplyr`, `purrr`, `stringr`, `readr`); that script skips itself with a message
if they are not installed.

The captured console output of each script is committed alongside it under
`outputs/` (e.g. `outputs/02_substrate_and_dual_pten.txt`), so the expected
results can be read without running anything.
