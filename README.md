# PhosphoEnzymes

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20732423-1682D4?logo=zenodo)](https://doi.org/10.5281/zenodo.20732423)

**A reproducible, substrate-aware reference of the human phosphorylation machinery — kinases and phosphatases — keyed on base Ensembl gene IDs and typed by what each enzyme acts on.**

Every gene is placed on two independent axes: **rigor** (how well its catalytic identity is established) and **substrate** (the chemical class it acts on). The two never couple — a strictly-curated lipid kinase is still a Gold-tier enzyme; it is removed from your results only if you explicitly ask for protein substrates.

---

## Quick start

```r
# install from GitHub (not yet on Bioconductor)
# remotes::install_github("MEladawi/PhosphoEnzymes")

library(PhosphoEnzymes)

# Three tables, one row per gene, typed by substrate:
get_kinases()         # 756 human kinases (556 act on protein, incl. 15 dual)
get_phosphatases()    # 298 human phosphatases (181 act on protein, incl. 23 dual)
get_phosphoenzymes()  # 1054-row class-agnostic summary spanning both
```

Two orthogonal filters that never couple — `mode` selects on rigor, `substrate` on what the enzyme acts on:

```r
get_kinases(mode = "strict")                          # curated core (>= 1 evidence dimension)
get_kinases(substrate = "protein")                    # protein-acting only, any rigor
get_phosphatases(mode = "strict", substrate = "nonprotein")
```

`mode` is `"comprehensive"` (default, all rows) or `"strict"`; `substrate` is `"any"` (default), `"protein"`, or `"nonprotein"`. Because the axes are independent, `mode = "strict"` keeps a well-curated lipid kinase — only `substrate = "protein"` drops it.

The EC/GO typing rules ship as cited, overridable data:

```r
get_term_set("kinase", "go")                          # read a rule table
validate_term_set()                                   # lint the shipped defaults
get_kinases(term_sets = my_rules)                     # re-type without rebuilding
provenance(get_kinases())                             # pinned source releases + term-set fingerprints
```

---

## The two-axis model

The single most consequential decision is to type genes by **what they phosphorylate or dephosphorylate**, not by which sequence superfamily they belong to. Sequence-family catalogs place the phosphoinositide kinases (PI3K/PI4K) inside the protein-kinase tree because they are *structurally* related to protein kinases — even though they act on lipids. The same trap recurs on the phosphatase side, where EC class 3.1.3 holds protein, lipid, and metabolite phosphatases side by side. Filtering for "protein enzymes" by family membership would wrongly keep PI4KA, SPHK1, DGKA, the myotubularins, and PTEN's lipid leg.

The organizing principle, then: **sequence-family membership is evidence of lineage, never an override of substrate** (the *PI3K / PTEN principle*). A gene in a classical protein-kinase fold is not forced to be typed as a protein kinase if its biochemistry says it acts on a lipid; PTEN, which dephosphorylates both a protein and the lipid PIP3, is kept intact as a single `dual` row rather than collapsed to one label.

| Axis | Question | Carried as |
|---|---|---|
| **Rigor** | How well-established is this gene as a real enzyme of its class? | `n_evidence_dimensions` (0–2) → `evidence_tier` (Gold / Silver / Bronze / Provisional) |
| **Substrate** | What does it act on? | co-equal booleans `acts_on_protein`, `acts_on_nonprotein` → `substrate_call` (protein / nonprotein / dual / untyped) |

---

## What the tables look like

One row per gene, keyed on a base Ensembl ID (columns trimmed for display — the full tables carry 55 / 53 / 15 columns).

**`get_kinases()`** — PIK3CA sits in a protein-kinase fold yet is typed `dual` with a `lipid` substrate, and the non-protein kinases carry no Manning `kinase_group`: lineage scores rigor but never overrides substrate.

| ensembl_gene_id | symbol | kinase_group | substrate_call | nonprotein_substrate_type | evidence_tier |
|---|---|---|---|---|---|
| ENSG00000146648 | EGFR | TK | protein | — | Silver |
| ENSG00000154229 | PRKCA | AGC | protein | — | Gold |
| ENSG00000121879 | PIK3CA | — | dual | lipid | Gold |
| ENSG00000176170 | SPHK1 | — | nonprotein | lipid | Silver |
| ENSG00000156515 | HK1 | — | nonprotein | carbohydrate | Silver |
| ENSG00000173334 | TRIB1 | Other | untyped | — | Bronze |

`TRIB1` is the worked case for `untyped`: a catalog-confirmed **pseudokinase** (`catalytic_status = pseudo`) with no kinase-activity GO term and no EC number. It clears the rigor gate on lineage but has no substrate signal in any source, so it is reported `untyped` rather than forced to "protein" by lineage.

**`get_phosphatases()`** — PTEN is the canonical bifunctional enzyme, kept as a single `dual` row (protein + lipid) rather than split across two.

| ensembl_gene_id | symbol | phosphatase_fold | substrate_call | nonprotein_substrate_type | evidence_tier |
|---|---|---|---|---|---|
| ENSG00000171862 | PTEN | CC1 | dual | lipid | Gold |
| ENSG00000196396 | PTPN1 | CC1 | protein | — | Gold |
| ENSG00000172531 | PPP1CA | PPPL | protein | — | Gold |
| ENSG00000168918 | INPP5D | DNase I | nonprotein | lipid | Bronze |
| ENSG00000143727 | ACP1 | CC2 | dual | other | Gold |

**`get_phosphoenzymes()`** — the class-agnostic summary spans both classes, carrying only the shared substrate/evidence columns plus `regulator_class`.

| ensembl_gene_id | symbol | regulator_class | substrate_call | evidence_tier | curated_core |
|---|---|---|---|---|---|
| ENSG00000146648 | EGFR | kinase | protein | Silver | TRUE |
| ENSG00000121879 | PIK3CA | kinase | dual | Gold | TRUE |
| ENSG00000171862 | PTEN | phosphatase | dual | Gold | TRUE |
| ENSG00000172531 | PPP1CA | phosphatase | protein | Gold | TRUE |
| ENSG00000156515 | HK1 | kinase | nonprotein | Silver | TRUE |

Class-specific taxonomy (kinase group/family, phosphatase fold/family) stays on the masters; a class-routed join against the summary recovers it.

---

## Evidence model (rigor)

`n_evidence_dimensions` (0–2) counts two **substrate-blind**, structurally distinct evidence classes:

1. **L — structural catalog:** the gene appears in a curated sequence-family catalog (kinase: pkinfam ∪ Manning ∪ KinHub; phosphatase: Chen ∪ HGNC phosphatase groups).
2. **E — biochemistry:** the gene carries at least one EC code matching **any** curated EC rule for its class — protein *or* non-protein alike. The "any class" scope is deliberate: EC membership is EC membership regardless of substrate, so a lipid kinase earns the same biochemical axis as a protein kinase.

These two are distinct in kind, not statistically independent (the catalogs and EC draw on one literature ecosystem). GO terms, UniProt keywords, and the IDG dark-kinome list are **not** counted as dimensions — they share that same provenance, so counting them would measure database fame, not independent confirmation. They contribute as *supplementary support*, which can break a Bronze/Silver tie but can never create a second axis or lift a zero-axis gene out of Provisional.

`evidence_tier` is a practical prioritization label over the two dimensions — **not a probability or a confidence score**:

| Tier | Condition |
|---|---|
| **Gold** | L and E both present (n = 2) |
| **Silver** | one axis **plus** supplementary support (experimental GO or a reviewed UniProt keyword) |
| **Bronze** | one axis, no supplementary support |
| **Provisional** | neither axis (`curated_core = FALSE`) |

Because GO/keyword never reach the second axis, no GO- or keyword-only gene is ever Gold by design. On the shipped build: kinases **Gold 146 / Silver 390 / Bronze 42 / Provisional 178**; phosphatases **Gold 136 / Silver 61 / Bronze 37 / Provisional 64**.

---

## Substrate typing

Substrate is two **co-equal** booleans with no inheritance hierarchy — neither defaults to the other, and lineage never authors a substrate call:

```
acts_on_protein    = go_protein    | ec_protein
acts_on_nonprotein = go_nonprotein | ec_nonprotein | chen_nonprotein   # chen flag: phosphatases only
```

`substrate_call` summarizes the pair as `protein`, `nonprotein`, `dual` (both TRUE), or `untyped` (both FALSE). **`untyped` is a first-class state, not an error** — a curated gene with no substrate evidence in any source is reported as such rather than guessed. `nonprotein_substrate_type` records the chemical class (`lipid`, `nucleotide`, `carbohydrate`, `metabolite`, or `other`), and three fields (`substrate_evidence`, `substrate_decider`, `substrate_concordance`) make every call auditable.

On the shipped build: kinase `substrate_call` **protein 541 / nonprotein 95 / dual 15 / untyped 105**; phosphatase **protein 158 / nonprotein 79 / dual 23 / untyped 38**.

For analyses that need a protein-enzyme background, two ready columns ship: `is_protein_catalytic_background` (active, curated, protein-acting — the right denominator for kinome/phosphatome enrichment) and the substrate-blind `is_catalytic_background` (all active, curated members).

---

## Worked examples

Five self-contained scripts live in `inst/examples/`, each printing its own output. From an installed package, locate and run one with:

```r
library(PhosphoEnzymes)
source(system.file("examples", "01_access_and_filter.R", package = "PhosphoEnzymes"))
```

| Script | Shows |
|---|---|
| `01_access_and_filter.R` | the three accessors and the two orthogonal `mode` / `substrate` knobs |
| `02_substrate_and_dual_pten.R` | the co-equal substrate columns, worked end to end on the dual enzyme PTEN with its provenance |
| `03_rigor_vs_substrate.R` | rigor (`evidence_tier`) and substrate as independent axes, via a tier-by-substrate cross-tab |
| `04_custom_term_sets.R` | reading, linting, and overriding the cited EC/GO term sets (`term_sets=`), including the round-trip identity |
| `05_unified_class_routed_join.R` | recovering class-specific taxonomy from the unified summary with a class-routed join |

The captured console output of each script is committed under `inst/examples/outputs/`, so the expected results can be read without running anything. The default accessor path uses only base R; only `04_custom_term_sets.R` needs the `Suggests` packages (`dplyr`, `purrr`, `stringr`, `readr`) and skips itself with a message if they are absent.

---

## Provenance & reproducibility

The tables ship as package data, regenerated by maintainers from **pinned, dated source snapshots** via `inst/scripts/make-data.R` — with no network access at install or load time. Every input's version and md5 is recorded in `inst/build_manifest.yaml` and stamped as an attribute on each shipped table; `provenance()` returns those release strings and the per-table term-set fingerprints.

Each release is a pinned, citable snapshot, not a live query: the same release always rebuilds to content-identical tables. Snapshots are refreshed through a reviewed pull request after passing the full QC sanity gate, so every release stays a reviewed, pinned snapshot — never an unattended auto-push.

---

## Notes & limitations

- **`untyped` is intentional.** A catalog-confirmed gene with no substrate evidence (TRIB1, the GO-silent pseudoenzymes) is reported `untyped`, not forced to a substrate by lineage. Downstream analyses must decide how to treat these explicitly.
- **Some substrate calls rest on GO** (e.g. NME1's and POMK's protein call). These are tracked as GO-dependent canaries so a flip is diagnosed as a GO-release change, not a code regression.
- **Frozen-vintage catalogs.** KinHub (2017), kinase.com/Manning (2002), and Chen (2017) carry fixed historical classifications; UniProt is the current taxonomy primary on the kinase side.
- **Base Ensembl IDs are stable but not immutable** — Ensembl occasionally retires/merges IDs across major releases, so joins are most robust when both sides use a comparable release.
- **The enzyme universe is defined by the chosen sources and term-set rules**, both isolated in declarative files so membership and typing decisions are auditable and reversible.

The design rationale — the sources, the identifier-resolution approach, the EC-subclass selection, the gate, and the QC gates — is narrated in the package vignette (`vignette("PhosphoEnzymes")`); per-source provenance lives in `inst/extdata/SOURCES.tsv` and the EC/GO typing rules in the cited term-set CSVs under `inst/extdata/`.

---

## License & attribution

Package **code** (`R/`, `inst/scripts/`) is released under the **MIT** license; the bundled **data** tables (`data/*.rda`) under **CC BY 4.0**. The split is possible because every bundled input is itself redistributable:

- **HGNC** — CC0 (identifier bridge, metadata, EC numbers).
- **UniProt** — pkinfam and keywords KW-0418 (Kinase) / KW-0904 (Protein phosphatase); CC BY 4.0.
- **Gene Ontology** — Ensembl-keyed gene sets; CC BY 4.0.
- **Manning 2002 / kinase.com** and the **IDG** dark-kinome list — MIT, via IDG DarkKinaseTools.
- **IUBMB EC numbers** — read as facts via HGNC `enzyme_id`.
- The curated human **pseudokinase** list (the kinase `catalytic_status` source) — compiled from the published pseudokinome literature with a per-row citation.

Two sources carry no redistribution license, so neither is reproduced verbatim — **KinHub** and the **Chen 2017** phosphatome are each bundled as reconstructed, HGNC-normalized gene-membership facts (`inst/extdata/kinhub_facts.tsv`, `inst/extdata/chen_phosphatome_facts.tsv`).

Per-source license, version, URL, and attribution are recorded in `inst/extdata/SOURCES.tsv` and mirrored in `inst/CITATION`. Cite those upstream sources when you publish results derived from the shipped tables.

---

## Citing

The software is archived on Zenodo under the release-independent DOI
[10.5281/zenodo.20732423](https://doi.org/10.5281/zenodo.20732423), which always
resolves to the latest version. See `CITATION.cff` for how to cite the software, and
`inst/extdata/SOURCES.tsv` for the upstream data sources to cite alongside it.
