# Human kinase reference table

A small, reproducible R pipeline that builds a **comprehensive human kinase gene list** by
integrating seven authoritative kinase resources, keyed on **base (unversioned) Ensembl gene
IDs**. Every gene is typed by enzymatic **function** (protein vs lipid / sugar / nucleotide /
etc.), and every source is exposed as an independent binary column so the kinome can be
filtered per source for any project.

Use it to select or annotate kinase genes in any analysis — filtering an expression matrix,
annotating a hit list, comparing assays, etc. Match your data on Ensembl ID (preferred) or
HGNC symbol.

## What it produces (`data_out/`)

| File | Contents |
|---|---|
| `human_kinases_master.csv` | the full table, one row per kinase gene |
| `human_kinases_master.xlsx` | same data (styled) + a README sheet with method, sources, counts |
| `kinases_ensembl_all.txt` | all Ensembl IDs, one per line |
| `kinases_ensembl_protein.txt` | Ensembl IDs of protein kinases only |
| `kinases_ensembl_highconf.txt` | Ensembl IDs with high confidence (≥2 sources or EC evidence) |
| `kinases_symbols_all.txt` / `_protein.txt` | the same selections as HGNC symbols |
| `kinases_unmapped.csv` | any source gene that could not be mapped (never dropped silently) |
| `source_versions.tsv` | the exact version/date of each source used in the run |

## Sources

| Source | Role |
|---|---|
| HGNC complete set | **identifier bridge** (Entrez / UniProt / symbol → Ensembl); also gene metadata |
| UniProt pkinfam | curated protein kinome |
| Manning kinome (kinase.com) | classical kinome + group/family/subfamily taxonomy |
| KinHub (kinhub.org) | Manning classification with current HGNC names (taxonomy fallback) |
| GO molecular-function gene sets (Bader Lab, Ensembl-keyed) | functional umbrella, the protein-kinase gate, and non-protein classes |
| HGNC EC 2.7 subclasses | phosphotransferase (kinase) enzyme classes |
| UniProt keyword KW-0418 ("Kinase") | broad, all-kinase-type membership **and** the primary group/subfamily taxonomy (UniProt `protein_families`) |
| IDG understudied ("dark") kinome | annotation of understudied kinases |

All sources auto-download (see `data_in/SOURCES.md` for URLs). They are cached in `data_in/`
with version-less filenames; the version actually used each run is recorded in
`source_versions.tsv` and the workbook README.

## Method (the functional gate)

Membership in the kinase list is the **union** of the seven source legs (mapped to Ensembl
via HGNC). Each gene is then typed by function: a gene with a *non-protein* GO kinase
activity is typed non-protein **only if** it is not annotated with GO protein-kinase activity
(`GO:0004672`). This single gate is what classifies e.g. PI4KA/SPHK1/DGKA as lipid kinases
but PRKCA/ATM/MTOR/PIK3CA as protein kinases. Non-protein typing priority is
lipid → inositol-phosphate → carbohydrate → nucleotide → creatine, with an EC-subclass
fallback.

Taxonomy is assigned per field, each "first non-blank" across the listed sources:

- `kinase_group` — UniProt → KinHub → kinase.com. UniProt's value is the leading token of its
  `protein_families` string (with `Tyr`→`TK`), accepted only if it is a real Manning group
  (AGC/CAMK/CK1/CMGC/NEK/RGC/STE/TK/TKL/Other/Atypical); otherwise it is left blank for
  KinHub/kinase.com to fill (non-protein kinases get no group).
- `kinase_family` — KinHub → kinase.com only: the Manning short tier (Akt, CDK, FGFR …). The
  verbose UniProt family string is **not** mixed into this column.
- `kinase_subfamily` — UniProt → KinHub → kinase.com.
- `uniprot_protein_family` — the raw UniProt `Protein families` string, kept verbatim.

UniProt is the actively-curated, up-to-date source; KinHub (2017) and kinase.com (2002) are
static, so they serve as fallbacks for the fixed Manning group/family/subfamily scheme.

Identifier resolution is guarded against stale source IDs: a candidate gene that is an RNA
locus, or whose HGNC symbol history (current + alias + previous) disagrees with the source
symbol, is rejected in favour of a consistent hit.

## Requirements

- R (≥ 4.1)
- Packages: `readr dplyr tidyr stringr purrr tibble rvest readxl openxlsx`
  (the script installs any that are missing).
- Internet access on first run (or whenever refreshing sources).

## Usage

From the shell:

```sh
Rscript build_kinases.R
```

This refreshes the auto-updatable sources (at most once per day), builds the table, writes
`data_out/`, and prints a QC report that asserts a set of sanity genes.

Everything runs through a single function, `build_kinase_list()`, so you can also drive it
from an R session:

```r
invisible(lapply(list.files("helper_code", pattern = "[.]R$", full.names = TRUE), source))

result <- build_kinase_list(
  refresh_data = TRUE,   # FALSE = offline, reproducible rerun of cached files in data_in/
  data_in_dir  = "data_in",
  output_dir   = "data_out",
  write_files  = TRUE,   # FALSE = build in memory only (no files written)
  quiet        = FALSE)  # TRUE = no progress/QC messages

result$kinases        # the table (data frame)
result$sanity_passed  # TRUE if all QC sanity genes passed
```

## Repo structure

```
build_kinases.R            thin entry point: loads helper_code/ and calls build_kinase_list()
helper_code/
  pipeline.R               build_kinase_list(): the end-to-end pipeline + messaging
  utils.R                  shared helpers
  source_registry.R        all sources + the updater + version manifest
  hgnc_bridge.R            HGNC reader, lookup maps, resolve_to_ensembl()
  go_functional_sets.R     GO gene sets (Ensembl-keyed), selected by GO accession
  source_pkinfam.R         one file per data source ...
  source_manning.R
  source_kinhub.R          (+ build_kinase_taxonomy)
  source_uniprot_kw.R
  source_idg.R
  source_ec.R
  classify.R               per-gene functional classification + table assembly
  outputs.R                writes files/workbook + QC report
data_in/                   cached source files (+ SOURCES.md)
data_out/                  generated outputs
```

## Extending it

- **Add a source:** add an entry to `SOURCE_REGISTRY` (`helper_code/source_registry.R`) and a
  small `load_*()` that returns `list(ensembl_ids = ..., unmapped = ...)`, then wire it into
  the `membership` list in `build_kinases.R`.
- **Adjust GO typing:** edit the GO accessions in `helper_code/go_functional_sets.R`.
- **Change kinase EC subclasses:** edit the vectors in `helper_code/source_ec.R`.

## Reproducibility & attribution

For a fixed result, keep the files in `data_in/` and run with `REFRESH_SOURCES <- FALSE`;
`source_versions.tsv` records what a given build used. If you publish results derived from
this table, cite the underlying data sources (HGNC, UniProt, Manning et al. 2002, KinHub /
Eid et al. 2017, the Gene Ontology, and the IDG program) per their terms.

## Citation

If you use this software or the kinase table it produces, please cite it. Citation metadata
is in [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button generated
from it). This citation covers the **pipeline**; also cite the upstream data sources listed
above for the underlying annotations.

## License

The source code is released under the [MIT License](LICENSE) — free to use, modify, and
redistribute, provided the copyright and license notice are retained. The MIT license applies
to the **code** only; the generated kinase table is derived from third-party data sources,
each of which carries its own terms that downstream users must respect.
