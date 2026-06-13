# Data sources

All inputs for `build_kinases.R`. Filenames are version-less; the exact version/date of
each file used in a run is recorded automatically in `data_out/source_versions.tsv` and in
the workbook README (along with a build version tag and an MD5 hash of the master table).
Most sources are auto-fetchable: `build_kinases.R` refreshes them (at most once per day) via
`helper_code/source_registry.R`. Call `build_kinase_list(refresh_data = FALSE)` for an
offline, reproducible rerun of the cached files.

| Local file | Source | Auto-fetch URL |
|---|---|---|
| `hgnc_complete_set.txt` | HGNC complete set (identifier bridge) | monthly archive via the `hgnc` package — latest by default (`hgnc::latest_monthly_url()`), or an exact release via `build_kinase_list(hgnc_archive_url=)`; the resolved archive URL is recorded in `hgnc_source_url.txt` and the manifest |
| `pkinfam.txt` | UniProt pkinfam (curated protein kinome) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/pkinfam.txt |
| `go_mf_genesets_with_iea_ensembl.gmt` (default) / `go_mf_genesets_no_iea_ensembl.gmt` | GO molecular-function gene sets, Ensembl-keyed (Bader Lab EM_Genesets). IEA-inclusive by default; `build_kinase_list(go_include_iea = FALSE)` selects the no-IEA variant. | https://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/GO/Human_GO_mf_with_GO_iea_ensembl.gmt · https://download.baderlab.org/EM_Genesets/current_release/Human/ensembl/GO/Human_GO_mf_no_GO_iea_ensembl.gmt |
| `kinase.com_manning_list.xls` | Manning human kinome (kinase.com, 2002) | https://raw.githubusercontent.com/IDG-Kinase/DarkKinaseTools/master/data-raw/dark_kinases/kinase.com_list.xls |
| `kinhub_facts.tsv` | KinHub human kinase list (kinhub.org; Eid et al. 2017) — **reconstructed facts**: HGNC-normalized gene memberships and the Manning group/family/subfamily labels, one row per gene. Pinned (not auto-fetched). The source web page itself is not redistributed. | http://www.kinhub.org/kinases.html (basis of the reconstruction) |
| `uniprot_kinase_KW-0418_human.tsv` | UniProtKB reviewed human, keyword KW-0418 "Kinase" (membership + `protein_families` taxonomy) | https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)AND(reviewed:true)AND(keyword:KW-0418)&fields=accession,gene_primary,protein_name,ec,protein_families |
| `IDG_dark_kinase_list.csv` | IDG understudied ("dark") kinome | https://github.com/IDG-Kinase/DarkKinaseTools (data-raw/dark_kinases/Dark Kinase List.csv) |

Notes
- The GO MF GMT must be **ancestor-propagated** (the protein-kinase gate relies on it); the
  loader asserts this at startup (member count + child-only canary kinases) and fails the build
  with a clear message otherwise. IEA evidence is included by default as a deliberate
  recall-over-precision choice. The variant used is recorded in the manifest.
- HGNC is the sole identifier authority: pkinfam / Manning / KinHub / UniProt-KW / IDG are
  mapped to a base Ensembl gene ID through HGNC (by Entrez, then UniProt accession, then
  current/alias/previous symbol). The GO gene sets are already Ensembl-keyed, so they need
  no mapping.
- Kinase taxonomy is assigned per field (each "first non-blank" across the listed sources):
  - `kinase_group`     : UniProt -> KinHub -> kinase.com. The UniProt value is the leading
                         token of its `protein_families` string (with `Tyr` -> `TK`), accepted
                         ONLY if it is a real Manning group (AGC/CAMK/CK1/CMGC/NEK/RGC/STE/TK/
                         TKL/Other/Atypical); otherwise NA so KinHub/kinase.com fill it. Non-
                         protein kinases get no group.
  - `kinase_family`    : KinHub -> kinase.com only (the Manning short tier: Akt, CDK, FGFR ...).
                         UniProt's verbose family string is NOT used here.
  - `kinase_subfamily` : UniProt -> KinHub -> kinase.com.
  - `uniprot_protein_family` : the raw UniProt `Protein families` string, verbatim.
  KinHub (2017) and kinase.com (2002) are static, so UniProt keeps current kinases assigned
  to the (fixed) Manning group/subfamily scheme.
- Resolution guards against stale source identifiers: a candidate that is an RNA-gene locus,
  or whose HGNC symbol history disagrees with the source symbol, is rejected in favour of a
  consistent hit (see `helper_code/hgnc_bridge.R`).
- Downstream: match your data on base Ensembl ID (preferred) or HGNC symbol. If matching
  on symbols, updating outdated symbols first (e.g. `HGNChelper::checkGeneSymbols()`) avoids
  misses; matching on Ensembl IDs avoids the symbol-versioning problem entirely.
