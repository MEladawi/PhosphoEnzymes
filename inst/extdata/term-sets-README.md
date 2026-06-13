# EC/GO Term Sets for Kinase and Phosphatase Classification

This directory contains four cited reference tables that define which EC numbers and GO
Molecular Function terms are used to classify human kinases and phosphatases. They are
the default term sets shipped with the package and are user-overridable via the package's
term-set API.

## Files

| File | Rows | Enzyme class |
|---|---|---|
| `kinase_ec_terms.csv` | 23 | Kinase — EC 2.7 subclasses and related |
| `kinase_go_terms.csv` | 22 | Kinase — GO:MF kinase activity subtree |
| `phosphatase_ec_terms.csv` | 14 | Phosphatase — EC 3 hydrolase entries |
| `phosphatase_go_terms.csv` | 17 | Phosphatase — GO:MF phosphatase activity subtree |

## Column schema

All four files share the same 8-column schema:

| Column | Type | Description |
|---|---|---|
| `term_id` | character | Canonical identifier — an EC number (e.g. `2.7.10`) or a GO accession (e.g. `GO:0004672`). |
| `class` | character | Enzyme class this row belongs to: `kinase` or `phosphatase`. |
| `substrate` | character | Broad substrate category the term implies. One of `protein`, `nonprotein`, or `na` (substrate is unknown or the term is class-membership-only). |
| `substrate_subtype` | character | Finer substrate category within `nonprotein`. One of `lipid`, `nucleotide`, `carbohydrate`, `metabolite`, `other`, or empty for rows where `substrate` is `protein` or `na`. |
| `role` | character | How the term is used in classification logic. See **Role values** below. |
| `scope` | character | How the term is matched against a gene's annotations. One of `subclass` (any EC term whose number starts with this prefix), `exact` (only this exact EC or GO accession), or `subtree` (this GO term and all of its descendants via `is_a`/`part_of`). |
| `citation` | character | Primary source(s) that define or motivate the term's inclusion or exclusion. Citations refer to IUBMB enzyme nomenclature, the GO Consortium, AmiGO, or primary literature as recorded in this column. |
| `note` | character | Free-text rationale, especially for non-obvious inclusions, exclusions, or asymmetries. |

### Role values

- **`rigor+substrate`** — the term contributes to the class rigor metric (i.e. it positively
  confirms that the gene is a kinase or phosphatase) *and* its `substrate` tag is used when
  typing the gene's substrate specificity. A lipid-kinase EC term with `role = rigor+substrate`
  and `substrate = nonprotein` still counts toward kinase rigor even though the substrate is
  not protein.

- **`rigor_umbrella`** — the term establishes class membership only. It contributes to the
  rigor count but does not vote on substrate type. Used for GO terms that cover broad
  phosphotransferase or phosphohydrolase activity without specifying substrate.

- **`exclude`** — a named trap that explicitly does not count toward membership, rigor, or
  substrate typing. Exclude rows prevent ambiguous terms from being treated as positive
  evidence. Examples: `EC 3.1.3.3` (PSPH, a phosphoserine phosphatase excluded from the
  phosphatase gate because it acts on a metabolite rather than a phosphoprotein) and
  `GO:0004647` (phosphoserine phosphatase activity, excluded for the same reason);
  `EC 2.7.1.37` (the legacy protein-kinase EC that has been superseded by 2.7.10–2.7.14 and
  would double-count).

## Structural asymmetries

### Kinase EC set: subclass matching for the protein-kinase block

EC subclasses 2.7.10–2.7.14 are matched by `scope = subclass`: any gene annotated to any
EC number whose prefix falls in this range is treated as a protein kinase. This is
appropriate because enzyme nomenclature defines 2.7.10 (receptor Tyr kinases), 2.7.11
(non-receptor Ser/Thr kinases), 2.7.12 (dual-specificity kinases), 2.7.13 (His kinases),
and 2.7.14 (Arg kinases) as protein-kinase subclasses by definition — every member of
these subclasses is a protein kinase.

### Phosphatase EC set: exact matching with out-of-class entries

The phosphatase EC set uses `scope = exact` for most entries because EC 3.1.3
(phosphoric-monoester hydrolases) interleaves protein-directed phosphatases (e.g. EC 3.1.3.16,
protein-serine phosphatase) with non-protein ones (e.g. glucose-6-phosphatase). Subclass
matching would pull in the entire heterogeneous EC 3.1.3 group.

Two entries deliberately fall outside EC 3.1 yet are included:

- **EC 3.9.1.3** (`scope = exact`) — the enzyme encoded by *PHPT1* (phosphohistidine
  phosphatase 1), a protein phosphatase that cleaves phosphohistidine. It is classified
  under EC 3.9 (acting on P–N bonds) rather than 3.1 and would be missed by any EC 3.1-only
  rule.

- **EC 3.6.1.1** (`scope = exact`, `role = rigor+substrate`, `substrate = nonprotein`) —
  *LHPP* (phospholysine phosphohistidine inorganic pyrophosphate phosphatase), classified
  under EC 3.6 (acting on acid anhydrides). Its annotation confirms a phosphatase, but the
  substrate is not a phosphoprotein, so it contributes to rigor without voting for protein
  substrate.

## Rigor versus substrate typing

**Rigor** (does this gene qualify as a kinase or phosphatase?) and **substrate typing**
(is its primary substrate protein or non-protein?) are computed independently:

- Any row with `role = rigor+substrate` or `role = rigor_umbrella` contributes to class
  rigor regardless of what its `substrate` column says.
- Substrate typing uses a precedence hierarchy: GO annotation wins over curated catalog
  flags, which win over EC. The `substrate` and `substrate_subtype` values in these tables
  are used at the EC tier of that hierarchy.

## Overriding the default term sets

These four CSVs are the defaults used when the package's classification functions are
called without explicit term-set arguments. The package's term-set API (shipped separately)
accepts any data frame matching this 8-column schema, allowing users to add custom terms,
narrow the scope of an existing entry, or add new exclusions without forking the package.
