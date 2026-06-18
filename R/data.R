#' Human kinase reference table
#'
#' One row per human kinase gene, keyed on base Ensembl gene ID and typed by
#' substrate class. Built by `inst/scripts/make-data.R` from pinned source
#' snapshots; see the vignette for the evidence model.
#'
#' @format A data frame with one row per gene and the following columns:
#' \describe{
#'   \item{ensembl_gene_id}{Base (unversioned) Ensembl gene ID. Primary key.}
#'   \item{symbol}{HGNC symbol.}
#'   \item{acts_on_protein}{Logical. TRUE if the enzyme acts on protein substrates.}
#'   \item{acts_on_nonprotein}{Logical. TRUE if the enzyme acts on a non-protein
#'     substrate. Co-equal with `acts_on_protein` -- a dual enzyme (e.g. PIK3CA) is
#'     TRUE for both; the two are filtered orthogonally by the `substrate` accessor knob.}
#'   \item{nonprotein_substrate_type}{Pipe-delimited non-protein substrate classes
#'     ("lipid", "nucleotide", "carbohydrate", "metabolite", "other"); empty for
#'     protein-only.}
#'   \item{substrate_subtype}{Coarse substrate-type bucket derived from
#'     `substrate_call` and `nonprotein_substrate_type`: one of "Protein kinase",
#'     "Lipid kinase", "Carbohydrate/sugar kinase", "Nucleotide/nucleoside kinase",
#'     "Metabolite kinase", "Other small-molecule kinase", or "Unclassified kinase".
#'     For the finer non-protein class see `nonprotein_substrate_type`.}
#'   \item{substrate_call}{One of "protein" / "nonprotein" / "dual" / "untyped" --
#'     the collapsed substrate verdict derived from `acts_on_protein` /
#'     `acts_on_nonprotein`.}
#'   \item{substrate_evidence}{Which evidence kinds drove the substrate call, joined
#'     with "+": any of "GO", "EC".}
#'   \item{substrate_concordance}{"concordant" (>1 substrate-evidence kind agrees),
#'     "single" (one kind), or "untyped".}
#'   \item{substrate_decider}{The precedence-winning evidence that set the substrate
#'     call ("GO-exp", "GO-elec", "EC", or "precedence-default").}
#'   \item{ec_protein}{Logical. A protein-substrate EC rule fired for this gene.}
#'   \item{ec_nonprotein}{Logical. A non-protein-substrate EC rule fired.}
#'   \item{go_protein}{Logical. A protein-substrate GO kinase-activity term fired.}
#'   \item{go_nonprotein}{Logical. A non-protein-substrate GO kinase-activity term fired.}
#'   \item{catalytic_status}{"active" / "pseudo" / "uncertain". `pseudo` is sourced from
#'     the curated pseudokinase set (`pseudokinases.csv`, 32 genes; TRIB1 etc.). `active`
#'     is the complement-default (a kinase not on the curated pseudo list), not an
#'     individually-verified catalytic assertion per gene: the active/pseudo distinction
#'     is properly sourced, so read "active" as "not curated-pseudo", not "catalytically
#'     confirmed". A soft signal (a pseudokinase may carry a legacy EC), not a veto.}
#'   \item{is_catalytic_background}{Logical. Substrate-blind enrichment background:
#'     `catalytic_status == "active"` AND `curated_core` (excludes pseudoenzymes and
#'     zero-dimension genes), spanning protein and non-protein substrates alike.}
#'   \item{is_protein_catalytic_background}{Logical. The default protein-kinase
#'     enrichment background: `catalytic_status == "active"` AND `curated_core` AND
#'     `acts_on_protein`. Use this when the universe is protein kinases.}
#'   \item{n_evidence_dimensions}{Integer. Count of distinct class-evidence kinds --
#'     structural/sequence catalog and any class-specific EC -- confirming a bona fide
#'     enzyme of the class. Range 0-2, substrate-agnostic. Distinct in kind, not
#'     statistically independent.}
#'   \item{evidence_tier}{"Gold" / "Silver" / "Bronze" / "Provisional". Practical
#'     prioritization heuristic over the two evidence classes plus supplementary
#'     GO/UniProt-keyword support. Not a probability, not an evidence count.}
#'   \item{curated_core}{Logical. TRUE if the gene has >= 1 evidence dimension
#'     (i.e. not Provisional / not comprehensive-only). Filtered by the accessor
#'     `mode = "strict"` knob, which never touches the substrate booleans.}
#'   \item{in_structural_catalog}{Logical. Axis 1: Manning / KinHub / kinase.com /
#'     pkinfam.}
#'   \item{is_protein_kinase_ec}{Logical. A protein-specific kinase EC
#'     (2.7.10-2.7.14) fired. The protein-only EC flag -- narrower than the EC rigor
#'     dimension counted in `n_evidence_dimensions`, which credits any class kinase EC
#'     (protein or small-molecule), so the two diverge for non-protein and dual EC
#'     kinases such as PIK3CA.}
#'   \item{membership_basis}{Deriving source of the Axis-1 (structural-catalog) call:
#'     "reconstructed:pkinfam" (cleanly-licensed anchor), "reconstructed:kinase.com",
#'     or "crosscheck:KinHub"; NA when the gene is in no structural catalog.}
#'   \item{go_experimental}{Logical. A non-electronic GO kinase-activity code is
#'     present (experimental/curated support).}
#'   \item{has_uniprot_kw}{Logical. Reviewed Swiss-Prot kinase keyword present.}
#'   \item{supplementary_support}{Logical. `go_experimental` OR `has_uniprot_kw`;
#'     splits Silver from Bronze among single-axis genes.}
#'   \item{kinase_family}{Manning family (catalog-pure); may be NA.}
#'   \item{classification_reason}{Human-readable rationale for the substrate call.}
#'   \item{kinase_group}{Manning group (AGC, CMGC, ...); may be NA.}
#'   \item{kinase_subfamily}{Manning subfamily; may be NA.}
#'   \item{derived_family}{Non-Manning family descriptor for genes lacking a
#'     `kinase_family` (from the UniProt family string or the GO class).}
#'   \item{uniprot_protein_family}{Raw UniProt "Protein families" string.}
#'   \item{dual_protein_nonprotein}{Logical. A protein kinase that also has a
#'     non-protein kinase activity (e.g. PI3K family, NME1, DGKQ).}
#'   \item{n_membership_sources}{Integer 0-7. Count of the per-source `is_*`
#'     flags; informational breadth-of-support (not the rigor metric).}
#'   \item{is_pseudogene}{Logical. HGNC locus type matches "pseudogene".}
#'   \item{hgnc_id}{HGNC identifier.}
#'   \item{gene_name}{HGNC approved gene name.}
#'   \item{entrez_id}{NCBI Entrez gene ID.}
#'   \item{uniprot_ids}{UniProt accession(s).}
#'   \item{prev_symbol}{HGNC previous symbol(s).}
#'   \item{alias_symbol}{HGNC alias symbol(s).}
#'   \item{enzyme_id_EC}{EC number(s) from HGNC `enzyme_id`.}
#'   \item{ec_kinase_subclass}{Matched EC 2.7 kinase subclass(es).}
#'   \item{hgnc_gene_group}{HGNC gene-group membership string.}
#'   \item{locus_type}{HGNC locus type.}
#'   \item{chromosomal_location}{Cytogenetic location.}
#'   \item{mane_select_transcript}{MANE Select transcript.}
#'   \item{iuphar_id}{IUPHAR/Guide to Pharmacology identifier.}
#'   \item{hgnc_kinase_gene_group}{Logical. HGNC gene-group names a kinase group
#'     (excluding non-catalytic terms).}
#'   \item{is_pkinfam}{Logical. Member of UniProt pkinfam.}
#'   \item{is_manning}{Logical. Member of the Manning / kinase.com kinome.}
#'   \item{is_kinhub}{Logical. Member of the KinHub kinome.}
#'   \item{is_go_kinase_activity}{Logical. In the GO kinase-activity umbrella.}
#'   \item{is_ec_kinase}{Logical. Carries an EC 2.7 kinase number (any subclass).}
#'   \item{is_idg_dark_kinase}{Logical. In the IDG understudied ("dark") kinome.}
#' }
#' @source Manning et al. Science 2002;298:1912-1934; KinHub/KinMap (Eid 2017);
#'   kinase.com; UniProt pkinfam; GO; IUBMB EC. See the package vignette.
"human_kinases"

#' Human protein-phosphatase reference table
#'
#' One row per human phosphatase gene, keyed on base Ensembl gene ID and typed by
#' substrate class. Parallel in schema to [human_kinases]. Built by
#' `inst/scripts/make-data.R` from pinned source snapshots; see the vignette.
#'
#' @format A data frame with one row per gene and the following columns:
#' \describe{
#'   \item{ensembl_gene_id}{Base (unversioned) Ensembl gene ID. Primary key.}
#'   \item{symbol}{HGNC symbol.}
#'   \item{acts_on_protein}{Logical. TRUE if the enzyme dephosphorylates protein substrates.}
#'   \item{acts_on_nonprotein}{Logical. Co-equal with `acts_on_protein` (a dual enzyme such as
#'     PTEN is TRUE for both); filtered orthogonally by the `substrate` accessor knob.}
#'   \item{nonprotein_substrate_type}{Pipe-delimited non-protein classes ("lipid", "nucleotide",
#'     "carbohydrate", "metabolite", "other"); empty for protein-only.}
#'   \item{substrate_subtype}{Coarse substrate-type bucket derived from
#'     `substrate_call` and `nonprotein_substrate_type`: one of "Protein phosphatase",
#'     "Lipid phosphatase", "Nucleotide phosphatase", "Carbohydrate/sugar phosphatase",
#'     "Other small-molecule phosphatase", or "Unclassified phosphatase". For the finer
#'     non-protein class see `nonprotein_substrate_type`.}
#'   \item{substrate_call}{One of "protein" / "nonprotein" / "dual" / "untyped" -- the collapsed
#'     substrate verdict from `acts_on_protein` / `acts_on_nonprotein`.}
#'   \item{substrate_evidence}{Which evidence kinds drove the substrate call, joined with "+":
#'     any of "GO", "Chen", "EC".}
#'   \item{substrate_concordance}{"concordant", "single", or "untyped".}
#'   \item{substrate_decider}{The precedence-winning evidence that set the substrate call
#'     ("GO-exp", "GO-elec", "Chen-flag", "EC", or "precedence-default").}
#'   \item{ec_protein}{Logical. A protein-substrate EC rule fired.}
#'   \item{ec_nonprotein}{Logical. A non-protein-substrate EC rule fired.}
#'   \item{go_protein}{Logical. A protein-substrate GO phosphatase-activity term fired.}
#'   \item{go_nonprotein}{Logical. A non-protein-substrate GO phosphatase-activity term fired.}
#'   \item{dual_protein_nonprotein}{Logical. Acts on both protein and non-protein substrates.}
#'   \item{catalytic_status}{"active" / "pseudo" / "uncertain", from Chen 2017.}
#'   \item{is_catalytic_background}{Logical. Substrate-blind enrichment background:
#'     `catalytic_status == "active"` AND `curated_core`, spanning protein and non-protein alike.}
#'   \item{is_protein_catalytic_background}{Logical. The default protein-phosphatase enrichment
#'     background: `catalytic_status == "active"` AND `curated_core` AND `acts_on_protein`.}
#'   \item{is_pseudophosphatase}{Logical. Predicted catalytically dead (Chen 2017).}
#'   \item{regulates}{Pipe-delimited symbol(s) of the catalytic phosphatase(s) this gene
#'     regulates as an adapter/activator; NA for non-regulatory genes. Populated for the
#'     catalytically inactive myotubularins, which stay in the catalytic master as untyped
#'     phosphatome members rather than being relocated to a separate companion.}
#'   \item{regulatory_role}{Free-text description of the regulatory relationship and its
#'     primary citation; NA when `regulates` is NA.}
#'   \item{n_evidence_dimensions}{Integer. Count of distinct class-evidence kinds -- structural/
#'     sequence catalog and any class-specific EC -- confirming a bona fide enzyme of the class.
#'     Range 0-2, substrate-agnostic. Distinct in kind, not statistically independent.}
#'   \item{evidence_tier}{"Gold" / "Silver" / "Bronze" / "Provisional"; heuristic over the two
#'     evidence classes plus supplementary support.}
#'   \item{curated_core}{Logical. >= 1 evidence dimension. Filtered by the accessor
#'     `mode = "strict"` knob, which never touches the substrate booleans.}
#'   \item{in_structural_catalog}{Logical. Axis 1: Chen 2017 phosphatome or an HGNC
#'     protein-phosphatase gene group.}
#'   \item{is_protein_phosphatase_ec}{Logical. A protein-specific phosphatase EC
#'     (3.1.3.16 Ser/Thr or 3.1.3.48 Tyr) fired. The protein-only EC flag -- narrower
#'     than the EC rigor dimension counted in `n_evidence_dimensions`, which credits any
#'     class phosphatase EC, so the two diverge for non-protein and dual EC phosphatases.}
#'   \item{go_experimental}{Logical. Non-electronic GO phosphatase-activity support.}
#'   \item{has_uniprot_kw}{Logical. Reviewed Swiss-Prot keyword KW-0904 present.}
#'   \item{supplementary_support}{Logical. `go_experimental` OR `has_uniprot_kw`.}
#'   \item{membership_basis}{Deriving source of the Axis-1 call ("reconstructed:Chen2017" /
#'     "reconstructed:HGNC_groups"); NA when in no structural catalog.}
#'   \item{classification_reason}{Human-readable rationale for the substrate call.}
#'   \item{phosphatase_fold}{Chen structural fold (CC1, HAD, PPM, ...); may be NA.}
#'   \item{phosphatase_family}{Chen family; may be NA.}
#'   \item{phosphatase_subfamily}{Chen subfamily; may be NA.}
#'   \item{n_membership_sources}{Integer. Count of TRUE per-source flags.}
#'   \item{is_pseudogene}{Logical. HGNC locus type matches "pseudogene".}
#'   \item{hgnc_id}{HGNC identifier.}
#'   \item{gene_name}{HGNC approved gene name.}
#'   \item{entrez_id}{NCBI Entrez gene ID.}
#'   \item{uniprot_ids}{UniProt accession(s).}
#'   \item{prev_symbol}{HGNC previous symbol(s).}
#'   \item{alias_symbol}{HGNC alias symbol(s).}
#'   \item{enzyme_id_EC}{EC number(s) from HGNC `enzyme_id`.}
#'   \item{hgnc_gene_group}{HGNC gene-group membership string.}
#'   \item{locus_type}{HGNC locus type.}
#'   \item{chromosomal_location}{Cytogenetic location.}
#'   \item{mane_select_transcript}{MANE Select transcript.}
#'   \item{iuphar_id}{IUPHAR/Guide to Pharmacology identifier.}
#'   \item{is_chen}{Logical. In the Chen 2017 protein phosphatome.}
#'   \item{is_hgnc_phosphatase_group}{Logical. In an HGNC protein-phosphatase gene group.}
#'   \item{is_go_phosphatase_activity}{Logical. In the GO phosphatase-activity umbrella.}
#'   \item{is_phosphatase_ec}{Logical. Carries an EC 3.1.3 number (any 4-digit).}
#'   \item{is_uniprot_kw_phosphatase}{Logical. Reviewed UniProt KW-0904.}
#' }
#' @source Chen, Dixon, Manning, Sci Signal 2017 (phosphatome.net); HGNC gene groups;
#'   UniProt KW-0904; GO; IUBMB EC. See the package vignette.
"human_phosphatases"

#' Unified human phospho-enzyme summary
#'
#' Thin, class-agnostic table spanning kinases and phosphatases, derived from [human_kinases]
#' and [human_phosphatases]. Shared substrate/evidence columns only; join to a master by
#' `ensembl_gene_id` to recover class-specific taxonomy.
#'
#' @format A data frame with one row per gene and the following columns:
#' \describe{
#'   \item{ensembl_gene_id}{Base Ensembl gene ID. Primary key.}
#'   \item{symbol}{HGNC symbol.}
#'   \item{regulator_class}{"kinase" or "phosphatase".}
#'   \item{acts_on_protein}{Logical. Acts on protein substrates.}
#'   \item{acts_on_nonprotein}{Logical. Acts on a non-protein substrate.}
#'   \item{nonprotein_substrate_type}{Pipe-delimited non-protein classes; empty = protein-only.}
#'   \item{dual_protein_nonprotein}{Logical. Acts on both.}
#'   \item{catalytic_status}{"active" / "pseudo" / "uncertain".}
#'   \item{n_evidence_dimensions}{Integer. Count of distinct class-evidence kinds -- structural/
#'     sequence catalog and any class-specific EC -- confirming a bona fide enzyme of the class.
#'     Range 0-2, substrate-agnostic. Distinct in kind, not statistically independent.}
#'   \item{evidence_sources}{Semicolon-joined names of the supporting sources.}
#'   \item{evidence_tier}{"Gold" / "Silver" / "Bronze" / "Provisional".}
#'   \item{curated_core}{Logical. >= 1 evidence dimension. Filtered by the accessor
#'     `mode = "strict"` knob.}
#'   \item{substrate_call}{One of "protein" / "nonprotein" / "dual" / "untyped".}
#'   \item{is_catalytic_background}{Logical. Substrate-blind: `catalytic_status == "active"` AND
#'     `curated_core`.}
#'   \item{is_protein_catalytic_background}{Logical. The default protein-enzyme enrichment
#'     background: active AND `curated_core` AND `acts_on_protein`.}
#' }
#' @source Derived from [human_kinases] and [human_phosphatases].
"human_phosphoenzymes"
