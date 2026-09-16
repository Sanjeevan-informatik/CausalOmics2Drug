# Inputs, exports and provenance

Start with the downloadable synthetic study JSON. It is the complete schema example, including missing data and negative controls. A whole-study JSON replaces the active in-memory study only after validation and explicit Apply. CSV/TSV tables replace one table at a time and must preserve references to all other tables. Use whole-study JSON when simultaneously changing related identifiers.

## Tables

| Table | Meaning and important fields |
|---|---|
| pairs | Coded pair/donor/recipient IDs; GRCh38; donor_hla and recipient_hla arrays; germline_source; source |
| samples | Pair, day relative to transplant, blood/marrow/graft/skin, bulk_tcr/single_cell, total denominator, unit, independently supplied origin and origin_method |
| candidates | Annotated biallelic missense variant, both GT/DP/GQ, REF/ALT target, 8–11 residue peptide pair and one-based variant_index, HLA-I, imported prediction and expression, provenance |
| clones | TRA/TRB CDR3 amino acid and V/J, paired/beta_only/ambiguous, productive, source |
| observations | Unique sample_id/clone_id, integer count, supplied phenotype, source |
| assays | Candidate/receptor/sample linkage, kind, outcome, controls, replicates, optional paired event counts/denominators, source and notes |
| escape | Sample, blast_percent, donor_chimerism, supplied hla_status and allele, antigen_tpm, somatic_clone, vaf 0–1, source |

CSV HLA arrays use semicolons. Numeric missing values are blank; zero remains zero. GT no-call is `./.`. Only normalized diploid biallelic genotype forms are accepted; haploid sex-chromosome calls require upstream handling and are not coerced. Source text is required and should identify a controlled accession, software version, assay run or analysis output. Do not put names or direct identifiers into free text.

The importer rejects unknown columns, wrong types, duplicate IDs, cross-pair assay links, duplicated observations, impossible count totals and inconsistent peptide annotations. File size is limited to 15 MB: this is for analyzed research tables, not sequencing matrices. Sample denominators must be the relevant repertoire population, not total unrelated cells. Missing clone rows are missing measurements, not assumed absences.

## Annotated VCF subset

The example file is downloadable in Data & provenance. It requires `##reference=GRCh38`, a normal `#CHROM` header and two explicitly selected sample columns. Only PASS biallelic SNVs are considered. Required INFO: GENE, TRANSCRIPT, PEP_REF, PEP_ALT, PEP_POS, HLA. Optional: TARGET (REF or ALT), RANK_EL, AFF_NM, LEUK_TPM, HEM_TPM, NORMAL_TPM, PREDICTOR, SOURCE. FORMAT fields GT, DP and GQ are retained. Missing values stay unknown. Unsupported records are counted and reasons displayed. Accepted IDs include pair and source line number.

Raw VCF annotation and peptide prediction are NOT performed. This custom INFO contract is not a universal VEP/AlloPipe export format. Prepare a reviewed adapter upstream. Applying VCF replaces that pair's candidates and removes their previously linked assays; the preview reports that removal before Apply.

## Single-cell AIRR subset

Required TSV fields: sequence_id, cell_id, sample_id, locus, v_call, j_call, junction_aa, productive. Samples must already exist as single_cell. Counts are unique cells grouped by exact paired receptor definition; one TRA plus one TRB is paired, TRB alone beta_only, multiple productive chains ambiguous. Nonproductive and non-TRA/TRB rows are excluded; alpha-only cells are reported as orphan. Full original chains must remain upstream because ambiguous cells retain a representative only.

Applying AIRR replaces the entire receptor/count tables and clears previous linked assays to prevent accidental reassignment. The preview explicitly reports this. Bulk repertoires use normalized clone/observation CSVs; bulk sequences cannot infer alpha–beta pairing.

## Saving

The app does not autosave to a server or browser storage. Export study JSON to retain input edits; it can be reimported. Export report includes full input, settings, selected-pair results, manual shortlist, references and session events. A report is a review artifact, not a study import: use its `input.project` member if reconstructing input. Its SHA-256 fingerprint covers canonicalized project and settings; it does not certify source authenticity. Session events are not an immutable audit log. Exported CSVs escape formula-like values; resulting quoted text may require review before round-trip import.
