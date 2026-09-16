# AlloTrace: transplant immunogenomics research prototype

AlloTrace is the transplant workspace within CausalOmics2Drug. It follows donor–recipient germline differences through antigen presentation, T-cell recognition, functional evidence and longitudinal immune escape. It is independent software, not an official Charité or Penter Lab product. The original seven-omics workspace remains at `/#omics`.

## Scientific model

A recipient-only allele at a callable donor–recipient germline locus is a **directional mismatch**, not an established minor histocompatibility antigen. This distinction is motivated by the mismatch/presentation separation in [AlloPipe, Dhuyser et al., 2026](https://doi.org/10.1111/tan.70590) and systematic miHA investigation in [Cieri et al., Nature Biotechnology](https://doi.org/10.1038/s41587-024-02348-3). The implementation is an original, simplified evidence-review model; it does not reproduce either paper's pipeline or validated outcome model.

The [Penter Lab research programme](https://www.penterlab.org/projects/) and [Penter et al. donor/recipient tracking work](https://doi.org/10.1158/2643-3230.BCD-23-0138) motivate separate sample origin, clonotype, compartment, time and leukemia-evolution records. This application does not infer cellular origin from mtDNA. Supplied origin labels require an independent method and source.

### Candidate gates and ordering

1. Donor and recipient calls must be present, with depth ≥20 and GQ ≥30 by default. Missing calls are unknown; they never prove antigen absence.
2. The supplied target allele (REF or ALT) must occur in the recipient and not in the donor. Reverse mismatches are flagged separately.
3. Restricting HLA-I must appear in both supplied typing lists. Only two-field HLA-A/B/C is supported. Recipient germline provenance cannot be unknown. Remission material still requires review for contamination and residual disease.
4. Imported EL percentile ≤2 and leukemia expression ≥1 TPM support a prediction tier. These are adjustable screening thresholds, not probabilities or validated clinical cutoffs. EL percentile is distinct from affinity in nM; [NetMHCpan documentation](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) explains those output types. No predictor executes in the browser.
5. Imported immunopeptidomics detection can support a presentation tier. A controlled tetramer positive can support binding evidence. Functional tier requires positive re-expression, leukemia killing and HLA-dependence records converging on the **same paired, productive receptor**, with ≥2 replicates and adequate controls. An adequate negative in the same assay category prevents that functional tier. These minimum evidence rules are transparent heuristics, not statistical proof.
6. Ranking orders eligibility first, then positive safety flags, evidence stage, count of eight evidence checks, EL rank and identifier. The eight-check count is completeness, not causal probability. Normal-tissue expression above the selected threshold is flagged, not proof of toxicity. Missing normal-tissue data is not safety. A positive normal-tissue or cross-reactivity record is a safety-review signal even when preliminary.

No score establishes clinical eligibility. [HA-1 TCR therapy research](https://doi.org/10.1182/blood.2024024105) illustrates a translational endpoint; the software's synthetic peptides are not HA-1, are not taken from that trial, and are not therapeutic constructs.

### Receptor dynamics and experimental evidence

Exact paired amino-acid CDR3 plus V/J defines imported single-cell clonotypes, within a donor–recipient pair. Beta-only receptors remain separate; multiple productive chains are marked ambiguous and represented for review. Sequence similarity does not establish antigen specificity, consistent with the QC emphasis of [Scirpy](https://scirpy.scverse.org/en/latest/tutorials/tutorial_3k_tcr.html). The tabular importer supports a documented subset of [AIRR](https://docs.airr-community.org/) fields; it does not claim full AIRR conformance.

Frequency = clone count / supplied sample denominator. Comparisons are restricted to the same assay and compartment. Missing observation rows remain missing and break chart lines. Explicit zero is a measured nondetection. Fold change is undefined for a zero baseline; percentage-point change is still available. Wilson 95% binomial intervals use z=1.95996398454. These intervals describe counting uncertainty only; bulk template PCR bias, library preparation, biological variability and repeated sampling are not modeled. No significance claim is made.

Tetramer event counts and negative-control counts remain visible. Assay categories distinguish positive, negative, conflicting, preliminary and missing evidence. The application records already-performed experiments; it neither conducts cell engineering nor prescribes clinical protocols. Shortlisting is manual for review and does not constitute suitability for infusion or a clinical study.

### Immune escape

Blast percentage, chimerism, antigen expression and somatic VAF are separate selectable longitudinal measurements. HLA loss/reduction labels are imported observations, not inferred from expression alone. VAF is not cancer-cell fraction: purity, copy number and zygosity affect it. Temporal coexistence does not establish immune selection. Samples with multiple somatic clones must be inspected in the table; no phylogeny is reconstructed.

## A reproducible upstream workflow

1. Retain coded donor/recipient identities and specimen provenance in approved institutional storage. Select independent germline material; do not treat post-transplant mixed marrow as uncontaminated recipient germline.
2. Run institution-validated WES alignment, duplicate/read QC, joint variant calling, contamination/relatedness and sex/chromosome-aware QC. Match reference assembly and normalize variants. Review HLA typing independently.
3. Annotate transcript/protein consequences with a pinned annotation release. Reconstruct both donor and recipient peptides, accounting for phase and adjacent variants upstream. The browser supports only supplied single-missense, biallelic SNV peptide pairs; indels, structural variants, noncanonical translation and complex haplotypes are outside scope.
4. Run licensed or institution-approved presentation prediction externally, preserving tool version, allele, EL rank, affinity and parameters. Record quantitative leukemia and normal-tissue expression with assay provenance. Transcript expression is neither peptide abundance nor assurance of safety.
5. Process scRNA/VDJ and bulk TCR data with established QC, doublet handling and assay-specific denominators. Import summarized clone/count tables or the supported single-cell AIRR subset. Keep bulk and single-cell measurements distinct.
6. Import experimental binding and functional evidence linked to the exact candidate, receptor and sample. Review controls, replicate independence, conflicting data and off-tissue findings. Export the complete study and report with settings and SHA-256 input fingerprint.

## Data and validation boundaries

All bundled examples are wholly synthetic, generated by `scripts/build_transplant_demo.py`; no patient-derived sequence or outcome is represented as a discovery. The software has no trained antigen-specificity model, performs no raw FASTQ/BAM processing and does not make causal therapeutic-target claims from correlations. Original multi-omics evidence exploration is preserved, but not silently merged with miHA scoring. No cell-omics-explorer dataset was supplied or incorporated in this update.

TypeScript is used for deterministic in-browser analysis and validated React interfaces; the existing Python FastAPI analysis backend is retained for the original workspace. Switching language alone does not improve biological inference. Tests cover missing genotypes, allele direction, conflicting functional evidence, pairing, totals, Wilson intervals and importer failures. External patient cohorts and prospective experimental validation are still required before drawing biological conclusions.
