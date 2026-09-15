# Scientific design and validation plan

## Research question

Which genes have context-matched, convergent evidence sufficient to justify a targeted validation experiment? Keep three distinct questions visible: association with disease, causal involvement in phenotype, and whether a selective intervention is feasible and beneficial. The current implementation addresses evidence organization and heuristic prioritization; it does not solve causal identification.

## Defensible causal evidence

Human genetic evidence needs credible fine-mapping, appropriate locus-to-gene assignment and colocalization. MR requires valid instruments, instrument-strength checks, sensitivity analyses and attention to horizontal pleiotropy. A perturbation needs biological replicates, appropriate controls, multiple independent reagents, off-target assessment and ideally rescue/orthogonal intervention. Evidence types recorded in the UI are claims supplied by an analyst, not automatically verified results.

Multiple evidence types can arise from the same underlying dataset. The score does not model this dependence and cannot be interpreted as a posterior probability. Equal modality weights and maximum aggregation are transparent starting heuristics, not optimized statistical estimators. Do not call a candidate disease-specific merely because it appears in one filtered context.

## EGFR example

The inspected teaching table contains EGFR amplification, accessibility, RNA, protein and phosphorylation observations. These motivate a receptor-signaling hypothesis. Reactome describes connections from EGFR to intracellular pathways controlling growth and survival. The application compresses intermediates in an explanatory diagram; it is not a new mechanistic network reconstruction. Lactate abundance is indirect context, and assigning its change to a particular enzyme or signaling driver requires further work.

The source data lacks subtype, donors and controls sufficient for inference. Therefore no approved-drug match, personalized intervention, differential-expression p-value or causal claim is fabricated. New spatial and single-cell demo observations exist only to exercise the data model and are marked synthetic.

## Study design for a real benchmark

1. Predefine one disease/subtype, relevant normal comparator, tissue, inclusion criteria and a measurable phenotype.
2. Assemble appropriate public or authorized cohorts. Record accessions, assay versions, sample/donor relationships, consent restrictions and preprocessing versions. Keep incompatible cohorts separate.
3. Use assay-specific preprocessing and replicated statistical contrasts. Aggregate scRNA by donor/cell type. Address batch, cell composition, confounders, spatial dependence and missingness before evidence integration.
4. Add context-matched genetic and perturbational evidence, with independent validation. Record contradictions instead of discarding them.
5. Partition train/tuning/validation by donor and study before any learned normalization or fitting. For future predictive models, hold out targets/families or entire studies as appropriate to prevent target leakage.
6. Compare against differential-expression-only, equal-weight and single-evidence-type baselines. Report PR-AUC or top-k precision against predeclared perturbation endpoints, rank stability under ablation/bootstrap, coverage and failure cases. Report calibration only if a probabilistic model is trained and validated.
7. Treat pharmacological tractability, on-target normal-tissue effects, chemical selectivity and exposure as separate evidence axes. Do not score absent safety information as safe.
8. Validate shortlisted targets experimentally. Only then proceed to hit identification/optimization and appropriate preclinical and clinical development.

## Method sources

- [Open Targets evidence documentation](https://platform-docs.opentargets.org/evidence): target–disease evidence organization and source-specific scoring.
- [Reactome EGFR signaling](https://reactome.org/content/detail/R-HSA-177929) and [EGFR signaling in cancer](https://reactome.org/content/detail/R-HSA-1643713): curated mechanistic context.
- [Reactome PI3K/AKT signaling in cancer](https://reactome.org/content/detail/R-HSA-2219528): pathway context.
- [Single-cell atlas protocol: differential expression](https://atlas-protocol.readthedocs.io/): biological-replicate pseudobulk analysis.

Consulted 2026-09-15. This repository does not bundle a full Reactome release or query Open Targets live. Gene-set source URLs and evidence URLs are references, not an automated provenance verification service.
