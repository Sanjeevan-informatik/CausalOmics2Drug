# Input contract

Use **Data workspace** to upload one of:

1. CSV/TSV with normalized evidence rows for any of the seven layers.
2. A dataset JSON with evidence, pathway gene sets, measured-gene universe and notes.
3. CellOmics expression CSV with `gene,tumor_tpm,control_tpm`.

An upload replaces the current in-memory dataset after successful server validation. Export the dataset or report before closing the page. No data is written to a server database.

## Required evidence fields

| Field | Meaning |
|---|---|
| id | Unique evidence record identifier |
| gene | Uppercase gene symbol, consistently mapped to a single organism/build upstream |
| layer | genomics, epigenomics, bulk_rna, single_cell, spatial, proteomics, metabolomics |
| disease, tissue | Exact context labels; use consistent ontology IDs/labels in a real study |
| study, source | Study ID and reference/file/accession providing this observation |
| observation | Human-readable description of the measured quantity and contrast |
| strength | Finite, supplied value in [0,1]; calibration is the uploader's responsibility |

Optional: `effect` (signed finite number), `unit`, `q_value` ([0,1]), `cell_type` (default all), `kind` (default association), `action` (default unknown), `synthetic` (default false), `x` and `y` (must occur together). Missing values are empty in CSV or null in JSON. Unknown columns are rejected to catch mapping errors.

`kind`: association / fine_mapping_coloc / perturbation / mendelian_randomization. `action`: inhibit / activate / unknown. Only causal-type rows with sufficient strength and no failed supplied q-value contribute intervention direction. These fields record upstream interpretation; selecting a label cannot validate an experiment.

See `data/evidence_template.csv` and `data/demo.json` for complete examples. The template contains synthetic rows on purpose. Replace all example context, observations, provenance and scores when using real evidence.

## Layer preparation

| Source assay | Required upstream work before upload |
|---|---|
| Genomics | QC and variant calling; functional annotation and defensible variant-to-gene mapping. Use copy number only as association unless independent causal analysis exists. |
| Epigenomics | QC, normalization, differential regions and regulatory-region-to-gene mapping. Methylation direction alone does not determine expression. |
| Bulk RNA | Count QC, replicated disease/control model, batch/covariate adjustment and multiple-testing correction. TPM ratios are descriptive only. |
| scRNA | Cell QC and annotation, aggregate raw counts by **donor and cell type**, then model biological replicates; cells are not independent patients. |
| Spatial | Assay-specific normalization, tissue/spot mapping, donor and spatial-dependence-aware inference; upload coordinates in a consistent system. |
| Proteomics | Peptide/protein QC, abundance/phosphosite normalization, missingness handling and explicit protein-to-gene mapping. |
| Metabolomics | Feature identification confidence, normalization and multiple testing; map metabolites to enzymes/pathways explicitly and label indirect evidence. |

The CellOmics expression adapter computes `log2((tumor_tpm+1)/(control_tpm+1))` and illustrative strength `min(abs(effect)/4,1)`. It cannot estimate variance, q-values or causality from one aggregate value per group.

## Pathway JSON

`pathways` is a list of `{id, name, genes, source}` objects. `universe` contains **all genes that could have been selected by the upstream experiment**, including non-significant measured genes, not just pathway members. It must include all evidence genes. Use a context-appropriate universe; if disease/tissue assays measure different backgrounds, analyze separate dataset files. Duplicate pathway IDs and evidence IDs are rejected. Duplicate gene names in a gene set/universe are treated as sets.

Pathway hits require both strength ≥ the minimum and a supplied q-value ≤ the threshold. Multiple assays selecting the same gene count as one hit. Provide complete, versioned gene sets for real analyses; the demo sets are deliberately small teaching subsets.

## API

- `GET /api/health`: health/version.
- `GET /api/demo`: complete typed teaching dataset.
- `POST /api/validate`: validate a dataset JSON and return all default fields.
- `POST /api/import`: `{text, name, disease, tissue}` → validated dataset.
- `POST /api/analyze`: `{dataset, disease, tissue, cell_type?, min_strength?, max_q?, weights?}` → targets, pathways, provenance, warnings, parameters and input SHA-256.

Limits: 8 MB request body, 20,000 evidence rows, 1,000 pathway sets, 5,000 genes/set and 60,000 universe genes. This contract targets summarized evidence, not raw single-cell matrices.
