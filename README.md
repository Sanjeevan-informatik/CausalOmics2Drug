# CausalOmics2Drug

**An explainable multi-omics workbench for therapeutic target hypotheses.**

Integrate genomics, epigenomics, bulk RNA, single-cell RNA, spatial transcriptomics, proteomics and metabolomics evidence in a single disease/tissue context. Inspect pathways, rank candidate targets, and trace every result to its source.

**Version 0.1 is a working research prototype for analyzed evidence tables.** Its scores summarize supplied evidence; they do not estimate causal effects or establish drug efficacy. The included EGFR example is explicitly teaching data. No validated drug is produced by the software.

## Run locally

Requires Python 3.11+ and Node.js 22.12+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
npm ci --prefix frontend
npm run build --prefix frontend
python -m uvicorn backend.app.main:app --host 127.0.0.1 --port 8000
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1` instead.

Open **http://127.0.0.1:8000**. API documentation: **http://127.0.0.1:8000/docs**. The FastAPI process serves the built React application and API together.

For development, start the Python API and run `npm run dev --prefix frontend` in a second terminal. Vite proxies `/api` to port 8000.

Or build the container:

```bash
docker build -t causalomics2drug .
docker run --rm -p 127.0.0.1:8000:8000 causalomics2drug
```

## What works

| Workspace | Function |
|---|---|
| Discovery | Target ranking, seven-layer coverage, gene search and context selectors |
| Data workspace | Validated CSV/TSV evidence imports, full dataset JSON and CellOmics expression CSV adapter |
| Target evidence | Per-layer contributions, missingness, source ledger, effect values, q-values, cell context and spatial coordinate plot |
| Pathway story | Explanatory EGFR branching diagram, observations-to-experiment walkthrough and active-dataset pathway overlap |
| Methods | Adjustable weights and thresholds, documented scoring, deterministic input fingerprint and full report export |

All seven modalities use one normalized evidence contract. Integration aligns evidence by **gene + disease + tissue + cell context**, retaining study, source and modality. It is evidence-level integration, not a joint latent model or automatic sample matching. Metabolite-to-enzyme and regulatory-region-to-gene mappings must be supplied explicitly. Rows are not inferred to come from the same patient.

## The CellOmics connection

Inspected and preserved four examples from [cell-omics-explorer/data](https://github.com/Sanjeevan-informatik/cell-omics-explorer/tree/main/data). The EGFR teaching narrative uses:

- Copy number approximately 6, high promoter accessibility, RNA 320 versus 30 TPM.
- Protein 2.8× control, phosphorylation at Y1068 3.4× control, lactate 2.2× control.
- A simplified EGFR → RAS/RAF/ERK and EGFR → PI3K/AKT/mTOR signaling story, linked to Reactome.

The source files have no disease subtype, donor replication or adjusted p-values. Generated single-cell, spatial and downstream pathway rows are marked as teaching extensions. The original six-cell immune-marker matrix has no EGFR measurements or donor identifiers and is preserved without being falsely joined to the tumor example. See [provenance](docs/PROVENANCE.md).

## Scoring and scientific limits

Omics score = `100 × sum(weight[layer] × max qualifying strength[layer]) / sum(weights)`.

- Strength is an explicitly supplied, uncalibrated number between 0 and 1. Different study types need scientifically justified upstream calibration.
- Missing layers contribute zero to this *evidence-support* score and remain visibly missing; zero is not biological evidence against a target.
- A supplied q-value above the threshold excludes that row from scoring. Missing q-values are permitted only as descriptive support.
- Repeated rows cannot increase the maximum for a layer. This is conservative deduplication, not a complete model of study dependence.
- Causal support is shown separately using three supplied evidence types: fine-mapping/colocalization, perturbation and Mendelian randomization. The software does not run those methods or check their assumptions. A high score is not proof of causality.
- Direction of intervention is shown only when qualifying causal-type evidence supplies it; conflicting directions are exposed.
- Ranking orders by causal support, then omics support. All scores, parameters and source rows are exported.

Pathway over-representation uses hypergeometric upper-tail probabilities and Benjamini–Hochberg adjustment over tested sets. It requires an explicitly supplied measured-gene universe and genes with supplied qualifying q-values. The teaching demo has neither and reports descriptive overlap without p-values. ORA can be biased by feature-selection and assay detection; it does not determine pathway activation or disease specificity. For unbiased disease-specific inference, use appropriate controls and independently validated disease/tissue contrasts.

Read [input format](docs/INPUTS.md) and [scientific design](docs/SCIENCE.md).

## Tests

```bash
pip install -r requirements-dev.txt
python -m pytest backend/tests -q
npm run build --prefix frontend
cd frontend
npx playwright install chromium
npm test
```

The browser suite starts the API itself. Stop any existing process on port 8000 first. Tests cover data import, duplicate invariance, contextual isolation, missing data, conflicting causal directions, FDR exclusion, known enrichment probabilities, deterministic exports, responsive layout and user navigation.

## Architecture and next research steps

Python/FastAPI/Pydantic/SciPy provide typed validation and auditable numerical analysis. TypeScript/React provides a responsive interface with typed view models. No language guarantees better science; this combination matches the available scientific tooling and keeps the analysis testable.

Raw FASTQ/VCF/H5AD/mzML processing, donor-level differential-expression models, fine-mapping/colocalization, MR, trained multi-view integration, live Open Targets/drug database queries, compound design and clinical decision support are **not implemented**. A defensible next study should select a disease and independent patient cohorts; prepare replicated contrasts with batch/confounder handling; validate target ranks against held-out perturbations; and compare against simple baselines before adding complex AI models. See the scientific design for concrete acceptance criteria.

This application is stateless and has no authentication or persistent database. Keep it on localhost for research development. A shared patient-data service would need an explicitly designed access and storage model. No public website has been deployed by this repository initialization.
