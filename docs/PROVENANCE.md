# Data provenance

Inspected on 2026-09-15 from the user's `Sanjeevan-informatik/cell-omics-explorer` repository. Original file contents are preserved under `data/cell_omics/`.

| File | Git blob SHA | Use |
|---|---|---|
| example_expression.csv | 80e3d7da9328553e317c3761dae0edf708e299e5 | Five descriptive TPM contrasts; +1 pseudocount log2 ratios |
| example_multiomics.tsv | 85435419dde5ecf09b0174d16712b329c64e9cfe | EGFR multi-layer and lactate teaching observations |
| example_methylation.tsv | 813c4baa2f5d218748155fae30cae86cc25b6c93 | Four methylation values without controls; descriptive context |
| example_single_cell.csv | a97a3218355d6ad8b484b496594bbd44a4b4e746 | Original immune-marker toy matrix; preserved, not integrated into tumor example |

Source base URL: https://github.com/Sanjeevan-informatik/cell-omics-explorer/tree/main/data . Blob SHAs identify the inspected versions even if main changes.

Run `python scripts/build_demo.py` to reproduce `data/demo.json` and the CSV template. Added EGFR single-cell/spatial entries, downstream pathway protein context and the LDHA mapping are generated teaching extensions. All demo records are marked `synthetic: true`; no real patient measurements or inferred disease subtype are asserted. No q-values, causal support records or clinical drug associations were invented for the demo.

The source teaching values and extensions are distinguished in each record's `source` and `observation`. Hand-set strengths are visible in the generator and dataset. They are UI/method demonstrations, not trained coefficients or empirical confidence estimates.
