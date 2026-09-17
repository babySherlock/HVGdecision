# HVGDecision

HVGDecision refines HVG panels for integration of multi-donor single-cell datasets from comparable biological settings. It combines donor-level and cell-level risk evidence with biological protection, records each gene decision, and returns an AnnData containing the retained raw-count features.

Version **0.13.0** implements continuous-v0.3 evidence fusion. It calculates replication adequacy from the discovery cohort and continuously weights both evidence components. It does not select one branch at a routing threshold.

[中文说明](README_CN.md) · [Mathematics](docs/METHODS.md) · [Reproducibility](docs/MIGRATION.md) · [Tests](docs/TEST_REPORT.md)

## Install

Use Python 3.10 or later. Install into the same environment as your notebook kernel, then restart that kernel.

```bash
git clone https://github.com/babySherlock/HVGdecision.git
cd HVGdecision/HVGDecision_v0.13.0
python -m pip install ./dist/hvgdecision-0.13.0-py3-none-any.whl
```

Alternatively, run `python -m pip install .` in this package directory. The import remains `import hvgdecision as hd`. This updates the same-named package in the selected environment; use a separate environment to retain an older installation.

## Reference / Query analysis

Read and filter your data first. These donor lists are the Kang2018 control split used in the supplied example, not general defaults.

```python
from datetime import datetime
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/filtered_counts.h5ad')
counts = hd.find_raw_counts(adata)
if not counts.valid:
    raise ValueError(counts.audit.to_string(index=False))

reference = [
    'patient_101', 'patient_1015', 'patient_1016',
    'patient_1039', 'patient_107', 'patient_1244',
]
query = ['patient_1256', 'patient_1488']

result = hd.refine(
    adata,
    batch_key='replicate', label_key='cell_type', counts=counts,
    reference=reference, query=query,
    n_hvg=2000, hvg_span=0.5, query_hvg_span=0.3,
    output_dir='./HVGDecision_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'),
    return_details=True,
)
print(result)
print('Removed:', result.removed_genes)
print(result.routing.audit.to_string(index=False))
adata_hvg = result.adata
```

Reference and Query must be disjoint and cover all input cells. Reference labels determine risk evidence, protection and replication adequacy. Query HVGs are fitted independently using Query counts and batch labels, without Query cell-type labels. The final panel is the Query panel minus its intersection with Reference risk calls. All input cells remain in the output.

The fixed spans above match the previously audited Kang feature-selection settings. The general default `hvg_span='auto'` tries 0.3, 0.5, 0.7 and 1.0 only after numerical fitting failure; each cohort starts independently. Attempts are saved. A numeric span fails rather than silently changing.

## Pooled discovery

For a labeled, same-tissue multi-donor cohort, omit both donor lists:

```python
result = hd.refine(
    adata, batch_key='donor_id', label_key='cell_type',
    counts=hd.find_raw_counts(adata), n_hvg=2000,
    return_details=True,
)
```

All input cells then contribute to discovery. This is pooled, label-informed analysis, not independent Query validation. If only `reference` is supplied, discovery uses that subset and the Reference panel is applied to all input cells.

## Counts and external panels

Counts must be finite, nonnegative and integer-like. Detection audits layers, X and raw.X; the name raw does not establish that a matrix contains counts. Numerical checks cannot establish data provenance.

Specify a layer by name with `hd.find_raw_counts(adata, source='counts')`. Cells-by-genes CSV files or DataFrames with cell and gene identifiers are also supported; see [counts usage](docs/COUNTS.md). Missing or duplicated identifiers are errors.

For frozen or externally ranked panels:

```python
import pandas as pd
reference_genes = pd.read_csv('reference_hvg_panel.csv')['gene'].astype(str).tolist()
query_genes = pd.read_csv('query_hvg_panel.csv')['gene'].astype(str).tolist()
result = hd.refine(
    adata, batch_key='replicate', label_key='cell_type', counts=counts,
    reference=reference, query=query,
    hvg_genes=reference_genes, query_hvg_genes=query_genes,
    return_details=True,
)
```

External panels bypass HVG fitting. Genes absent from counts cause an error; they are not renamed or silently dropped. `n_hvg` sets the fitted base-panel size, not a deletion target or a budget search. `protected_genes=[...]` prevents deletion of user-specified genes.

## Outputs

With `return_details=True`, inspect `adata`, `removed_genes`, `decision_table`, `reference_decision_table`, `component_evidence`, `routing` and `run_info`. The name `routing` is retained for API compatibility; its current mode is `continuous`.

A new output folder contains:

- `adata_hvg.h5ad`: all input cells and retained genes; raw counts in X and layers['counts'].
- `gene_decisions.csv`, `removed_genes.csv`, `retained_genes.csv`: target-panel decisions.
- `reference_gene_decisions.csv`, `reference_risk_genes.csv`: discovery-panel decisions.
- `donor_evidence.csv`, `cell_evidence.csv`: component evidence and protection flags.
- `design_audit.csv`, `donor_celltype_coverage.csv`, `discovery_cell_ids.csv`: cohort/design audit.
- Reference/Query HVG panels, fitting audits and span attempts.
- `run_manifest.json`: parameters, seeds, panel/source hashes and dependency versions.

Existing nonempty output folders are never overwritten. Zero deletions are allowed. Query genes outside the Reference test panel are marked untested, not safe by inference.

## Downstream use and limitations

`refine` performs feature refinement only. It does not fit integration models, select seeds or optimize ARI. Use the counts layer for scVI; apply the appropriate normalization and representation steps for Harmony, BBKNN, Scanorama or Seurat CCA. Old embeddings and graphs are not carried into the output.

The continuous score is not a calibrated probability or an FDR. Cell-level permutation tests are conditional diagnostics, not independent-donor population inference. Confounding is audited separately from A. Set `design_policy='error'` to stop on a non-identifiable supported design; the default `'report'` warns and records the limitation. At least two observed discovery donors and two cell types are required.

Historical PBMC and Lung deletion counts must be checked against frozen inputs and gene-level evidence, not imposed as targets. Release tests do not establish improved performance on every dataset.

## Examples and tests

[PBMC notebook](examples/01_PBMC_reference_query.ipynb) · [Human Lung notebook](examples/02_HumanLung_pooled.ipynb)

```bash
python -m pip install '.[dev]'
python -m pytest -q
```

MIT license. Research use; no clinical interpretation is implied.- [Weighted A and continuous fusion](HVGDecision_v0.13.0/docs/METHODS.md)
- [Counts sources](HVGDecision_v0.13.0/docs/COUNTS.md)
- [Migration and historical replay](HVGDecision_v0.13.0/docs/MIGRATION.md)
- [Release verification](HVGDecision_v0.13.0/docs/TEST_REPORT.md)

Twelve local tests passed for source and installed-wheel execution. These checks are not a rerun of the full server benchmark and do not guarantee any fixed deletion count. The continuous score is not a calibrated probability or final-set FDR. Cell-level permutation tests do not replace independent donor replication.

Old software versions remain recoverable from Git history. MIT license.
