# HVGDecision

HVGDecision refines highly variable gene (HVG) panels before integration of multi-donor single-cell RNA-seq datasets. It jointly evaluates donor-aware and cell-level conditional donor-leakage evidence, converts both components to standardized signed margins, applies component-specific biological protection, and combines them with centered log-mean-exp fusion.

Version **0.14.0** implements the manuscript v0.4 decision rule:

\[
M_g^{(u)} = \tau_f\log\left[\frac{\exp(\widetilde M_g^{(d)}/\tau_f)+\exp(\widetilde M_g^{(c)}/\tau_f)}{2}\right],
\qquad \tau_f=0.10,
\]

with discovery-risk call

\[
M_g^{(u)}>0.
\]

Replication adequacy and design identifiability are reported as dataset-level audits and **do not enter the gene-level fusion calculation**.

[中文说明](README_CN.md) · [Mathematics](docs/METHODS.md) · [Migration from v0.13](docs/MIGRATION.md) · [Counts](docs/COUNTS.md) · [Test report](docs/TEST_REPORT.md)

## Install

Use Python 3.10 or later. Install into the same environment as your notebook kernel, then restart the kernel.

From the supplied wheel:

```bash
python -m pip install ./dist/hvgdecision-0.14.0-py3-none-any.whl
```

or from source:

```bash
python -m pip install .
```

The public import remains:

```python
import hvgdecision as hd
print(hd.__version__)
```

Expected version:

```text
0.14.0
```

## Reference / Query analysis

Read and QC-filter the data first. The donor lists below are an example, not package defaults.

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
    batch_key='replicate',
    label_key='cell_type',
    counts=counts,
    reference=reference,
    query=query,
    n_hvg=2000,
    hvg_span=0.5,
    query_hvg_span=0.3,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    output_dir='./HVGDecision_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'),
    return_details=True,
)

print(result)
print('Reference discovery-risk genes:', result.reference_risk_genes)
print('Final removed genes:', result.removed_genes)
print(result.design.audit.to_string(index=False))

adata_hvg = result.adata
```

Reference and Query must be disjoint and cover all input cells. Risk evidence is learned from the Reference cohort. Query HVGs are fitted independently using Query expression and batch labels; Query cell-type labels are not used for feature refinement. The final removed set is

\[
R_{\mathrm{final}}
=
R_{\mathrm{Reference\ discovery}}
\cap
HVG_{\mathrm{Query}}.
\]

Therefore, the number of Reference discovery-risk genes can be larger than the number of genes finally removed from the Query panel.

## Pooled discovery

For a labeled same-tissue multi-donor cohort, omit `reference` and `query`:

```python
result = hd.refine(
    adata,
    batch_key='donor_id',
    label_key='cell_type',
    counts=hd.find_raw_counts(adata),
    n_hvg=2000,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    return_details=True,
)
```

All input cells contribute to discovery. This is pooled label-informed refinement, not independent external validation.

## Frozen / external HVG panels

You can bypass Seurat v3 HVG fitting by supplying an ordered list of gene IDs:

```python
import pandas as pd

reference_genes = (
    pd.read_csv('reference_hvg_panel.csv')['gene']
    .astype(str)
    .tolist()
)

query_genes = (
    pd.read_csv('query_hvg_panel.csv')['gene']
    .astype(str)
    .tolist()
)

result = hd.refine(
    adata,
    batch_key='replicate',
    label_key='cell_type',
    counts=counts,
    reference=reference,
    query=query,
    hvg_genes=reference_genes,
    query_hvg_genes=query_genes,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    return_details=True,
)
```

Genes absent from the count matrix are errors; they are not renamed or silently dropped.

## The two evidence components

### Cell-level conditional donor-leakage evidence

This component evaluates donor-associated variation from individual-cell expression after conditioning on cell type. It includes conditional donor leakage, within-cell-type donor instability, cell-type biological support, permutation FDR and bootstrap stability.

### Donor-aware evidence

This component evaluates technical non-replication across donors, including within-cell-type donor nuisance, single-donor dominance, cross-donor heterogeneity, leave-one-donor-out instability, directional agreement and biological support.

Both components are evaluated for every candidate gene.

## Signed-margin fusion

Each risk criterion is converted into a standardized signed margin. Positive values satisfy the corresponding risk condition. AND conditions are represented by the minimum margin and OR conditions by the maximum margin.

After component-specific protection:

\[
\widetilde M_g^{(b)}=
\begin{cases}
M_g^{(b)}, & \pi_g^{(b)}=0,\\
-\infty, & \pi_g^{(b)}=1,
\end{cases}
\qquad b\in\{c,d\}.
\]

The two effective margins are fused by centered log-mean-exp:

\[
M_g^{(u)} = \tau_f\log\left[\frac{\exp(\widetilde M_g^{(d)}/\tau_f)+\exp(\widetilde M_g^{(c)}/\tau_f)}{2}\right].
\]

Properties:

- symmetric in the two evidence components;
- if both margins equal `m`, the unified margin is exactly `m`;
- when one component dominates, `M_u ≈ max(M_d, M_c) - tau_f*log(2)`;
- smaller `tau_f` approaches the larger effective margin;
- larger `tau_f` imposes a stronger disagreement penalty;
- no dataset-dependent donor-versus-cell fusion weight is used.

## Fusion temperature

Default:

```python
cfg = hd.UnifiedFusionConfig(tau=0.10)
```

The default is a method-development setting and should be evaluated on independent datasets rather than tuned separately for each downstream integration algorithm.

For direct mathematical sensitivity analysis:

```python
import numpy as np
import hvgdecision as hd

mu = hd.centered_logmeanexp(
    donor_margin=np.array([0.20, 0.05]),
    cell_margin=np.array([0.20, -0.40]),
    tau=0.10,
)
print(mu)
```

## Counts

Counts must be finite, nonnegative and integer-like. Detection audits `layers`, `X` and `raw.X`; the name `raw` alone does not prove raw-count provenance.

Specify a layer explicitly if needed:

```python
counts = hd.find_raw_counts(adata, source='counts')
```

See [docs/COUNTS.md](docs/COUNTS.md).

## Outputs

With `return_details=True`, the returned `RefinementResult` provides:

- `adata`: all input cells and retained HVGs, with raw counts in `X` and `layers['counts']`;
- `removed_genes`: genes removed from the final target panel;
- `reference_risk_genes`: Reference/pooled discovery-risk genes;
- `decision_table`: decisions on the final target panel;
- `reference_decision_table`: complete discovery-panel gene-level evidence and calls;
- `component_evidence['donor']`: donor-aware evidence table;
- `component_evidence['cell']`: cell-level evidence table;
- `design`: replication-adequacy and design-identifiability audit;
- `run_info`: parameters, hashes, dependency versions and fusion manifest.

A saved output directory contains:

- `adata_hvg.h5ad`
- `gene_decisions.csv`
- `reference_gene_decisions.csv`
- `reference_risk_genes.csv`
- `removed_genes.csv`
- `retained_genes.csv`
- `donor_evidence.csv`
- `cell_evidence.csv`
- `design_audit.csv`
- `donor_celltype_coverage.csv`
- `design_by_donor.csv`
- `design_by_celltype.csv`
- Reference / Query HVG panel files and HVG fitting audits
- `run_manifest.json`

Existing non-empty output folders are never overwritten. Zero deletions are permitted.

## Design audit

```python
report = hd.audit_design(
    adata,
    batch_key='donor',
    label_key='cell_type',
)

print(report.audit)
```

Replication adequacy and design identifiability are metadata audits. They do not select an evidence component and do not enter the centered log-mean-exp fusion equation.

## Command-line interface

Audit only:

```bash
hvgdecision audit \
  --input data.h5ad \
  --batch-key donor \
  --label-key cell_type \
  --output audit_out
```

Refinement:

```bash
hvgdecision refine \
  --input data.h5ad \
  --batch-key donor \
  --label-key cell_type \
  --n-hvg 2000 \
  --tau-f 0.10 \
  --output refinement_out
```

## Examples

- `examples/01_PBMC_reference_query.py`
- `examples/02_HumanLung_pooled.py`
- corresponding Jupyter notebooks

## Tests

```bash
python -m pip install '.[dev]'
python -m pytest -q
```

The release tests cover the fusion mathematics, design audit, Reference–Query transport, protection, raw-count handling, repeatability and HVG fitting behavior.

## Version migration

Version 0.13 used replication-adequacy-dependent confidence weighting. Version 0.14 replaces that final fusion layer with centered log-mean-exp fusion of effective signed margins. The donor-aware and cell-level evidence engines and signed-gate reconstruction remain the basis of the two component margins.

See [docs/MIGRATION.md](docs/MIGRATION.md).

## Citation

A manuscript describing HVGDecision is in preparation. Replace this section with the final publication citation after acceptance.

## License

MIT license. Research software; no clinical interpretation is implied.
