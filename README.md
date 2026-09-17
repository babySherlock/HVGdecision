# HVGDecision

HVGDecision refines HVG panels for multi-donor single-cell integration. It combines donor-level and cell-level risk evidence with biological protection and returns retained raw-count features together with gene-level decisions.

The current distribution is **v0.13.0**, implementing **continuous-v0.3 evidence fusion**. Replication adequacy A is calculated from donor-by-cell-type coverage using a weighted harmonic mean. Both evidence components are computed and continuously weighted; no 0.75 hard routing threshold is used.

## Install

```bash
git clone https://github.com/babySherlock/HVGdecision.git
cd HVGdecision/HVGDecision_v0.13.0
python -m pip install ./dist/hvgdecision-0.13.0-py3-none-any.whl
```

Install in the same environment as your notebook and restart its kernel. The import is `import hvgdecision as hd`. Source installation with `python -m pip install .` is also supported.

## Reference / Query example

The donor IDs below are the Kang2018 control example. Filter the input to the intended condition and donors first.

```python
from datetime import datetime
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/filtered_counts.h5ad')
raw_counts = hd.find_raw_counts(adata)
if not raw_counts.valid:
    raise ValueError(raw_counts.audit.to_string(index=False))

result = hd.refine(
    adata,
    batch_key='replicate', label_key='cell_type', counts=raw_counts,
    reference=['patient_101', 'patient_1015', 'patient_1016',
               'patient_1039', 'patient_107', 'patient_1244'],
    query=['patient_1256', 'patient_1488'],
    n_hvg=2000, hvg_span=0.5, query_hvg_span=0.3,
    output_dir='./HVGDecision_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'),
    return_details=True,
)
print(result.removed_genes)
print(result.routing.audit.to_string(index=False))
adata_hvg = result.adata
```

Risk evidence and A use Reference only. Query HVGs are fitted independently, without Query cell-type labels, and only genes overlapping the Reference risk list are removed. Output retains all input cells. Omit both `reference` and `query` for pooled discovery; that is label-informed analysis of the full cohort, not independent Query validation.

`n_hvg=2000` sets the Scanpy Seurat v3 base-panel size. It is not an HVG-budget search or a deletion target. External ordered panels and user-protected genes are supported. The output contains raw counts; downstream integration still requires the appropriate preprocessing.

## Documentation

- [Complete usage and outputs](HVGDecision_v0.13.0/README.md)
- [中文说明](HVGDecision_v0.13.0/README_CN.md)
- [PBMC Reference/Query notebook](HVGDecision_v0.13.0/examples/01_PBMC_reference_query.ipynb)
- [Human Lung pooled notebook](HVGDecision_v0.13.0/examples/02_HumanLung_pooled.ipynb)
- [Weighted A and continuous fusion](HVGDecision_v0.13.0/docs/METHODS.md)
- [Counts sources](HVGDecision_v0.13.0/docs/COUNTS.md)
- [Migration and historical replay](HVGDecision_v0.13.0/docs/MIGRATION.md)
- [Release verification](HVGDecision_v0.13.0/docs/TEST_REPORT.md)

Twelve local tests passed for source and installed-wheel execution. These checks are not a rerun of the full server benchmark and do not guarantee any fixed deletion count. The continuous score is not a calibrated probability or final-set FDR. Cell-level permutation tests do not replace independent donor replication.

Old software versions remain recoverable from Git history. MIT license.
