# HVGDecision

HVGDecision refines HVG panels for integration of multi-donor single-cell RNA-seq datasets. The current distribution is **v0.14.0**, implementing **dual-evidence signed-margin fusion**.

For every candidate gene, donor-aware and cell-level conditional donor-leakage evidence are converted to standardized signed margins, subjected to component-specific biological protection, and combined by centered log-mean-exp fusion:

\[
M_g^{(u)}=\tau_f\log\left[\frac{\exp(\widetilde M_g^{(d)}/\tau_f)+\exp(\widetilde M_g^{(c)}/\tau_f)}{2}\right],
\qquad \tau_f=0.10.
\]

A discovery-risk gene satisfies `M_u > 0`. Replication adequacy and design identifiability are reported as audits and do not enter the gene-level fusion calculation.

## Install

```bash
git clone https://github.com/babySherlock/HVGdecision.git
cd HVGdecision/HVGDecision_v0.14.0
python -m pip install ./dist/hvgdecision-0.14.0-py3-none-any.whl
```

or:

```bash
python -m pip install .
```

Then:

```python
import hvgdecision as hd
print(hd.__version__)
```

## Minimal use

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/data.h5ad')
counts = hd.find_raw_counts(adata)

result = hd.refine(
    adata,
    batch_key='donor',
    label_key='cell_type',
    counts=counts,
    n_hvg=2000,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    return_details=True,
)

print(result.removed_genes)
print(result.design.audit.to_string(index=False))
```

## Documentation

- [Complete usage](HVGDecision_v0.14.0/README.md)
- [中文说明](HVGDecision_v0.14.0/README_CN.md)
- [Mathematics](HVGDecision_v0.14.0/docs/METHODS.md)
- [v0.13 → v0.14 migration](HVGDecision_v0.14.0/docs/MIGRATION.md)
- [Counts](HVGDecision_v0.14.0/docs/COUNTS.md)
- [Release tests](HVGDecision_v0.14.0/docs/TEST_REPORT.md)
- [PBMC Reference/Query example](HVGDecision_v0.14.0/examples/01_PBMC_reference_query.ipynb)
- [Human Lung pooled example](HVGDecision_v0.14.0/examples/02_HumanLung_pooled.ipynb)

MIT License. Research software; no clinical interpretation is implied.
