# HVGDecision

HVGDecision refines highly variable gene (HVG) panels in multi-donor single-cell
datasets from comparable biological settings. It balances donor-associated risk
against cell-type support, reports each gene decision, and returns an AnnData
object containing the retained raw-count features.

[中文用法](README_CN.md) | [Methods and defaults](docs/METHODS.md) |
[Output schema](docs/OUTPUT_SCHEMA.md) | [Migration](docs/MIGRATION.md)

## Installation

Python 3.10 or newer is required. From the downloaded source folder:

```bash
python -m pip install .
```

Or install the matching release wheel:

```bash
python -m pip install hvgdecision-0.10.0-py3-none-any.whl
```

Restart the notebook kernel after upgrading. This release does not require R or
fit Harmony, BBKNN, Scanorama, scVI or Seurat CCA during feature refinement.
Install the downstream method separately when you need it.

## Quickstart

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/data.h5ad')  # QC-filtered cells, full count gene axis
counts = hd.find_raw_counts(adata)
print(counts.audit)
if not counts.valid:
    raise ValueError(counts.error)

result = hd.refine(
    adata,
    batch_key='donor',
    label_key='cell_type',
    counts=counts,
    n_hvg=2000,
    output_dir='./HVGDecision_run01',  # new or empty directory
    return_details=True,
)
print(result)
print('Removed genes:', result.removed_genes)
print(result.routing.audit)
adata_hvg = result.adata
```

The one-call form returns AnnData directly:

```python
adata_hvg = hd.refine(adata, batch_key='donor', label_key='cell_type')
```

`adata_hvg.X` and `adata_hvg.layers['counts']` contain raw counts for the retained
genes. All input cells and their observation metadata are retained. Old embeddings,
neighbor graphs and normalized expression layers are not copied. Recompute the
preprocessing required by your downstream method. For scVI, use the counts layer:

```python
import scvi
scvi.model.SCVI.setup_anndata(adata_hvg, layer='counts', batch_key='donor')
model = scvi.model.SCVI(adata_hvg)
model.train()
```

## Automatic workflow

1. Check the count source and cell/gene alignment.
2. Audit donor-by-cell-type support and design identifiability.
3. Select a fixed Seurat v3 batch-aware base panel, or use a supplied gene list.
4. Route identifiable designs to `donor_aware` when the routing score is at least
   0.75, otherwise to `cell_level`. Insufficient or confounded designs retain the
   entire base panel and receive a warning.
5. Apply branch-specific risk tests and protection rules, then remove only the
   passing genes. No genes are replenished; removing zero genes is valid.

`audit_design()` runs step 2 without accessing expression values. The routing
score is a design heuristic, not a probability or a learned optimal threshold.

By default all input cells supply selection evidence. `reference=['D1', 'D2']`
restricts HVG selection, routing and risk evidence to those donors, while the
output still contains all input cells. This function performs no held-out
evaluation. If labels were used for feature selection, evaluating label transfer
on those same cells does not constitute an independent transfer test.

## External counts and feature panels

```python
import pandas as pd

# Choose a source explicitly. Invalid sources are never silently replaced.
counts = hd.find_raw_counts(adata, source='raw')  # or 'X', or a layer name
# External CSV/TSV: cell IDs and gene IDs are required; either orientation works.
counts = hd.find_raw_counts(adata, source='/path/to/counts.csv')

# R-derived scran/scry/SCT or any ordered feature list; matching gene IDs required.
panel = pd.read_csv('/path/to/features.csv')['gene'].tolist()
result = hd.refine(
    adata, batch_key='donor', label_key='cell_type', counts=counts,
    hvg_genes=panel, return_details=True,
)
```

External panels bypass Seurat v3. Their length determines the base size; `n_hvg`
does not truncate them. Gene IDs are not silently renamed, dropped or mapped
between symbols and Ensembl IDs. Optional `protected_genes=[...]` protects exact
IDs; the default explicit protection list is empty.

## What the output means

`gene_decisions.csv` contains every base-panel gene, its rank, branch-specific
evidence, protection flags and final decision. `removed_genes.csv` has a header
even when no genes were removed. Full evidence is also stored in
`adata_hvg.uns['hvgdecision']['gene_decisions']`; removed genes cannot be recovered
from the retained-only `.var` table.

The package flags integration-risk candidates. It does not establish that a gene
is biologically harmful or guarantee improved performance for every integration
method. Compare the refined panel with Base2000, direct same-N and rank-matched
random deletion controls using paired seeds. Do not tune deletion thresholds or
select seeds using the final evaluation outcomes.

## Reproducibility and testing

Every run records parameter values, selection donor IDs, panel hashes, engine
hashes and installed dependency versions. [Engine provenance](docs/ENGINE_PROVENANCE.json)
identifies the two notebook implementations. [Testing](docs/TESTING.md) describes
release checks and their limits.

```bash
python -m pip install '.[dev]'
python -m pytest -q
hvgdecision audit --input data.h5ad --batch-key donor --label-key cell_type --output audit_run
hvgdecision refine --input data.h5ad --batch-key donor --label-key cell_type --output refine_run
```

The software does not include patient data. Check dataset access terms and
privacy requirements before sharing `.h5ad` files or donor metadata.
