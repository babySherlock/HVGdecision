# HVGDecision 0.9.1

Integration-aware HVG refinement before single-cell data integration.

## Install from the wheel

Normal users should install the bundled wheel rather than use an editable source install:

```bash
python -m pip install ./hvgdecision-0.9.1-py3-none-any.whl
```

or, if using the copy under `dist/`:

```bash
python -m pip install ./dist/hvgdecision-0.9.1-py3-none-any.whl
```

Verify:

```bash
python -c "import hvgdecision as hd; print(hd.__version__); print(hd.VALID_MODES)"
```

`pip install -e .` is only for editable/development installs.

## Explicit mode selection

Every analysis must explicitly choose one mode:

```text
within_domain
cross_domain
```

- `within_domain`: shared-domain multi-donor/multi-batch refinement.
- `cross_domain`: Reference→Query refinement across datasets or technologies.

## Raw counts are detected automatically

The default is:

```python
counts_layer="auto"
```

Inspect the selected source with:

```python
import hvgdecision as hd
check = hd.find_raw_counts(adata)
print(check.location)
print(check.audit)
```

HVGDecision audits common count layers, `adata.raw.X`, `adata.X`, and unusually named layers. An explicit source should only be supplied when automatic detection fails or when the user intentionally wants a specific source.

## Within-domain quick start

```python
study = hd.setup_reference_query(
    adata,
    mode="within_domain",
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2", "D3"],
    query=["D4"],
    counts_layer="auto",
)

result = study.run(return_details=True)
```

Within-domain risk:

$$
R_g^{within}=Z(L_g)+0.75Z(I_g)-Z(B_g)
$$

## Cross-domain: two separate AnnData objects

Reference and Query may keep raw counts in different AnnData locations. Prepare them safely with:

```python
combined, input_audit = hd.prepare_cross_domain_inputs(
    reference_adata,
    query_adata,
    reference_counts="auto",
    query_counts="auto",
    domain_key="hvgdecision_domain",
    return_audit=True,
)

study = hd.setup_reference_query(
    combined,
    mode="cross_domain",
    batch_key="donor",
    label_key="cell_type",
    split_key="hvgdecision_domain",
    reference=["reference"],
    query=["query"],
    counts_layer="auto",
)

result = study.run(
    n_hvg=2000,
    cross_domain_delete_budget=5,
    return_details=True,
)
```

The helper independently audits Reference and Query raw counts, intersects raw-count gene identifiers, and builds a combined AnnData with `layers['counts']`.

Dataset/domain shift:

$$
S_g=0.70E_g^*+0.30D_g^*
$$

Reference biology protection:

$$
B_g=0.60B_g^{celltype}+0.25M_g+0.15P_g
$$

Final Cross-domain Rule V3:

$$
R_g^{cross}
=
\min(R_g,S_g)
\left[
0.75+0.25\left(1-\left|R_g-S_g\right|\right)
\right]
\left(1-0.75B_g\right)
$$

Query true biological labels are never used in cross-domain feature selection.

## Use a local HVG file

A local CSV/TSV can contain a `gene` column. A `.txt` file may contain one gene per line with no header.

```python
result = study.refine_hvg(
    "my_hvg.csv",  # or my_hvg.tsv / my_hvg.txt
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

In `within_domain`, the imported file is the frozen within-domain base panel. In `cross_domain`, it is the frozen Query base HVG panel.

For the complete Chinese tutorial, see `README_CN.md` and `docs/CROSS_DOMAIN_TUTORIAL_CN.md`.
