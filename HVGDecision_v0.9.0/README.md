# HVGDecision 0.9.0

Integration-aware HVG refinement before single-cell data integration.

**Mode selection is explicit in 0.9.0:**

- `within_domain`: shared-domain multi-donor/multi-batch Reference-only refinement.
- `cross_domain`: cross-dataset/cross-technology Reference→Query Rule V3.

```python
import hvgdecision as hd

study = hd.setup_reference_query(
    adata,
    mode="within_domain",  # or "cross_domain"
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="counts",
    output_dir="HVGDecision_results/example",
)

result = study.run(return_details=True)
```

Cross-domain example:

```python
study = hd.setup_reference_query(
    adata,
    mode="cross_domain",
    batch_key="donor",
    label_key="cell_type",
    split_key="technology",
    reference=["10X"],
    query=["Smart-seq2"],
    counts_layer="counts",
)

result = study.run(
    n_hvg=2000,
    cross_domain_delete_budget=5,
    return_details=True,
)
```

## Final within-domain rule

\[
R_g=Z(L_g)+0.75Z(I_g)-Z(B_g)
\]

A gene is removed only when BH-FDR, leakage/risk effect floors, bootstrap recurrence and biological-protection gates all pass.

## Final cross-domain Rule V3

\[
S_g=0.70E_g^*+0.30D_g^*
\]

\[
CrossRisk_g=\min(R_g,S_g)[0.75+0.25(1-|R_g-S_g|)](1-0.75B_g)
\]

Query true labels are never used by cross-domain feature selection.

## Optional benchmark aggregation

`hvgdecision.three_domain_scores()` computes a scIB-inspired three-domain summary (Transfer / Batch correction / Biological fidelity) after panel-wise min-max normalization within each dataset × integration method task.

See `README_CN.md` and `docs/METHODS.md` for the full method definition.
