
# HVGDecision v0.10.0 Reference→Query patch (core-preserving)

This companion patch does **not** reimplement the HVGDecision v0.10.0 donor-aware risk engine. Reference-side risk discovery is delegated to the installed `hvgdecision.refine()`.

It adds the missing Reference→Query workflow semantics:

- every independent Seurat-v3 fit starts at `span=0.3` and retries `0.5`, `0.7`, `1.0` only after numerical failure;
- Query fitting is a new independent fit and therefore starts at `0.3` again;
- the Reference risk list is transferred to Query by intersection with Query-specific HVGs;
- final panel: `H_final = H_query \\ (R_reference ∩ H_query)`.

Install alongside HVGDecision v0.10.0:

```bash
python -m pip install ./hvgdecision_reference_query_patch-0.10.1-py3-none-any.whl --no-deps
```

Usage:

```python
import hvgdecision_reference_query_patch as hrq

result = hrq.refine_reference_query(
    adata,
    batch_key="replicate",
    label_key="cell_type",
    counts=raw_counts,
    reference=["patient_101", "patient_1015", "patient_1016", "patient_1039", "patient_107", "patient_1244"],
    query=["patient_1256", "patient_1488"],
    n_hvg=2000,
    output_dir="./HVGDecision_Kang_reference6_to_query2",
)

print(result.reference_risk_genes)
print(result.removed_genes)
print(result.final_n_hvg)
adata_hvg = result.adata
```

For the reproduced Kang2018 control split, Reference used span 0.5 after 0.3 failed; Query independently succeeded at span 0.3; 19 Reference risk genes yielded three Query deletions: `CSF2`, `GJB2`, `AC147651.3`, giving 1997 final genes. These names are not hard-coded.
