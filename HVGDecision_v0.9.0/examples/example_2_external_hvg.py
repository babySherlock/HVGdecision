"""Refine an external frozen HVG list in the selected mode."""
import pandas as pd
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")

study = hd.setup_reference_query(
    adata,
    mode="within_domain",  # change to cross_domain when appropriate
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="counts",
    output_dir="HVGDecision_results/external_hvg",
)

hvg = pd.read_csv("hvg_table.csv")
result = study.refine_hvg(
    hvg,
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)

print(result.harmful_genes)
print(result.final_hvg_genes[:20])
