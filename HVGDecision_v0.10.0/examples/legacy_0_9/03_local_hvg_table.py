"""Use a local HVG CSV/TSV/TXT without re-running HVG selection."""
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")
study = hd.setup_reference_query(
    adata,
    mode="within_domain",  # or cross_domain
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="auto",
    output_dir="HVGDecision_results/local_hvg",
)

result = study.refine_hvg(
    "my_hvg.txt",  # CSV/TSV also supported
    method_name="external",
    initial_n_hvg=2000,
    return_details=True,
)
print(result.final_hvg_genes[:20])
