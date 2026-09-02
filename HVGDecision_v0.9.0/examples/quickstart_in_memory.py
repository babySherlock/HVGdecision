import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")

# REQUIRED: choose the algorithm mode explicitly.
study = hd.setup_reference_query(
    adata,
    mode="within_domain",
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2"],
    query=["D3"],
    counts_layer="counts",
    output_dir="HVGDecision_results/quickstart",
)

result = study.run(return_details=True)
print("mode:", study.mode)
print("harmful genes:", result.harmful_genes)
print("final N:", result.final_n_hvg)
