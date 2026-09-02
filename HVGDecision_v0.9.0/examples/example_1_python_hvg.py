"""Cross-domain example using internally selected Query HVGs."""
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("combined_reference_query.h5ad")

study = hd.setup_reference_query(
    adata,
    mode="cross_domain",
    batch_key="donor",
    label_key="cell_type",
    split_key="technology",
    reference=["10X"],
    query=["smartseq2"],
    counts_layer="counts",
    output_dir="HVGDecision_results/cross_domain",
)

result = study.run(
    n_hvg=2000,
    cross_domain_delete_budget=5,
    return_details=True,
)

print(result.harmful_genes)
print(result.final_n_hvg)
