"""Cross-domain example starting from two independent h5ad files."""
import scanpy as sc
import hvgdecision as hd

ref = sc.read_h5ad("reference_10x.h5ad")
qry = sc.read_h5ad("query_smartseq2.h5ad")

combined, audit = hd.prepare_cross_domain_inputs(
    ref,
    qry,
    reference_counts="auto",
    query_counts="auto",
    domain_key="hvgdecision_domain",
    return_audit=True,
)
print(audit.to_string(index=False))

study = hd.setup_reference_query(
    combined,
    mode="cross_domain",
    batch_key="donor",
    label_key="cell_type",
    split_key="hvgdecision_domain",
    reference=["reference"],
    query=["query"],
    counts_layer="auto",
    output_dir="HVGDecision_results/cross_domain",
)

result = study.run(
    n_hvg=2000,
    cross_domain_delete_budget=5,
    return_details=True,
)
print("removed:", result.harmful_genes)
print("final N:", result.final_n_hvg)
