import hvgdecision as hd
import hvgdecision_reference_query_patch as hrq

# adata and raw_counts should already be prepared.
reference_donors = [
    "patient_101",
    "patient_1015",
    "patient_1016",
    "patient_1039",
    "patient_107",
    "patient_1244",
]
query_donors = ["patient_1256", "patient_1488"]

result = hrq.refine_reference_query(
    adata,
    batch_key="replicate",
    label_key="cell_type",
    counts=raw_counts,
    reference=reference_donors,
    query=query_donors,
    n_hvg=2000,
    output_dir="./HVGDecision_Kang_reference6_to_query2",
)

print(result)
print("Reference risk:", result.reference_risk_genes)
print("Final Query deletion:", result.removed_genes)
print("Final n:", result.final_n_hvg)
adata_hvg = result.adata
