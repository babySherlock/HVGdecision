from pathlib import Path
from datetime import datetime
import scanpy as sc
import hvgdecision as hd

INPUT_H5AD = Path("/data/gj/project/20260615/HVGDecision_data/human_lung_travaglini_10x.h5ad")
BATCH_KEY, LABEL_KEY = "donor_id", "cell_type"
OUTPUT_DIR = Path("/data/gj/project/20260615/HVGDecision_results") / ("HumanLung_pooled_continuous_013_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
adata = sc.read_h5ad(INPUT_H5AD)
if "tissue_general" in adata.obs:
    mask = adata.obs["tissue_general"].astype(str).str.lower().eq("lung")
elif "tissue" in adata.obs:
    mask = adata.obs["tissue"].astype(str).str.lower().str.contains("lung", regex=False)
else:
    raise KeyError("Expected tissue_general or tissue for explicit lung-only filtering")
if not mask.any():
    raise ValueError("No lung tissue cells matched")
adata = adata[mask].copy()
assert adata.n_obs == 60993, f"Historical cohort mismatch: {adata.n_obs} cells"
assert adata.obs[BATCH_KEY].nunique() == 3
assert adata.obs[LABEL_KEY].nunique() == 46
print("Version:", hd.__version__, "Cells:", adata.n_obs)
counts = hd.find_raw_counts(adata)
print(counts.audit.to_string(index=False))
if not counts.valid:
    raise ValueError(counts.error)

design = hd.audit_design(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY)
print(design.audit.to_string(index=False))
# Pooled: all retained lung cells contribute to discovery.
# For historical panel replay, pass hvg_genes=<frozen ordered gene list>.
result = hd.refine(
    adata, batch_key=BATCH_KEY, label_key=LABEL_KEY, counts=counts,
    reference=None, query=None, n_hvg=2000, hvg_span="auto",
    seed=20260829, output_dir=OUTPUT_DIR, return_details=True,
)
print(result)
print("Actually removed:", result.removed_genes)
print("Output:", OUTPUT_DIR)
adata_hvg = result.adata

