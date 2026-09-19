"""Travaglini Human Lung pooled example for HVGDecision v0.14.0."""

from datetime import datetime
from pathlib import Path

import pandas as pd
import scanpy as sc
import hvgdecision as hd


INPUT_H5AD = Path("/path/to/human_lung_travaglini_10x.h5ad")
OUTPUT_ROOT = Path("./results")

BATCH_KEY = "donor_id"
LABEL_KEY = "cell_type"

# Optional frozen pooled HVG panel. Set to None to refit Seurat v3 HVGs.
POOLED_HVG_CSV = None  # e.g. Path("pooled_hvg_2000.csv")

OUTPUT_DIR = OUTPUT_ROOT / (
    "HumanLung_HVGDecision_v014_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")
)

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

counts = hd.find_raw_counts(adata)
print(counts.audit.to_string(index=False))
if not counts.valid:
    raise ValueError(counts.error)

pooled_genes = None
if POOLED_HVG_CSV is not None:
    pooled_genes = pd.read_csv(POOLED_HVG_CSV)["gene"].astype(str).tolist()

design = hd.audit_design(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
)
print(design.audit.to_string(index=False))

result = hd.refine(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    counts=counts,
    reference=None,
    query=None,
    n_hvg=2000,
    hvg_genes=pooled_genes,
    hvg_span="auto",
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    seed=20260829,
    output_dir=OUTPUT_DIR,
    return_details=True,
)

print(result)
print("Pooled discovery-risk / final removed genes:", result.removed_genes)
print("Output:", OUTPUT_DIR)

adata_hvg = result.adata
