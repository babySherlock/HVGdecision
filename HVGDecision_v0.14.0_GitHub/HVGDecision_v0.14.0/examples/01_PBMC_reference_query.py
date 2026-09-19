"""Kang2018 PBMC Reference -> Query example for HVGDecision v0.14.0."""

from datetime import datetime
from pathlib import Path

import pandas as pd
import scanpy as sc
import hvgdecision as hd


INPUT_H5AD = Path("/path/to/kang_2018.h5ad")
OUTPUT_ROOT = Path("./results")

REFERENCE_DONORS = [
    "patient_101",
    "patient_1015",
    "patient_1016",
    "patient_1039",
    "patient_107",
    "patient_1244",
]
QUERY_DONORS = ["patient_1256", "patient_1488"]

BATCH_KEY = "replicate"
LABEL_KEY = "cell_type"

# Optional frozen panels. Set to None to refit Seurat v3 HVGs.
REFERENCE_HVG_CSV = None  # e.g. Path("reference_hvg_2000.csv")
QUERY_HVG_CSV = None      # e.g. Path("query_hvg_2000.csv")

OUTPUT_DIR = OUTPUT_ROOT / (
    "PBMC_HVGDecision_v014_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")
)

adata = sc.read_h5ad(INPUT_H5AD)

# Kang2018 example: retain unstimulated/control cells explicitly.
if "label" in adata.obs:
    adata = adata[
        adata.obs["label"].astype(str).str.lower().isin(["ctrl", "control"])
    ].copy()

adata = adata[
    adata.obs[BATCH_KEY].astype(str).isin(REFERENCE_DONORS + QUERY_DONORS)
].copy()

counts = hd.find_raw_counts(adata)
print(counts.audit.to_string(index=False))
if not counts.valid:
    raise ValueError(counts.error)

reference_genes = None
query_genes = None
if REFERENCE_HVG_CSV is not None:
    reference_genes = pd.read_csv(REFERENCE_HVG_CSV)["gene"].astype(str).tolist()
if QUERY_HVG_CSV is not None:
    query_genes = pd.read_csv(QUERY_HVG_CSV)["gene"].astype(str).tolist()

# Dataset-level audit only; it does not enter gene-level fusion.
design = hd.audit_design(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    reference=REFERENCE_DONORS,
)
print(design.audit.to_string(index=False))

result = hd.refine(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    counts=counts,
    reference=REFERENCE_DONORS,
    query=QUERY_DONORS,
    n_hvg=2000,
    hvg_genes=reference_genes,
    query_hvg_genes=query_genes,
    hvg_span=0.5,
    query_hvg_span=0.3,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    seed=20260829,
    output_dir=OUTPUT_DIR,
    return_details=True,
)

print(result)
print("Reference discovery-risk genes:", result.reference_risk_genes)
print("Final Query removed genes:", result.removed_genes)
print("Output:", OUTPUT_DIR)

adata_hvg = result.adata
