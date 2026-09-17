from pathlib import Path
from datetime import datetime
import scanpy as sc
import hvgdecision as hd

INPUT_H5AD = Path("/data/gj/project/20260615/R_HVG_六套Benchmark_ReferenceOnly_正式初筛_不含Pancreas_anndataR_v1修复/00_input_cache/kang_2018.h5ad")
REFERENCE_DONORS = ["patient_101", "patient_1015", "patient_1016", "patient_1039", "patient_107", "patient_1244"]
QUERY_DONORS = ["patient_1256", "patient_1488"]
BATCH_KEY, LABEL_KEY = "replicate", "cell_type"
OUTPUT_DIR = Path("/data/gj/project/20260615/HVGDecision_results") / ("PBMC_continuous_013_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))

adata = sc.read_h5ad(INPUT_H5AD)
if "label" not in adata.obs:
    raise KeyError("Expected condition column 'label'; inspect adata.obs before choosing another field.")
adata = adata[adata.obs["label"].astype(str).str.lower().isin(["ctrl", "control"])].copy()
adata = adata[adata.obs[BATCH_KEY].astype(str).isin(REFERENCE_DONORS + QUERY_DONORS)].copy()
assert adata.n_obs == 12315, f"Historical cohort mismatch: {adata.n_obs} cells"
assert adata.obs[BATCH_KEY].nunique() == 8
assert adata.obs[LABEL_KEY].nunique() == 8
print("Version:", hd.__version__, "Cells:", adata.n_obs)
counts = hd.find_raw_counts(adata)
print(counts.audit.to_string(index=False))
if not counts.valid:
    raise ValueError(counts.error)

design = hd.audit_design(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY, reference=REFERENCE_DONORS)
print(design.audit.to_string(index=False))

result = hd.refine(
    adata, batch_key=BATCH_KEY, label_key=LABEL_KEY, counts=counts,
    reference=REFERENCE_DONORS, query=QUERY_DONORS,
    n_hvg=2000, hvg_span=0.5, query_hvg_span=0.3,
    seed=20260829, output_dir=OUTPUT_DIR, return_details=True,
)
print(result)
print("Reference risk genes:", result.reference_risk_genes)
print("Actually removed:", result.removed_genes)
print("Output:", OUTPUT_DIR)
adata_hvg = result.adata

