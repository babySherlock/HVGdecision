"""Edit the file path and metadata keys, then run in your analysis environment."""
import scanpy as sc
import hvgdecision as hd

DATA_PATH = '/path/to/data.h5ad'
BATCH_KEY = 'donor'
LABEL_KEY = 'cell_type'
OUTPUT_DIR = './HVGDecision_run01'

adata = sc.read_h5ad(DATA_PATH)
counts = hd.find_raw_counts(adata)
print(counts.audit.to_string(index=False))
if not counts.valid:
    raise ValueError(counts.error)

design = hd.audit_design(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY)
print(design.audit.to_string(index=False))

result = hd.refine(
    adata, batch_key=BATCH_KEY, label_key=LABEL_KEY,
    counts=counts, n_hvg=2000, output_dir=OUTPUT_DIR, return_details=True,
)
print(result)
print('Removed genes:', result.removed_genes)
print(result.decision_table.loc[result.decision_table.final_action.eq('remove')].to_string(index=False))
adata_hvg = result.adata
