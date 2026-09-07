"""Step 0: inspect the automatically selected raw-count source."""
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")
check = hd.find_raw_counts(adata)
print("valid:", check.valid)
print("selected:", check.location)
print(check.audit.to_string(index=False))
