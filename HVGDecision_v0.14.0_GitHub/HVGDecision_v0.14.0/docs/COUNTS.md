# Counts inputs

```python
counts = hd.find_raw_counts(adata)                        # automatic audit
counts = hd.find_raw_counts(adata, source='counts')       # named layer
counts = hd.find_raw_counts(adata, source='raw.X')        # explicitly raw.X
counts = hd.find_raw_counts(adata, source='X')            # explicitly X
counts = hd.find_raw_counts(adata, source='counts.csv')   # first column is cell ID
```

Check `counts.valid`, `counts.location` and `counts.audit`. An explicit invalid source is not replaced automatically. Integer-like normalized data could pass a numerical check; establish raw-count provenance independently.

```python
import pandas as pd
frame = pd.read_csv('counts.csv', index_col=0)
counts = hd.find_raw_counts(adata, source=frame)
if not counts.valid:
    raise ValueError(counts.audit.to_string(index=False))
```

Named tables align by cell identifiers; genes remain those of the supplied count table. A gene-by-cell table can be transposed automatically when column IDs match the input cells. Duplicated IDs or missing input cells are errors. For an unnamed matrix, provide explicit axes:

```python
counts = hd.find_raw_counts(adata, source={
    'matrix': count_matrix,
    'obs_names': count_cell_ids,
    'gene_names': count_gene_ids,
})
```

Pass the resulting object to `hd.refine(counts=counts, ...)`. Recompute the count-source audit after subsetting or reordering AnnData. The refinement API validates counts again and rejects zero-library cells. It does not add a pseudocount or round normalized expression into counts.
