# Migration from 0.9.1

Use `hd.refine(...)` for the current fixed-panel, automatically routed workflow.
The default is a 2000-gene independent Seurat v3 fit, not a search over budgets.
`hd.audit_design(...)` exposes the metadata-only routing audit.

The previous `setup_reference_query(...).find_best_hvg()` and `refine_hvg()`
methods remain available for historical analyses. They retain the 0.9.1
within-domain/cross-domain semantics and do NOT invoke the current automatic
branches. Their manuals are under `docs/legacy_0_9` and their example scripts
under `examples/legacy_0_9`. The old YAML CLI is `hvgdecision-legacy`.
Do not use it to describe the current manuscript algorithm.

The current CLI is `hvgdecision audit/refine --input ... --batch-key ...
--label-key ... --output ...`. The current `REFINEMENT_MODES` are auto,
cell_level and donor_aware. `VALID_MODES` and `normalize_mode` are retained only
for the legacy API.

Count validation was strengthened for both APIs: it now checks every stored
value, rejects non-finite or fractional values, and requires named external
tables to match cell IDs. DataFrames are no longer aligned by dimension alone.
Explicit matrices without IDs still require the caller to supply the correct
cell order. No existing result or source input is overwritten by the new API.

The current AnnData contains retained genes only. Replace code that searches
the final `.var` for deleted genes with `result.removed_genes` or the full
`adata.uns['hvgdecision']['gene_decisions']` audit. No `best_n_hvg` is reported
because a fixed base size is not an optimized recommendation.

The two notebook engines were extracted into self-contained modules. No R
installation, benchmark runner import or notebook monkeypatch is required.
The donor engine uses a per-call configuration and writable copies of arrays
for pandas 3 compatibility; these changes do not alter its formulas or gates.

One inherited unit test expected the legacy budget-search observed peak to
equal the recommended budget. It failed on the unmodified 0.9.1 source too.
The fixture now distinguishes the observed peak (6) from the smaller
biology-eligible recommendation (4); the legacy selection algorithm was not
changed to satisfy the test.
