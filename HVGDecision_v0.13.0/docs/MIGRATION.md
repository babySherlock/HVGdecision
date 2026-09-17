# Migration and reproducibility

Version 0.13.0 packages continuous-v0.3 fusion. It does not imply that historical benchmarks were run with this release.

- Use `hd.refine` for continuous fusion and `hd.audit_design` for weighted A and design diagnostics.
- `mode='auto'` is an alias for continuous fusion, not hard routing.
- Explicit endpoint tests may use `hd.ContinuousFusionConfig(weight_override=0.0)` or 1.0. Normal analysis leaves this unset.
- `setup_reference_query(...).find_best_hvg()` and the legacy CLI remain compatibility APIs, not the current continuous workflow.
- The cell component uses marker log effect 0.75 and replication fraction 0.90, matching the continuous-v0.3 wrapper.
- Both component engines are included. No historical package installation is required.
- The vendored v0.3 mathematics exclude external adapters. The percentile array is copied for pandas compatibility, without changing rank calculations.

## Historical replay

Keep the exact filtered cells, counts gene axis and ordering, Reference/Query panels, protected genes, settings and seed. Pass frozen ordered panels via `hvg_genes` and `query_hvg_genes`. Compare full evidence and gene sets, not deletion counts alone.

Audited Kang HVG fits used Reference span 0.5 and Query span 0.3. This package does not force three PBMC or ten Lung deletions. New inputs or changed settings require comparison before claiming equivalence.

Synthetic tests and wheel installation checks are not a rerun of the five-method, five-seed server benchmark. Record the distribution version, mathematical engine version and original analysis source separately.

