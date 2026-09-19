# Changelog

## 0.14.0

- Replaced replication-adequacy-dependent confidence weighting with dual-evidence signed-margin fusion.
- Added centered log-mean-exp fusion with default `tau_f = 0.10`.
- Changed the final decision boundary from a 0–1 fused confidence threshold to `unified_margin > 0`.
- Replication adequacy and design identifiability are now dataset-level audits only and do not enter gene-level fusion.
- Preserved component-specific biological protection by mapping a protected component to an effective margin of `-inf`.
- Added `UnifiedFusionConfig` / `DualEvidenceFusionConfig`.
- Added `centered_logmeanexp()` and `unified_margin_fusion()` to the public API.
- Kept the high-level `hd.refine(...)` workflow and Reference → Query panel transport.
- Added `result.design` as the current audit interface; `result.routing` remains only as a deprecated compatibility alias.
- Removed obsolete confidence-weighting fields from the public v0.14 fusion table.
- Added updated PBMC and Human Lung examples, tests and migration documentation.
