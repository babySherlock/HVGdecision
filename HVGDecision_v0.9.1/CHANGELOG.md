# Changelog

## 0.9.1 — explicit dual-mode release

- Added mandatory algorithm `mode` selection: `within_domain` or `cross_domain`.
- Replaced the legacy Stage-2 maxT/consistency risk with the manuscript Within-domain V1 rule: `Z(leakage) + 0.75 Z(interaction) - Z(biology)`.
- Added per-gene conditional donor-label permutation with BH-FDR and bootstrap recurrence gates.
- Added Cross-domain Rule V3: Reference risk percentile + label-free Reference/Query mean/detection shift + Reference biology/marker protection + fixed top-k removal.
- Added explicit `cross_domain_delete_budget` (primary audited values 5/10/20; default 5 in the Python convenience workflow).
- Added `study.run()` to route the selected mode while retaining `find_best_hvg()` for within-domain budget search and `refine_hvg()` for frozen external panels.
- Query true labels remain excluded from cross-domain feature selection; Query expression is used only for Query HVG selection and label-free shift evidence.
- Added `three_domain_scores()` for optional scIB-inspired Transfer / Batch correction / Biological fidelity / Overall benchmark aggregation.
- Added mode-specific YAML examples, Methods and output-schema documentation.
- Added wheel-first installation instructions for normal users.
- Added automatic two-AnnData cross-domain input preparation with independent raw-count detection and shared-gene alignment.
- Added direct local HVG CSV/TSV/TXT tutorials; headerless one-gene-per-line TXT is supported.
- Reworked GitHub Markdown equations to use display-math blocks that do not get parsed as hyperlinks.

## 0.8.0 - 2026-08-29

- Changed HVG-budget selection to a biology-first minimum-sufficient rule.
- Added donor-held-out rare-cell macro-F1 and minimum per-class recall.
- Added a composite `biology_sufficiency_score`.
- Replaced stability-near-peak gating with a robust stability quality floor.
- Replaced fixed-fraction donor rejection with a robust donor-instability upper-outlier guardrail.
- Added a noise-adaptive biology plateau with conservative 0.01–0.015 default bounds.
- Kept Query expression, Query labels, and downstream integration fully excluded from selection.
- Kept old public arguments accepted for compatibility; the old 0.03 biology tolerance can no longer widen the v0.8 default plateau.

## 0.8.0 - 2026-08-28

- Replaced stability-only HVG-budget recommendation with a joint Reference-only minimum-sufficient rule.
- Added nested leave-one-Reference-donor cell-type transfer using macro-F1 and balanced accuracy.
- The HVG panel used for each held-out Reference donor is refit without that donor.
- Added a within-biology donor-centroid instability guardrail that conditions on all declared `biology_keys`.
- A small budget such as 500 is now eligible only when stability, donor-held-out biology preservation, and the donor guardrail all pass.
- Added explicit `selection_status` and a minimum-constraint-violation fallback when the three gates have no common plateau.
- Added `BudgetSearchResult.budget_summary` plus per-budget/per-holdout audit CSVs.
- Query expression, Query labels, and downstream integration performance remain excluded from HVG-budget selection.

## 0.6.0 - 2026-08-28

- Refocused the public workflow on HVG recommendation and harmful-gene detection before integration.
- `find_best_hvg()` no longer runs Harmony/BBKNN/Scanorama/scVI/scANVI internally.
- HVG budget selection now uses Reference-only leave-one-batch-out HVG reproducibility.
- Query expression and Query labels are excluded from budget selection and harmful-gene estimation.
- Final AnnData now retains all genes and marks the final panel with `var["highly_variable"]`.
- Raw counts are guaranteed in `layers["counts"]` in the returned object.
- Added separate `risk_flagged`, marker-protection, `harmful`, `removed`, and final-HVG annotations.
- Added `result.harmful_genes`, `result.final_hvg_genes`, and `result.recommended_n_hvg` convenience properties.
- Added `biology_keys`/`protected_genes` examples for CD45+/CD45- and other protected biological contrasts.

## 0.3.0 — 2026-08-28

- changed `find_best_hvg()` default to auditable `selection_mode="auto"`; Seurat-v3 budgets are independently refit instead of being approximated by tail truncation;
- added budget-composition audits showing overlap, replacements, and exact prefix nesting;
- added public `study.refine_hvg(...)` for frozen R/Python HVG lists and tables without re-running feature selection;
- added automatic nested vs non-nested routing, `control_mode=auto/independent/truncate`, and `run_refinement=True/False` override;
- added external independent same-N controls without fabricating a Scanpy surrogate for an R selector;
- added gene-level decisions for removed, retained, marker-protected, missing, and finally adopted genes;
- expanded default HVG budgets through 5000;
- returned a downstream-ready AnnData with raw `layers["counts"]` and compact `uns["hvgdecision"]`;
- reduced default output clutter by disabling decorative budget plots and large diagnostic long tables unless requested.

## 0.2.1 — 2026-08-28

- added a complete raw-count source audit across all layers, `adata.raw.X`, and `adata.X`;
- added explicit user-source validation without silent fallback;
- added external DataFrame, CSV/TSV, AnnData/Raw, dense/sparse matrix, and `(matrix, gene_names)` support;
- made the selected AnnData directly consumable by scVI and CPU integration engines.

## 0.2.0 — 2026-08-28

- Added a Scanpy-style in-memory `AnnData` API.
- Added one-call Reference/Query setup and one-call HVG-budget screening.
- Added automatic raw-count source detection and strict validation.
- Added resumable budget-by-engine benchmarking and smallest-near-optimal selection.

## 0.1.0 — 2026-08-28

- Added input and raw-count audits.
- Added smallest-near-optimal HVG budget recommendation.
- Added Reference-only conditional maxT risk estimation.
- Added stratified stability bootstrap and marker hard protection.
- Added original, refined, direct-same-N, random-matched, and rank-tail panels.
- Added optional Harmony, BBKNN, Scanorama, scVI, and scANVI benchmark engines.
- Added incremental recovery, paired statistics, and complete audit outputs.

## 0.4.0 - 2026-08-28

- `find_best_hvg()` now continues from budget selection into Reference-only marker-protected risk refinement by default.
- Python-route results expose `removed_genes`, `retained_genes`, `final_genes`, `final_n_hvg`, and `decision_table`; zero deletion is explicitly supported.
- The exact recommended budget panel is frozen before risk scoring.
- Internal Seurat-v3 workflows can independently refit matched-N controls after targeted deletion; external R panels are never replaced by a Scanpy surrogate.
- Default integration engines are Harmony + BBKNN; scVI/scANVI remain optional.
- Added two explicit runnable examples for Python-selected and externally imported HVGs.