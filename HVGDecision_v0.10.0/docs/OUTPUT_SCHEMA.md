# Output schema

`refine(..., return_details=True)` returns `RefinementResult`; otherwise it
returns the same `.adata` directly. `result.save(path)` requires a new/empty
directory and preserves existing results.

| File | Content |
| --- | --- |
| adata_hvg.h5ad | All input cells, retained genes, raw X and counts layer |
| gene_decisions.csv | All base-panel genes, ranks, original branch evidence, protection, final action |
| removed_genes.csv | Removed rows, with headers even if empty |
| retained_genes.csv | Retained rows in input-panel rank order |
| counts_audit.csv | Candidate-source numeric checks and selected source |
| routing_summary.csv | A, w, rank, iota, graph components, automatic mode and reason |
| donor_celltype_coverage.csv | All donor/type strata, cell counts, support and supported-edge flag |
| routing_by_donor.csv | Donor cell counts and mean support |
| routing_by_celltype.csv | Type counts, effective donor support and type weights |
| run_manifest.json | Actual mode, override status, settings, hashes, cohort and software versions |

Common gene columns: `gene`, `input_rank` (one-based), `selected_mode`, `harmful`
(the actual removal call), `explicitly_protected`, `final_action` (`keep`/`remove`),
`decision_reason`, `in_final_panel`. `gene` is the exact input count-axis ID, not
an inferred symbol. Cell-level evidence includes `biology_z`, `donor_leakage_z`,
`interaction_range_z`, `risk_score_z`, permutation p/FDR/maxT p and bootstrap
fraction. Donor-aware evidence includes `celltype_eta2`, technical and biological
support scores, Q risk score, direction agreement, dominance, LODO metrics and
the realized quantile thresholds. See METHODS.md for definitions.

For `insufficient_confounded`, genes are kept with
`not_tested_insufficient_confounded_design`. Absent evidence means untested,
not a zero risk value. `routing_summary.csv` always records the automatic route;
`run_manifest.json` records the actual mode if a branch override was requested.

`pass_*` flags and `failed_rules` list the removal conditions a gene did not
meet, including protective conditions. Common symbol/Ensembl annotation columns
are copied from the selected AnnData gene metadata when available. No online
lookup or guessed ID mapping is performed; missing annotations remain empty.

AnnData `.var` describes retained genes only. `.uns['hvgdecision']` includes
`gene_decisions`, `routing_summary`, `base_genes`, `removed_genes`,
`selection_cell_ids`, `selected_mode` and `manifest_json`. The latter is JSON
text so nested parameter structures round-trip across h5ad versions.
Selection cell IDs are local audit data: review them before public sharing.

No PCA, UMAP, integration graph, accuracy score or benchmark claim is stored.
Original observation metadata is retained; old embeddings and layers are not.
