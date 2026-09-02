# Output schema — HVGDecision 0.9.0

## Common outputs

`01_refine/base_hvg_ranking.csv`
: Base panel that will be refined. In `within_domain` this is normally a Reference panel; in `cross_domain` it is the Query panel.

`01_refine/gene_risk_evidence.csv`
: One row per base-panel gene with auditable risk fields.

`01_refine/selected_risk_genes.csv`
: Genes selected for deletion by the active mode.

`01_refine/refinement_audit.json`
: Method version, selected mode, split, parameters, leakage discipline, removed genes and final panel size.

Python API additionally writes `hvg_decision_table.csv`, `final_hvg_genes.csv`, `harmful_genes.csv` and `decision.json` in the decision directory.

## Within-domain key columns

```text
gene
base_hvg_rank
donor_leakage
interaction_instability
biology_eta2
donor_leakage_z
interaction_instability_z
biology_z
risk_score_raw
risk_z
permutation_p
permutation_fdr
bootstrap_pass_fraction
passes_permutation_fdr
passes_effect_floors
passes_stability
risk_flagged_before_protection
marker_replication_fraction
hard_replicated_marker_protection
selected_risk_gene
```

## Cross-domain additional outputs

```text
QUERY_hvg_ranking.csv
REFERENCE_hvg_ranking.csv
REFERENCE_RULE_V1_full_audit_for_crossdomain.csv
CROSSDOMAIN_reference_query_shift_audit.csv
CROSSDOMAIN_RULE_V3_ranking.csv
CROSSDOMAIN_RULE_V3_budget_gene_membership.csv
CROSSDOMAIN_RULE_V3_selection_manifest.csv
CROSSDOMAIN_RULE_V3_removal_audit.csv
```

Cross-domain key columns include:

```text
reference_rule_percentile
mean_shift_effect
mean_shift_percentile
detection_rate_shift
detection_shift_percentile
dataset_shift_score
reference_biology_percentile
marker_replication_fraction_v3
biology_protection_score_v3
crossdomain_hard_protected
technical_consensus
rule_shift_agreement
cross_risk_score
cross_risk_rank
```

## Final AnnData annotations

The full count-source gene space is retained. The final feature panel is marked by:

```text
var["highly_variable"]
var["hvgdecision_candidate"]
var["hvgdecision_risk_score"]
var["hvgdecision_risk_flagged"]
var["hvgdecision_marker_protected"]
var["hvgdecision_harmful"]
var["hvgdecision_removed"]
var["hvgdecision_final"]
var["hvgdecision_reason"]
```

Raw counts are available in `layers["counts"]` in the returned object.
