from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import hvgdecision as hd
from hvgdecision import api


def _fake_refine_no_deletion(config):
    output = Path(config["output_dir"]) / "01_refine"
    output.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(config["base_gene_table"]).sort_values("gene_rank")
    risk = pd.DataFrame(
        {
            "gene": base["gene"].astype(str),
            "risk_score": np.zeros(len(base), dtype=float),
            "passes_maxT": np.zeros(len(base), dtype=bool),
            "passes_stability": np.zeros(len(base), dtype=bool),
            "passes_effect_floors": np.zeros(len(base), dtype=bool),
            "selected_risk_gene": np.zeros(len(base), dtype=bool),
            "hard_replicated_marker_protection": np.zeros(len(base), dtype=bool),
        }
    )
    risk.to_csv(output / "gene_risk_evidence.csv", index=False)
    return {
        "empirical_maxT_threshold": 1.0,
        "stability_design": "test",
    }




def _fake_reference_biology_transfer(adata, sampled_positions, genes, config, reference_batches):
    return pd.DataFrame(
        [
            {
                "held_out_reference_batch": str(batch),
                "biology_eval_cells_total": 6,
                "biology_eval_cells_known": 6,
                "biology_label_coverage": 1.0,
                "biology_macro_f1": 0.90,
                "biology_balanced_accuracy": 0.90,
                "biology_rare_macro_f1": 0.88,
                "biology_rare_eval_cells": 2,
                "biology_rare_label_count": 1,
                "biology_min_class_recall": 0.85,
                "biology_transfer_score": 0.90,
                "biology_sufficiency_score": 0.8825,
            }
            for batch in reference_batches
        ]
    )


def _fake_donor_instability(adata, sampled_positions, genes, config):
    return 0.10

def synthetic_adata() -> ad.AnnData:
    rng = np.random.default_rng(12)
    counts = rng.poisson(2.0, size=(18, 10)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "donor": np.repeat(["D1", "D2", "D3"], 6),
            "cell_type": ["A", "A", "A", "B", "B", "B"] * 3,
        },
        index=[f"cell_{index}" for index in range(18)],
    )
    var = pd.DataFrame(index=[f"gene_{index}" for index in range(10)])
    return ad.AnnData(X=counts, obs=obs, var=var)


def test_setup_reference_query_accepts_scanpy_anndata(tmp_path: Path):
    adata = synthetic_adata()
    study = hd.setup_reference_query(
        adata,
        mode="within_domain",
        batch_key="donor",
        label_key="cell_type",
        reference=["D1", "D2"],
        query="D3",
        output_dir=tmp_path,
    )
    assert study.reference_mask.sum() == 12
    assert study.query_mask.sum() == 6
    assert adata.obs["_hvgdecision_role"].value_counts()["reference"] == 12
    assert study.config["counts_layer"] is None


def test_find_best_hvg_uses_reference_stability_not_integration(tmp_path: Path, monkeypatch):
    study = hd.setup_reference_query(
        synthetic_adata(),
        mode="within_domain",
        batch_key="donor",
        label_key="cell_type",
        reference=["D1", "D2"],
        query="D3",
        output_dir=tmp_path,
    )

    def fake_ranking(adata, mask, config, n_top):
        # Full Reference = 12 cells. Leave-one-Reference-out = 6 cells.
        genes = [f"gene_{index}" for index in range(n_top)]
        if int(mask.sum()) == 6 and n_top == 4:
            genes = ["gene_0", "gene_1", "gene_8", "gene_9"]
        return pd.DataFrame({"gene": genes, "gene_rank": np.arange(1, n_top + 1)})

    monkeypatch.setattr(api, "_hvg_ranking", fake_ranking)
    monkeypatch.setattr(api, "refine_panel", _fake_refine_no_deletion)
    monkeypatch.setattr(api, "_reference_biology_transfer", _fake_reference_biology_transfer)
    monkeypatch.setattr(api, "_within_biology_donor_instability", _fake_donor_instability)

    result = study.find_best_hvg(
        budgets=[4, 6],
        hvg_method="seurat_v3",
        selection_mode="auto",
        return_details=True,
    )
    assert result.best_n_hvg == 6
    assert result.harmful_genes == []
    assert result.final_n_hvg == 6
    assert result.adata.n_vars == 10  # full gene space is retained
    assert int(result.adata.var["highly_variable"].sum()) == 6
    assert "counts" in result.adata.layers
    assert (tmp_path / "00_budget_search" / "reference_hvg_stability_summary.csv").exists()
    assert (tmp_path / "01_hvg_decision" / "hvg_decision_table.csv").exists()


def test_budget_composition_audit_distinguishes_prefix_and_replacement():
    prefix = {
        4: pd.DataFrame({"gene": ["a", "b", "c", "d"], "gene_rank": [1, 2, 3, 4]}),
        3: pd.DataFrame({"gene": ["a", "b", "c"], "gene_rank": [1, 2, 3]}),
    }
    audit = api._budget_composition_audit(prefix)
    row = audit.loc[audit["n_hvg"].eq(3)].iloc[0]
    assert bool(row["strict_prefix_truncation"])

    replaced = {
        4: prefix[4],
        3: pd.DataFrame({"gene": ["a", "x", "c"], "gene_rank": [1, 2, 3]}),
    }
    audit = api._budget_composition_audit(replaced)
    row = audit.loc[audit["n_hvg"].eq(3)].iloc[0]
    assert not bool(row["strict_prefix_truncation"])
    assert row["new_vs_largest_prefix"] == 1


def test_find_best_hvg_auto_refits_seurat_v3_for_full_and_holdouts(tmp_path: Path, monkeypatch):
    study = hd.setup_reference_query(
        synthetic_adata(),
        mode="within_domain",
        batch_key="donor",
        label_key="cell_type",
        reference=["D1", "D2"],
        query="D3",
        output_dir=tmp_path,
    )
    calls = []

    def fake_ranking(adata, mask, config, n_top):
        calls.append((int(mask.sum()), n_top))
        genes = [f"gene_{index}" for index in range(n_top)]
        return pd.DataFrame({"gene": genes, "gene_rank": np.arange(1, n_top + 1)})

    monkeypatch.setattr(api, "_hvg_ranking", fake_ranking)
    monkeypatch.setattr(api, "refine_panel", _fake_refine_no_deletion)
    monkeypatch.setattr(api, "_reference_biology_transfer", _fake_reference_biology_transfer)
    monkeypatch.setattr(api, "_within_biology_donor_instability", _fake_donor_instability)

    result = study.find_best_hvg(
        budgets=[4, 6],
        hvg_method="seurat_v3",
        selection_mode="auto",
        return_details=True,
    )
    # 1 full-Reference + 2 leave-one-Reference-batch-out fits per budget.
    assert sum(n_top == 4 for _, n_top in calls) == 3
    assert sum(n_top == 6 for _, n_top in calls) == 3
    assert result.adata.uns["hvgdecision"]["actual_selection_mode"] == "independent"
    assert result.adata.uns["hvgdecision"]["integration_run_during_selection"] is False


def test_external_hvg_input_preserves_frozen_order():
    frame = pd.DataFrame(
        {
            "method": ["Seurat_v3"] * 3,
            "n_hvg": [3, 3, 3],
            "gene_rank": [2, 1, 3],
            "gene": ["B", "A", "C"],
        }
    )
    normalized = api._read_external_hvg_input(
        frame,
        gene_column="gene",
        rank_column="gene_rank",
        budget_column="n_hvg",
        method_column="method",
        method_name="Seurat_v3",
    )
    assert normalized.sort_values("gene_rank")["gene"].tolist() == ["A", "B", "C"]


def test_budget_search_result_exposes_harmful_and_final_aliases():
    result = api.BudgetSearchResult(
        recommendation=pd.DataFrame(
            [{"recommended_minimum_sufficient_hvg": 3, "observed_peak_hvg": 3}]
        ),
        results=pd.DataFrame([{"stability_score": 1.0}]),
        gene_table=pd.DataFrame(
            {
                "n_hvg": [3, 3, 3],
                "gene_rank": [1, 2, 3],
                "gene": ["A", "B", "C"],
            }
        ),
        composition_audit=pd.DataFrame(),
        decision_table=pd.DataFrame(
            {
                "gene": ["A", "B", "C"],
                "input_rank": [1, 2, 3],
                "harmful_gene": [False, True, False],
                "final_action": ["retain", "remove", "retain"],
                "in_final_panel": [True, False, True],
            }
        ),
        refinement_audit={},
        output_dir=Path("."),
        adata=None,
    )
    assert result.recommended_n_hvg == 3
    assert result.harmful_genes == ["B"]
    assert result.final_hvg_genes == ["A", "C"]
    assert result.final_n_hvg == 2


def test_joint_budget_rule_rejects_stable_but_biology_insufficient_500():
    summary = pd.DataFrame(
        {
            "n_hvg": [500, 1000, 2000],
            "stability_score_mean": [0.99, 0.98, 0.97],
            "biology_transfer_score_mean": [0.60, 0.90, 0.91],
            "within_biology_donor_instability_mean": [0.08, 0.07, 0.075],
        }
    )
    audited, recommendation = api._recommend_joint_reference_budget(
        summary,
        stability_tolerance=0.03,
        biology_tolerance=0.03,
        donor_guardrail_quantile=0.75,
    )
    assert int(recommendation.iloc[0]["recommended_minimum_sufficient_hvg"]) == 1000
    row500 = audited.loc[audited["n_hvg"].eq(500)].iloc[0]
    assert bool(row500["stability_sufficient"])
    assert not bool(row500["biology_sufficient"])
    assert not bool(row500["joint_sufficient"])


def test_v08_budget_rule_can_allow_500_when_it_is_truly_on_biology_plateau():
    summary = pd.DataFrame(
        {
            "n_hvg": [500, 1000, 2000],
            "stability_score_mean": [0.99, 0.985, 0.98],
            "biology_transfer_score_mean": [0.905, 0.91, 0.91],
            "biology_sufficiency_score_mean": [0.905, 0.91, 0.91],
            "within_biology_donor_instability_mean": [0.06, 0.07, 0.08],
        }
    )
    _, recommendation = api._recommend_joint_reference_budget(
        summary,
        stability_tolerance=0.02,
        biology_tolerance=0.03,
        donor_guardrail_quantile=0.75,
    )
    assert int(recommendation.iloc[0]["recommended_minimum_sufficient_hvg"]) == 500


def test_v08_pbmc_profile_moves_recommendation_from_500_to_1500():
    summary = pd.DataFrame(
        {
            "n_hvg": [
                500, 750, 1000, 1250, 1500, 1750, 1800, 1900, 2000,
                2250, 2500, 2750, 3000, 3500, 4000, 4500, 5000,
            ],
            "stability_score_mean": [
                0.897252, 0.880622, 0.877783, 0.872864, 0.869367,
                0.869111, 0.868992, 0.869956, 0.869927, 0.871321,
                0.873218, 0.876915, 0.879776, 0.882031, 0.889050,
                0.895905, 0.898581,
            ],
            "biology_transfer_score_mean": [
                0.881896, 0.887738, 0.889063, 0.887047, 0.894269,
                0.894561, 0.894030, 0.897123, 0.897493, 0.895770,
                0.904194, 0.902512, 0.878144, 0.899899, 0.872435,
                0.874518, 0.873082,
            ],
            "within_biology_donor_instability_mean": [
                0.041636, 0.045249, 0.049899, 0.053590, 0.050393,
                0.051216, 0.051830, 0.051866, 0.051570, 0.051332,
                0.052901, 0.053988, 0.050852, 0.052828, 0.056308,
                0.056231, 0.058103,
            ],
        }
    )
    audited, recommendation = api._recommend_joint_reference_budget(
        summary,
        stability_tolerance=0.02,
        biology_tolerance=0.03,
        donor_guardrail_quantile=0.75,
    )
    assert int(recommendation.iloc[0]["recommended_minimum_sufficient_hvg"]) == 1500
    row500 = audited.loc[audited["n_hvg"].eq(500)].iloc[0]
    row1500 = audited.loc[audited["n_hvg"].eq(1500)].iloc[0]
    assert not bool(row500["biology_sufficient"])
    assert bool(row1500["biology_sufficient"])
    assert bool(row1500["stability_sufficient"])
