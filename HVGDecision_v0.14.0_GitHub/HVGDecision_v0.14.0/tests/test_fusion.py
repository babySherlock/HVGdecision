import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import hvgdecision as hd


class DummyAdata:
    def __init__(self, obs):
        self.obs = obs
        self.n_obs = len(obs)


def test_centered_logmeanexp_symmetry():
    md = np.array([0.2, -0.4, 1.0, -np.inf])
    mc = np.array([-0.1, 0.3, 1.0, 0.2])
    a = hd.centered_logmeanexp(md, mc, tau=0.10)
    b = hd.centered_logmeanexp(mc, md, tau=0.10)
    assert np.allclose(a, b, equal_nan=True)


def test_equal_margins_are_preserved():
    m = np.array([-2.0, -0.2, 0.0, 0.4, 2.0])
    u = hd.centered_logmeanexp(m, m, tau=0.10)
    assert np.allclose(u, m)


def test_single_component_penalty():
    u = hd.centered_logmeanexp(
        np.array([0.2]),
        np.array([-np.inf]),
        tau=0.10,
    )
    expected = 0.2 - 0.10 * np.log(2.0)
    assert np.allclose(u, expected)


def test_invalid_tau():
    with pytest.raises(ValueError):
        hd.UnifiedFusionConfig(tau=0.0)
    with pytest.raises(ValueError):
        hd.centered_logmeanexp(np.array([0.0]), np.array([0.0]), tau=-1.0)


def test_design_audit_has_no_fusion_weight():
    rows = []
    for donor in ["d1", "d2", "d3"]:
        for celltype in ["A", "B", "C"]:
            rows.extend([(donor, celltype)] * 20)
    obs = pd.DataFrame(rows, columns=["donor", "celltype"])
    obs.index = [f"c{i}" for i in range(len(obs))]
    result = hd.audit_design(
        DummyAdata(obs),
        batch_key="donor",
        label_key="celltype",
    )
    assert result.summary["identifiable"] is True
    assert result.summary["replication_adequacy"] > 0
    assert result.summary["replication_adequacy_role"] == "audit_metadata_only_not_used_in_fusion"
    assert result.summary["component_selection_used"] is False
    assert result.summary["dataset_dependent_component_weight_used"] is False
    assert "donor_aware_weight" not in result.summary


def _synthetic_evidence(n=60):
    rng = np.random.default_rng(14)
    genes = [f"g{i}" for i in range(n)]

    donor = pd.DataFrame(
        {
            "gene": genes,
            "nuisance_eta2_within_celltype": rng.uniform(0, 1, n),
            "single_group_effect_dominance": rng.uniform(0, 1, n),
            "celltype_effect_heterogeneity": rng.uniform(0, 1, n),
            "lodo_instability": rng.uniform(0, 1, n),
            "celltype_effect_direction_agreement": rng.uniform(0, 1, n),
            "celltype_eta2": rng.uniform(0, 1, n),
            "replicate_count": rng.integers(3, 7, n),
            "lodo_positive_fraction": rng.uniform(0, 1, n),
            "replicated_marker_protection": np.zeros(n, dtype=bool),
            "rare_marker_protection": np.zeros(n, dtype=bool),
            "user_protected": np.zeros(n, dtype=bool),
        }
    )

    cell = pd.DataFrame(
        {
            "gene": genes,
            "donor_leakage": rng.uniform(0, 1, n),
            "interaction_range": rng.uniform(0, 2, n),
            "biology_eta2": rng.uniform(0, 1, n),
            "permutation_fdr": rng.uniform(0, 0.2, n),
            "bootstrap_risk_fraction": rng.uniform(0, 1, n),
            "marker_protected": np.zeros(n, dtype=bool),
            "explicitly_protected": np.zeros(n, dtype=bool),
        }
    )
    return donor, cell


def test_unified_fusion_public_columns_and_audit_only_A():
    donor, cell = _synthetic_evidence()
    table, manifest = hd.unified_margin_fusion(
        donor,
        cell,
        replication_adequacy=3.56,
        config=hd.UnifiedFusionConfig(tau=0.10),
    )
    assert len(table) == len(donor)
    for col in [
        "donor_branch_margin_raw",
        "cell_branch_margin_raw",
        "donor_branch_margin_effective",
        "cell_branch_margin_effective",
        "unified_margin",
        "discovery_remove",
        "final_action",
    ]:
        assert col in table
    assert not any("confidence" in c for c in table.columns)
    assert manifest["replication_adequacy_role"] == "audit_metadata_only_not_used_in_fusion"
    assert manifest["dataset_dependent_component_weight_used"] is False
    assert manifest["component_selection_used"] is False


def test_component_specific_protection_sets_effective_margin_to_minus_inf():
    donor, cell = _synthetic_evidence()
    donor.loc[0, "user_protected"] = True
    table, _ = hd.unified_margin_fusion(donor, cell, config=hd.UnifiedFusionConfig())
    row = table.set_index("gene").loc["g0"]
    assert np.isneginf(row["donor_branch_margin_effective"])


