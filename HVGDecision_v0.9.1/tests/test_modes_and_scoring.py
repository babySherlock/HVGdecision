import numpy as np
import pandas as pd
import pytest

from hvgdecision.modes import normalize_mode
from hvgdecision.scoring import three_domain_scores
from hvgdecision.statistics import compose_risk


def test_mode_aliases_are_explicit():
    assert normalize_mode("within") == "within_domain"
    assert normalize_mode("cross-technology") == "cross_domain"
    with pytest.raises(ValueError):
        normalize_mode(None)


def test_within_risk_formula_has_expected_signs():
    leakage = np.array([0.0, 1.0, 2.0])
    interaction = np.array([0.0, 1.0, 2.0])
    biology = np.array([2.0, 1.0, 0.0])
    raw, parts = compose_risk(leakage, interaction, biology)
    assert raw[-1] > raw[0]
    assert "risk_z" in parts


def test_three_domain_scoring_prefers_jointly_better_panel():
    frame = pd.DataFrame(
        {
            "dataset": ["D"] * 3,
            "integration_method": ["M"] * 3,
            "panel": ["base", "random", "hvgdecision"],
            "macro_f1": [0.5, 0.55, 0.70],
            "balanced_accuracy": [0.5, 0.52, 0.68],
            "rare_cell_macro_f1": [0.3, 0.35, 0.60],
            "dataset_mixing_score": [0.4, 0.45, 0.65],
            "donor_mixing_within_celltype": [0.5, 0.52, 0.66],
            "leiden_ari": [0.4, 0.45, 0.7],
            "leiden_nmi": [0.5, 0.53, 0.72],
            "celltype_silhouette_query": [0.2, 0.25, 0.5],
        }
    )
    scored = three_domain_scores(frame)
    top = scored.sort_values("overall_score").iloc[-1]
    assert top["panel"] == "hvgdecision"
    assert np.isclose(top["overall_score"], 1.0)
