"""Signed-margin gate reconstruction used by HVGDecision v0.14.0.

This private module reconstructs the validated donor-aware and cell-level
risk gates and converts them to standardized signed margins. It does not
perform the final dual-evidence fusion.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence
import hashlib
import inspect
import math

import numpy as np
import pandas as pd


__version__ = "0.14.0"


# =====================================================================
# Configuration
# =====================================================================

@dataclass(frozen=True)
class DonorGateConfig:
    # Historical donor-aware gates.
    technical_quantile: float = 0.975
    risk_quantile: float = 0.985
    minimum_absolute_risk: float = 0.30
    maximum_biological_support_quantile: float = 0.45
    minimum_replicates: int = 3
    minimum_single_donor_dominance: float = 0.55
    lodo_positive_fraction_ceiling: float = 0.80
    minimum_direction_agreement: float = 0.75


@dataclass(frozen=True)
class CellGateConfig:
    # Historical cell-level gates.
    alpha: float = 0.05
    leakage_z_floor: float = 1.0
    risk_z_floor: float = 1.0
    biology_z_ceiling: float = 0.0
    bootstrap_pass_fraction: float = 0.80



# =====================================================================
# Basic utilities
# =====================================================================

def strict_bool(values) -> np.ndarray:
    s = pd.Series(values).copy()

    def one(v):
        if pd.isna(v):
            raise ValueError(
                "Missing Boolean evidence; not silently converted to False"
            )
        if v is True or v == 1 or str(v).strip().lower() == "true":
            return True
        if v is False or v == 0 or str(v).strip().lower() == "false":
            return False
        raise ValueError(f"Invalid Boolean value: {v!r}")

    return s.map(one).to_numpy(dtype=bool)


def _numeric(table: pd.DataFrame, name: str) -> np.ndarray:
    if name not in table.columns:
        raise KeyError(f"Missing evidence column: {name}")
    return pd.to_numeric(
        table[name],
        errors="raise",
    ).to_numpy(dtype=float)


def percentile_rank(values) -> np.ndarray:
    x = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    )
    r = x.rank(
        method="average",
        pct=True,
    ).to_numpy(dtype=float, copy=True)
    r[~np.isfinite(r)] = 0.0
    return r


def robust_z(values) -> np.ndarray:
    x = np.asarray(values, dtype=float)

    med = np.nanmedian(x)
    scale = 1.4826 * np.nanmedian(
        np.abs(x - med)
    )

    if not np.isfinite(scale) or scale < 1e-8:
        scale = np.nanstd(x)

    if not np.isfinite(scale) or scale < 1e-8:
        scale = 1.0

    return (x - med) / scale


def robust_scale(values, floor: float = 1e-8) -> float:
    """
    Robust scale used ONLY to put gate distances on comparable units.
    It does not alter the sign of a gate margin.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) == 0:
        return 1.0

    med = np.median(x)
    scale = 1.4826 * np.median(
        np.abs(x - med)
    )

    if not np.isfinite(scale) or scale < floor:
        scale = np.std(x)

    if not np.isfinite(scale) or scale < floor:
        scale = 1.0

    return float(scale)


def quantile(values, q: float) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) == 0:
        raise ValueError(
            "Cannot calculate a gate quantile from zero finite values"
        )

    return float(
        np.quantile(x, float(q))
    )


def sigmoid(values):
    x = np.asarray(values, dtype=float)
    out = np.empty_like(x, dtype=float)

    pos = x >= 0
    out[pos] = 1.0 / (
        1.0 + np.exp(-x[pos])
    )

    ex = np.exp(x[~pos])
    out[~pos] = ex / (1.0 + ex)

    return out



def _signed_margin(
    raw_margin,
    hard_pass,
    *,
    scale: float,
    epsilon: float,
) -> np.ndarray:
    """
    Convert a raw gate distance into a standardized signed margin.

    Positive = historical gate passes.
    Negative = historical gate fails.

    Equality handling is corrected with +/- epsilon so that the sign
    reproduces the exact historical >= / <= / < semantics.
    """
    raw = np.asarray(
        raw_margin,
        dtype=float,
    )
    hp = np.asarray(
        hard_pass,
        dtype=bool,
    )

    s = float(scale)
    if not np.isfinite(s) or s <= 0:
        s = 1.0

    m = raw / s

    near_zero = (
        np.isfinite(m)
        & (np.abs(m) <= float(epsilon))
    )

    m[
        near_zero & hp
    ] = float(epsilon)

    m[
        near_zero & (~hp)
    ] = -float(epsilon)

    # Any non-finite gate is treated as failing.
    m[~np.isfinite(m)] = -np.inf

    return m


def _soft_and(
    margins: Sequence[np.ndarray],
) -> np.ndarray:
    """
    Continuous AND surrogate with exact sign semantics:
    min(margins) > 0 iff all signed margins > 0.
    """
    return np.min(
        np.vstack(margins),
        axis=0,
    )


def _soft_or(
    margins: Sequence[np.ndarray],
) -> np.ndarray:
    """
    Continuous OR surrogate with exact sign semantics:
    max(margins) > 0 iff at least one signed margin > 0.
    """
    return np.max(
        np.vstack(margins),
        axis=0,
    )


# =====================================================================
# Donor-aware evidence: hard gates -> standardized signed margin
# =====================================================================

def donor_soft_gate(
    evidence: pd.DataFrame,
    *,
    config: Optional[DonorGateConfig] = None,
    soft_temperature: float = 1.0,
    boundary_epsilon: float = 1e-9,
    scale_floor: float = 1e-8,
) -> pd.DataFrame:
    cfg = config or DonorGateConfig()

    if float(soft_temperature) <= 0:
        raise ValueError(
            "soft_temperature must be > 0"
        )

    t = evidence.copy().reset_index(
        drop=True
    )

    if "gene" not in t.columns:
        raise KeyError(
            "donor evidence lacks gene"
        )

    if (
        t["gene"].isna().any()
        or t["gene"].duplicated().any()
    ):
        raise ValueError(
            "donor evidence must contain one unique row per gene"
        )

    nuisance = percentile_rank(
        _numeric(
            t,
            "nuisance_eta2_within_celltype",
        )
    )

    dominance_pct = percentile_rank(
        _numeric(
            t,
            "single_group_effect_dominance",
        )
    )

    heterogeneity = percentile_rank(
        _numeric(
            t,
            "celltype_effect_heterogeneity",
        )
    )

    lodo_instability = percentile_rank(
        _numeric(
            t,
            "lodo_instability",
        )
    )

    technical = (
        0.35 * nuisance
        + 0.25 * dominance_pct
        + 0.20 * heterogeneity
        + 0.20 * lodo_instability
    )

    direction = _numeric(
        t,
        "celltype_effect_direction_agreement",
    )

    biology = (
        0.55
        * percentile_rank(
            _numeric(
                t,
                "celltype_eta2",
            )
        )
        + 0.45 * direction
    )

    risk = technical - biology

    technical_threshold = quantile(
        technical,
        cfg.technical_quantile,
    )

    risk_threshold = max(
        float(
            cfg.minimum_absolute_risk
        ),
        quantile(
            risk,
            cfg.risk_quantile,
        ),
    )

    biology_ceiling = quantile(
        biology,
        cfg.maximum_biological_support_quantile,
    )

    replicate_count = _numeric(
        t,
        "replicate_count",
    )

    dominance = _numeric(
        t,
        "single_group_effect_dominance",
    )

    lodo_positive = _numeric(
        t,
        "lodo_positive_fraction",
    )

    # Exact historical gate Booleans.
    pass_technical = (
        technical >= technical_threshold
    )

    pass_risk = (
        risk >= risk_threshold
    )

    pass_biology = (
        biology <= biology_ceiling
    )

    pass_replicates = (
        replicate_count
        >= int(cfg.minimum_replicates)
    )

    pass_dominance = (
        dominance
        >= float(
            cfg.minimum_single_donor_dominance
        )
    )

    pass_lodo = (
        lodo_positive
        < float(
            cfg.lodo_positive_fraction_ceiling
        )
    )

    pass_dominance_or_lodo = (
        pass_dominance
        | pass_lodo
    )

    pass_direction = (
        direction
        < float(
            cfg.minimum_direction_agreement
        )
    )

    # Signed / scaled distances from each hard gate.
    m_technical = _signed_margin(
        technical - technical_threshold,
        pass_technical,
        scale=robust_scale(
            technical,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_risk = _signed_margin(
        risk - risk_threshold,
        pass_risk,
        scale=robust_scale(
            risk,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_biology = _signed_margin(
        biology_ceiling - biology,
        pass_biology,
        scale=robust_scale(
            biology,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_replicates = _signed_margin(
        replicate_count
        - float(
            cfg.minimum_replicates
        ),
        pass_replicates,
        scale=max(
            1.0,
            robust_scale(
                replicate_count,
                floor=scale_floor,
            ),
        ),
        epsilon=boundary_epsilon,
    )

    m_dominance = _signed_margin(
        dominance
        - float(
            cfg.minimum_single_donor_dominance
        ),
        pass_dominance,
        scale=robust_scale(
            dominance,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_lodo = _signed_margin(
        float(
            cfg.lodo_positive_fraction_ceiling
        )
        - lodo_positive,
        pass_lodo,
        scale=robust_scale(
            lodo_positive,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_or = _soft_or(
        [
            m_dominance,
            m_lodo,
        ]
    )

    m_direction = _signed_margin(
        float(
            cfg.minimum_direction_agreement
        )
        - direction,
        pass_direction,
        scale=robust_scale(
            direction,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    # Branch AND = weakest signed gate.
    branch_margin = _soft_and(
        [
            m_technical,
            m_risk,
            m_biology,
            m_replicates,
            m_or,
            m_direction,
        ]
    )

    hard_pre_protection = (
        pass_technical
        & pass_risk
        & pass_biology
        & pass_replicates
        & pass_dominance_or_lodo
        & pass_direction
    )

    protection_cols = [
        c
        for c in [
            "replicated_marker_protection",
            "rare_marker_protection",
            "user_protected",
        ]
        if c in t.columns
    ]

    protected = np.zeros(
        len(t),
        dtype=bool,
    )

    for c in protection_cols:
        protected |= strict_bool(
            t[c]
        )

    hard_remove = (
        hard_pre_protection
        & (~protected)
    )

    confidence_pre_protection = sigmoid(
        branch_margin
        / float(soft_temperature)
    )

    # Branch-specific protection only.
    confidence = np.where(
        protected,
        0.0,
        confidence_pre_protection,
    )

    # Strong internal invariant:
    # > 0.5 must exactly reproduce the historical donor decision.
    endpoint_call = confidence > 0.5

    if not np.array_equal(
        endpoint_call,
        hard_remove,
    ):
        bad = np.where(
            endpoint_call != hard_remove
        )[0]

        raise RuntimeError(
            "Donor endpoint-consistency invariant failed for "
            f"{len(bad)} genes. First indices: {bad[:10].tolist()}"
        )

    out = pd.DataFrame(
        {
            "gene": t["gene"].astype(str),
            "donor_technical_evidence": technical,
            "donor_biological_support": biology,
            "donor_net_risk": risk,
            "donor_technical_threshold": technical_threshold,
            "donor_risk_threshold": risk_threshold,
            "donor_biology_ceiling": biology_ceiling,
            "donor_pass_technical": pass_technical,
            "donor_pass_risk": pass_risk,
            "donor_pass_biology": pass_biology,
            "donor_pass_replicates": pass_replicates,
            "donor_pass_dominance": pass_dominance,
            "donor_pass_lodo": pass_lodo,
            "donor_pass_dominance_or_lodo": pass_dominance_or_lodo,
            "donor_pass_direction": pass_direction,
            "donor_margin_technical": m_technical,
            "donor_margin_risk": m_risk,
            "donor_margin_biology": m_biology,
            "donor_margin_replicates": m_replicates,
            "donor_margin_dominance": m_dominance,
            "donor_margin_lodo": m_lodo,
            "donor_margin_dominance_or_lodo": m_or,
            "donor_margin_direction": m_direction,
            "donor_branch_margin": branch_margin,
            "donor_confidence_pre_protection": confidence_pre_protection,
            "donor_protected": protected,
            "donor_confidence": confidence,
            "donor_hard_pre_protection": hard_pre_protection,
            "donor_hard_remove": hard_remove,
            "donor_endpoint_call_from_confidence": endpoint_call,
        }
    )

    # Post-hoc compatibility audit only; never used for selection.
    if "risk_gate_passed" in t.columns:
        packaged = strict_bool(
            t["risk_gate_passed"]
        )
        out[
            "donor_packaged_risk_gate_match"
        ] = (
            packaged
            == hard_pre_protection
        )

    if "final_action" in t.columns:
        packaged_remove = (
            t["final_action"]
            .astype(str)
            .eq("remove")
            .to_numpy()
        )

        out[
            "donor_packaged_final_call_match"
        ] = (
            packaged_remove
            == hard_remove
        )

    return out


# =====================================================================
# Cell-level evidence: hard gates -> standardized signed margin
# =====================================================================

def cell_soft_gate(
    evidence: pd.DataFrame,
    *,
    config: Optional[CellGateConfig] = None,
    soft_temperature: float = 1.0,
    boundary_epsilon: float = 1e-9,
    scale_floor: float = 1e-8,
) -> pd.DataFrame:
    cfg = config or CellGateConfig()

    if float(soft_temperature) <= 0:
        raise ValueError(
            "soft_temperature must be > 0"
        )

    t = evidence.copy().reset_index(
        drop=True
    )

    if "gene" not in t.columns:
        raise KeyError(
            "cell evidence lacks gene"
        )

    if (
        t["gene"].isna().any()
        or t["gene"].duplicated().any()
    ):
        raise ValueError(
            "cell evidence must contain one unique row per gene"
        )

    leakage_z = robust_z(
        _numeric(
            t,
            "donor_leakage",
        )
    )

    interaction_z = robust_z(
        _numeric(
            t,
            "interaction_range",
        )
    )

    biology_z = robust_z(
        _numeric(
            t,
            "biology_eta2",
        )
    )

    technical = (
        leakage_z
        + 0.75 * interaction_z
    )

    risk_raw = (
        technical
        - biology_z
    )

    risk_z = robust_z(
        risk_raw
    )

    fdr = _numeric(
        t,
        "permutation_fdr",
    )

    bootstrap = _numeric(
        t,
        "bootstrap_risk_fraction",
    )

    # Exact historical gate Booleans.
    pass_fdr = (
        fdr <= float(cfg.alpha)
    )

    pass_leakage = (
        leakage_z
        >= float(
            cfg.leakage_z_floor
        )
    )

    pass_risk = (
        risk_z
        >= float(
            cfg.risk_z_floor
        )
    )

    pass_biology = (
        biology_z
        <= float(
            cfg.biology_z_ceiling
        )
    )

    pass_bootstrap = (
        bootstrap
        >= float(
            cfg.bootstrap_pass_fraction
        )
    )

    # Signed / scaled distances.
    m_fdr = _signed_margin(
        float(cfg.alpha) - fdr,
        pass_fdr,
        scale=robust_scale(
            fdr,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_leakage = _signed_margin(
        leakage_z
        - float(
            cfg.leakage_z_floor
        ),
        pass_leakage,
        scale=robust_scale(
            leakage_z,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_risk = _signed_margin(
        risk_z
        - float(
            cfg.risk_z_floor
        ),
        pass_risk,
        scale=robust_scale(
            risk_z,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_biology = _signed_margin(
        float(
            cfg.biology_z_ceiling
        )
        - biology_z,
        pass_biology,
        scale=robust_scale(
            biology_z,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    m_bootstrap = _signed_margin(
        bootstrap
        - float(
            cfg.bootstrap_pass_fraction
        ),
        pass_bootstrap,
        scale=robust_scale(
            bootstrap,
            floor=scale_floor,
        ),
        epsilon=boundary_epsilon,
    )

    branch_margin = _soft_and(
        [
            m_fdr,
            m_leakage,
            m_risk,
            m_biology,
            m_bootstrap,
        ]
    )

    hard_pre_protection = (
        pass_fdr
        & pass_leakage
        & pass_risk
        & pass_biology
        & pass_bootstrap
    )

    protection_cols = [
        c
        for c in [
            "marker_protected",
            "explicitly_protected",
        ]
        if c in t.columns
    ]

    protected = np.zeros(
        len(t),
        dtype=bool,
    )

    for c in protection_cols:
        protected |= strict_bool(
            t[c]
        )

    hard_remove = (
        hard_pre_protection
        & (~protected)
    )

    confidence_pre_protection = sigmoid(
        branch_margin
        / float(soft_temperature)
    )

    # Branch-specific protection only.
    confidence = np.where(
        protected,
        0.0,
        confidence_pre_protection,
    )

    endpoint_call = confidence > 0.5

    if not np.array_equal(
        endpoint_call,
        hard_remove,
    ):
        bad = np.where(
            endpoint_call != hard_remove
        )[0]

        raise RuntimeError(
            "Cell endpoint-consistency invariant failed for "
            f"{len(bad)} genes. First indices: {bad[:10].tolist()}"
        )

    out = pd.DataFrame(
        {
            "gene": t["gene"].astype(str),
            "cell_donor_leakage_z": leakage_z,
            "cell_interaction_range_z": interaction_z,
            "cell_technical_evidence": technical,
            "cell_biological_support": biology_z,
            "cell_risk_score_raw": risk_raw,
            "cell_risk_score_z": risk_z,
            "cell_permutation_fdr": fdr,
            "cell_bootstrap_fraction": bootstrap,
            "cell_pass_fdr": pass_fdr,
            "cell_pass_leakage": pass_leakage,
            "cell_pass_risk": pass_risk,
            "cell_pass_biology": pass_biology,
            "cell_pass_bootstrap": pass_bootstrap,
            "cell_margin_fdr": m_fdr,
            "cell_margin_leakage": m_leakage,
            "cell_margin_risk": m_risk,
            "cell_margin_biology": m_biology,
            "cell_margin_bootstrap": m_bootstrap,
            "cell_branch_margin": branch_margin,
            "cell_confidence_pre_protection": confidence_pre_protection,
            "cell_protected": protected,
            "cell_confidence": confidence,
            "cell_hard_pre_protection": hard_pre_protection,
            "cell_hard_remove": hard_remove,
            "cell_endpoint_call_from_confidence": endpoint_call,
        }
    )

    # Post-hoc compatibility audit only.
    if "risk_flagged_before_protection" in t.columns:
        packaged = strict_bool(
            t["risk_flagged_before_protection"]
        )
        out[
            "cell_packaged_risk_gate_match"
        ] = (
            packaged
            == hard_pre_protection
        )

    if "final_action" in t.columns:
        packaged_remove = (
            t["final_action"]
            .astype(str)
            .eq("remove")
            .to_numpy()
        )

        out[
            "cell_packaged_final_call_match"
        ] = (
            packaged_remove
            == hard_remove
        )

    return out


