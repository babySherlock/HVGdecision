"""Statistical primitives used by HVGDecision.

The public 0.9 method has two feature-refinement modes:

* ``within_domain``: Reference-only donor/batch risk conditional on biology.
* ``cross_domain``: Reference risk is combined with label-free Reference/Query shift.

This module intentionally contains only deterministic statistical primitives;
no machine-learning model is fitted for gene selection.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np


def robust_z(values: np.ndarray) -> np.ndarray:
    """Median/MAD z score with an SD fallback for degenerate MADs."""
    values = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(values)
    fill = float(np.nanmedian(values[finite])) if finite.any() else 0.0
    values = np.where(finite, values, fill)
    center = float(np.median(values))
    mad = float(np.median(np.abs(values - center)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 1e-8:
        sd = float(np.std(values, ddof=0))
        scale = sd if np.isfinite(sd) and sd > 1e-8 else 1.0
    return (values - center) / scale


def eta_squared(matrix: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Per-gene fraction of total variance explained by categorical labels."""
    matrix = np.asarray(matrix)
    grand = matrix.mean(axis=0, dtype=np.float64)
    total = np.square(matrix.astype(np.float64) - grand).sum(axis=0)
    between = np.zeros(matrix.shape[1], dtype=np.float64)
    for label in sorted(np.unique(labels)):
        positions = np.flatnonzero(labels == label)
        if positions.size:
            mean = matrix[positions].mean(axis=0, dtype=np.float64)
            between += positions.size * np.square(mean - grand)
    return np.divide(between, total, out=np.zeros_like(between), where=total > 1e-12)


def benjamini_hochberg(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def donor_risk_components(
    matrix: np.ndarray,
    donor_labels: np.ndarray,
    celltype_labels: np.ndarray,
    min_group_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return cell-type-conditional donor leakage and interaction instability.

    Leakage is donor eta-squared after subtracting the corresponding biology-
    stratum mean from every cell. Interaction instability is the median,
    across eligible biology strata, of the donor mean range divided by the
    within-stratum SD.
    """
    matrix = np.asarray(matrix, dtype=np.float32)
    residual = matrix.copy()
    for celltype in sorted(np.unique(celltype_labels)):
        positions = np.flatnonzero(celltype_labels == celltype)
        residual[positions] -= matrix[positions].mean(axis=0, keepdims=True)
    leakage = eta_squared(residual, donor_labels)

    instability_rows = []
    for celltype in sorted(np.unique(celltype_labels)):
        cell_positions = np.flatnonzero(celltype_labels == celltype)
        if cell_positions.size == 0:
            continue
        denominator = matrix[cell_positions].std(axis=0) + 1e-6
        donor_means = []
        for donor in sorted(np.unique(donor_labels)):
            positions = np.flatnonzero((celltype_labels == celltype) & (donor_labels == donor))
            if positions.size >= min_group_cells:
                donor_means.append(matrix[positions].mean(axis=0))
        if len(donor_means) >= 2:
            instability_rows.append(np.ptp(np.vstack(donor_means), axis=0) / denominator)
    instability = (
        np.nanmedian(np.vstack(instability_rows), axis=0)
        if instability_rows
        else np.zeros(matrix.shape[1], dtype=np.float64)
    )
    return leakage, instability


def biology_components(
    matrix: np.ndarray,
    donor_labels: np.ndarray,
    celltype_labels: np.ndarray,
    min_group_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return biological explained variance and an audit-only consistency score.

    ``biology`` is the direct cell-type/biology eta-squared used by the final
    within-domain rule. ``consistency`` is retained for audit compatibility but
    is *not* part of the 0.9 risk equation.
    """
    matrix = np.asarray(matrix, dtype=np.float32)
    biology = eta_squared(matrix, celltype_labels)

    donors = sorted(np.unique(donor_labels))
    celltypes = sorted(np.unique(celltype_labels))
    group_means: dict[tuple[str, str], np.ndarray] = {}
    for donor in donors:
        for celltype in celltypes:
            positions = np.flatnonzero((donor_labels == donor) & (celltype_labels == celltype))
            if positions.size >= min_group_cells:
                group_means[(donor, celltype)] = matrix[positions].mean(axis=0)
    correlations = []
    for donor_a, donor_b in combinations(donors, 2):
        common = [
            celltype
            for celltype in celltypes
            if (donor_a, celltype) in group_means and (donor_b, celltype) in group_means
        ]
        if len(common) < 3:
            continue
        a = np.vstack([group_means[(donor_a, celltype)] for celltype in common]).astype(float)
        b = np.vstack([group_means[(donor_b, celltype)] for celltype in common]).astype(float)
        a -= a.mean(axis=0, keepdims=True)
        b -= b.mean(axis=0, keepdims=True)
        denominator = np.sqrt(np.square(a).sum(axis=0) * np.square(b).sum(axis=0))
        correlations.append(
            np.divide(
                (a * b).sum(axis=0),
                denominator,
                out=np.full(matrix.shape[1], np.nan),
                where=denominator > 1e-12,
            )
        )
    consistency = (
        np.nanmedian(np.vstack(correlations), axis=0)
        if correlations
        else np.zeros(matrix.shape[1], dtype=np.float64)
    )
    return biology, consistency


def compose_risk(
    leakage: np.ndarray,
    instability: np.ndarray,
    biology: np.ndarray,
    consistency: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Final within-domain raw risk.

    R_g = Z(leakage_g) + 0.75 Z(interaction_g) - Z(biology_g)

    The returned ``risk_z`` is the robust z score of this raw risk and is used
    by the final gate. ``consistency`` is accepted only for backward API
    compatibility and never changes the score.
    """
    leakage_z = robust_z(leakage)
    interaction_z = robust_z(instability)
    biology_z = robust_z(biology)
    raw = leakage_z + 0.75 * interaction_z - biology_z
    risk_z = robust_z(raw)
    components = {
        "donor_leakage_z": leakage_z,
        "interaction_instability_z": interaction_z,
        "interaction_range_z": interaction_z,  # old notebook alias
        "protected_biology_z": biology_z,
        "celltype_biology_z": biology_z,
        "biology_z": biology_z,
        "risk_z": risk_z,
    }
    if consistency is not None:
        components["cross_donor_consistency_z"] = robust_z(consistency)
    return raw, components


def holm_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = np.maximum.accumulate((len(values) - np.arange(len(values))) * ranked)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted
