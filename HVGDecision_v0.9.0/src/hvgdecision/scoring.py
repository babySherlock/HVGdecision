"""scIB-inspired three-domain benchmark aggregation for HVGDecision papers.

This is an evaluation helper, not part of gene selection. It does not make
HVGDecision a machine-learning method and is not an exact reproduction of
scIB. Metrics are normalized panel-wise within each fixed benchmark task.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

DEFAULT_TRANSFER_METRICS = (
    "macro_f1",
    "balanced_accuracy",
    "rare_cell_macro_f1",
)
DEFAULT_BATCH_METRICS = (
    "dataset_mixing_score",
    "donor_mixing_within_celltype",
)
DEFAULT_BIOLOGY_METRICS = (
    "leiden_ari",
    "leiden_nmi",
    "celltype_silhouette_query",
)


def _minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite = values[np.isfinite(values)]
    if finite.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)
    low = float(finite.min())
    high = float(finite.max())
    if high - low <= 1e-12:
        # No panel differentiates on this metric, so it should not favor any panel.
        return pd.Series(0.5, index=series.index, dtype=float)
    return (values - low) / (high - low)


def three_domain_scores(
    frame: pd.DataFrame,
    *,
    task_keys: Sequence[str] = ("dataset", "integration_method"),
    panel_key: str = "panel",
    transfer_metrics: Sequence[str] = DEFAULT_TRANSFER_METRICS,
    batch_metrics: Sequence[str] = DEFAULT_BATCH_METRICS,
    biology_metrics: Sequence[str] = DEFAULT_BIOLOGY_METRICS,
    transfer_weight: float = 0.30,
    batch_weight: float = 0.40,
    biology_weight: float = 0.30,
) -> pd.DataFrame:
    """Calculate Transfer, Batch, Biological-fidelity and Overall scores.

    Each metric is min-max normalized across feature panels *inside a fixed
    task*, normally ``dataset × integration_method``. Missing metrics are
    omitted from that domain's mean. The default overall weighting preserves
    a 0.6 biological side (0.3 transfer + 0.3 biological fidelity) and a 0.4
    batch-correction side.
    """
    data = frame.copy()
    required = [*task_keys, panel_key]
    missing = [key for key in required if key not in data]
    if missing:
        raise KeyError(f"Missing required benchmark columns: {missing}")
    weights = np.asarray([transfer_weight, batch_weight, biology_weight], dtype=float)
    if np.any(weights < 0) or not np.isclose(weights.sum(), 1.0):
        raise ValueError("transfer_weight + batch_weight + biology_weight must equal 1")

    domains: dict[str, Sequence[str]] = {
        "transfer_score": tuple(transfer_metrics),
        "batch_correction_score": tuple(batch_metrics),
        "biological_fidelity_score": tuple(biology_metrics),
    }
    available = {
        domain: [metric for metric in metrics if metric in data.columns]
        for domain, metrics in domains.items()
    }
    if not available["transfer_score"]:
        raise KeyError("No requested Transfer metrics were found")
    if not available["batch_correction_score"]:
        raise KeyError("No requested Batch-correction metrics were found")
    if not available["biological_fidelity_score"]:
        raise KeyError("No requested Biological-fidelity metrics were found")

    normalized_parts = []
    group_keys: Any = list(task_keys) if len(task_keys) > 1 else task_keys[0]
    for _, group in data.groupby(group_keys, sort=False, dropna=False):
        part = group.copy()
        normalized_columns: dict[str, list[str]] = {key: [] for key in domains}
        for domain, metrics in available.items():
            for metric in metrics:
                column = f"__norm__{metric}"
                part[column] = _minmax(part[metric])
                normalized_columns[domain].append(column)
            part[domain] = part[normalized_columns[domain]].mean(axis=1, skipna=True)
        part["overall_score"] = (
            transfer_weight * part["transfer_score"]
            + batch_weight * part["batch_correction_score"]
            + biology_weight * part["biological_fidelity_score"]
        )
        normalized_parts.append(part)

    result = pd.concat(normalized_parts, ignore_index=True)
    keep = [
        *required,
        "transfer_score",
        "batch_correction_score",
        "biological_fidelity_score",
        "overall_score",
    ]
    return result[keep].copy()
