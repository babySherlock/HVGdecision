"""Scanpy-style in-memory API for Reference/Query HVG experiments."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import balanced_accuracy_score, f1_score, recall_score
from tqdm.auto import tqdm

from .inspect import inspect_dataset
from .io import (
    audit_counts,
    count_source_label,
    count_var_names,
    external_count_matrix,
    looks_like_count_table_path,
    resolve_counts_source,
    role_masks,
    write_csv,
    write_json,
)
from .refine import _hvg_ranking, refine_panel
from .modes import normalize_mode

DEFAULT_HVG_BUDGETS = (
    500,
    750,
    1000,
    1250,
    1500,
    1750,
    1800,
    1900,
    2000,
    2250,
    2500,
    2750,
    3000,
    3500,
    4000,
    4500,
    5000,
)

_UNSET = object()

def _normalized_method_name(value: str) -> str:
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")


def _is_seurat_v3_like(value: str) -> bool:
    name = _normalized_method_name(value)
    return name in {"seurat_v3", "seuratv3", "scanpy_seurat_v3", "scanpy_seuratv3"} or (
        "seurat" in name and "v3" in name
    )


def _budget_composition_audit(ranking_tables: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Compare every budget with the largest-budget prefix.

    A strict truncation means the independently obtained N-gene panel is
    exactly the first N genes of the largest panel, in the same order.
    """
    if not ranking_tables:
        return pd.DataFrame()
    largest = max(ranking_tables)
    base = (
        ranking_tables[largest]
        .sort_values("gene_rank", kind="mergesort")["gene"]
        .astype(str)
        .tolist()
    )
    base_set = set(base)
    rows = []
    for budget in sorted(ranking_tables):
        genes = (
            ranking_tables[budget]
            .sort_values("gene_rank", kind="mergesort")["gene"]
            .astype(str)
            .tolist()
        )
        prefix = base[: len(genes)]
        prefix_set = set(prefix)
        gene_set = set(genes)
        rows.append(
            {
                "largest_n_hvg": int(largest),
                "n_hvg": int(budget),
                "observed_n_genes": int(len(genes)),
                "declared_budget_matches_size": bool(len(genes) == int(budget)),
                "overlap_with_largest": int(len(gene_set & base_set)),
                "overlap_with_largest_prefix": int(len(gene_set & prefix_set)),
                "new_vs_largest_prefix": int(len(gene_set - prefix_set)),
                "missing_from_largest_prefix": int(len(prefix_set - gene_set)),
                "same_set_as_largest_prefix": bool(gene_set == prefix_set),
                "strict_prefix_truncation": bool(genes == prefix),
            }
        )
    return pd.DataFrame(rows)


def _read_external_hvg_input(
    hvg: Any,
    *,
    gene_column: str,
    rank_column: str,
    budget_column: str,
    method_column: str,
    method_name: str,
) -> pd.DataFrame:
    """Normalize an external R/Python HVG table while preserving order."""
    if isinstance(hvg, pd.DataFrame):
        frame = hvg.copy()
    elif isinstance(hvg, pd.Series):
        frame = pd.DataFrame({gene_column: hvg.astype(str).tolist()})
    elif isinstance(hvg, (str, Path)):
        path = Path(hvg).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        suffixes = {suffix.lower() for suffix in path.suffixes}
        sep = "\t" if {".tsv", ".txt"} & suffixes else ","
        frame = pd.read_csv(path, sep=sep)
    elif isinstance(hvg, Sequence) and not isinstance(hvg, (str, bytes)):
        frame = pd.DataFrame({gene_column: [str(value) for value in hvg]})
    else:
        raise TypeError(
            "hvg must be a gene sequence, pandas Series/DataFrame, or CSV/TSV path"
        )
    if gene_column not in frame:
        if frame.shape[1] == 1:
            frame = frame.rename(columns={frame.columns[0]: gene_column})
        else:
            raise KeyError(f"External HVG table is missing gene column {gene_column!r}")
    frame = frame.copy()
    frame[gene_column] = frame[gene_column].astype(str)
    frame = frame.loc[frame[gene_column].ne("")].copy()
    if rank_column not in frame:
        if budget_column in frame:
            frame[rank_column] = frame.groupby(budget_column, sort=False).cumcount() + 1
        else:
            frame[rank_column] = np.arange(1, len(frame) + 1)
    frame[rank_column] = pd.to_numeric(frame[rank_column], errors="raise").astype(int)
    if budget_column not in frame:
        frame[budget_column] = int(frame[rank_column].max()) if len(frame) else 0
    frame[budget_column] = pd.to_numeric(frame[budget_column], errors="raise").astype(int)
    if method_column not in frame:
        frame[method_column] = str(method_name)
    frame[method_column] = frame[method_column].astype(str)
    frame = frame.rename(
        columns={
            gene_column: "gene",
            rank_column: "gene_rank",
            budget_column: "n_hvg",
            method_column: "method",
        }
    )
    frame["_input_order"] = np.arange(len(frame), dtype=int)
    frame = frame.sort_values(["n_hvg", "gene_rank", "_input_order"], kind="mergesort")
    duplicated = frame.duplicated(["method", "n_hvg", "gene"], keep="first")
    frame = frame.loc[~duplicated].copy()
    return frame[["method", "n_hvg", "gene_rank", "gene", "_input_order"]].reset_index(drop=True)



def _coerce_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return (
        series.fillna(False)
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )

def _values(value: str | Sequence[str], name: str) -> list[str]:
    if isinstance(value, str):
        result = [value]
    else:
        result = [str(item) for item in value]
    result = list(dict.fromkeys(result))
    if not result:
        raise ValueError(f"{name} cannot be empty")
    return result



def _weighted_panel_recovery(base_genes: list[str], comparison_genes: list[str]) -> tuple[int, float, float, float]:
    """Return overlap, overlap fraction, Jaccard, and rank-weighted recovery.

    The weighted score gives more weight to genes near the top of the full-
    Reference ranking.  It is used only for Reference-only HVG-budget
    stability; no integration model or Query expression is involved.
    """
    base = list(map(str, base_genes))
    comparison = set(map(str, comparison_genes))
    base_set = set(base)
    overlap = int(len(base_set & comparison))
    denominator = max(len(base), 1)
    overlap_fraction = float(overlap / denominator)
    union = len(base_set | comparison)
    jaccard = float(overlap / union) if union else 1.0
    ranks = np.arange(1, len(base) + 1, dtype=float)
    weights = 1.0 / np.log2(ranks + 1.0)
    hit = np.fromiter((gene in comparison for gene in base), dtype=float, count=len(base))
    weighted = float(np.sum(weights * hit) / np.sum(weights)) if len(weights) else 1.0
    return overlap, overlap_fraction, jaccard, weighted


def _balanced_reference_sample_indices(
    adata,
    reference: np.ndarray,
    config: dict[str, Any],
    *,
    max_cells_per_batch_label: int,
) -> np.ndarray:
    """Deterministically cap Reference cells per batch × label for budget evaluation."""
    batch_key = config["batch_key"]
    label_key = config["label_key"]
    rng = np.random.default_rng(int(config.get("selection_seed", 20260828)))
    reference_positions = np.flatnonzero(reference)
    meta = adata.obs.iloc[reference_positions][[batch_key, label_key]].copy()
    selected: list[int] = []
    for _, group in meta.groupby([batch_key, label_key], observed=True, sort=True):
        local = meta.index.get_indexer(group.index)
        available = reference_positions[local]
        take = min(len(available), int(max_cells_per_batch_label))
        chosen = rng.choice(available, size=take, replace=False)
        selected.extend(np.sort(chosen).tolist())
    return np.asarray(sorted(set(selected)), dtype=int)


def _log_normalized_subset(
    adata,
    obs_positions: np.ndarray,
    genes: list[str],
    config: dict[str, Any],
):
    """Return sparse log1p(1e4-normalized counts) for selected cells/genes."""
    from .io import subset_count_matrix

    counts = subset_count_matrix(
        adata,
        np.asarray(obs_positions, dtype=int),
        list(map(str, genes)),
        config.get("counts_layer"),
    )
    matrix = (
        counts.tocsr().astype(np.float32)
        if sparse.issparse(counts)
        else sparse.csr_matrix(counts, dtype=np.float32)
    )
    library = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.divide(
        1e4,
        library,
        out=np.zeros_like(library, dtype=np.float32),
        where=library > 0,
    )
    logged = (sparse.diags(scale) @ matrix).tocsr()
    logged.data = np.log1p(logged.data)
    return logged


def _l2_normalize_sparse_rows(matrix):
    matrix = matrix.tocsr().astype(np.float32, copy=True)
    norms = np.sqrt(np.asarray(matrix.multiply(matrix).sum(axis=1)).ravel())
    inverse = np.divide(
        1.0,
        norms,
        out=np.zeros_like(norms, dtype=np.float32),
        where=norms > 1e-12,
    )
    return sparse.diags(inverse) @ matrix


def _reference_biology_transfer(
    adata,
    sampled_positions: np.ndarray,
    genes: list[str],
    config: dict[str, Any],
    reference_batches: list[str],
) -> pd.DataFrame:
    """Reference-only donor-held-out nearest-centroid cell-type transfer.

    v0.9 adds two safeguards against tiny panels that classify abundant cell
    types well but lose weaker biology: a training-defined rare-cell macro-F1
    and the minimum per-class recall in the held-out Reference donor.
    """
    label_key = config["label_key"]
    batch_key = config["batch_key"]
    rare_quantile = float(config.get("budget_rare_reference_quantile", 0.25))
    matrix = _log_normalized_subset(adata, sampled_positions, genes, config)
    obs = adata.obs.iloc[sampled_positions]
    donors = obs[batch_key].astype(str).to_numpy()
    labels = obs[label_key].astype(str).to_numpy()

    rows: list[dict[str, Any]] = []
    for held_out in reference_batches:
        train = donors != str(held_out)
        test = donors == str(held_out)
        train_labels = sorted(np.unique(labels[train]).tolist())
        if not train_labels or not np.any(test):
            continue

        train_counts = pd.Series(labels[train], dtype="object").value_counts()
        rare_cutoff = (
            float(train_counts.quantile(rare_quantile))
            if len(train_counts)
            else np.nan
        )
        rare_labels = set(
            train_counts[train_counts <= rare_cutoff].index.astype(str).tolist()
        ) if pd.notna(rare_cutoff) else set()

        centroids = []
        centroid_labels = []
        for label in train_labels:
            positions = np.flatnonzero(train & (labels == label))
            if positions.size == 0:
                continue
            centroid = np.asarray(matrix[positions].mean(axis=0)).ravel().astype(np.float32)
            centroids.append(centroid)
            centroid_labels.append(label)
        if not centroids:
            continue

        centroid_matrix = np.vstack(centroids).astype(np.float32)
        centroid_norm = np.linalg.norm(centroid_matrix, axis=1)
        centroid_matrix = np.divide(
            centroid_matrix,
            centroid_norm[:, None],
            out=np.zeros_like(centroid_matrix),
            where=centroid_norm[:, None] > 1e-12,
        )

        test_positions = np.flatnonzero(test)
        y_all = labels[test_positions]
        known = np.isin(y_all, np.asarray(centroid_labels, dtype=object))
        if not np.any(known):
            rows.append(
                {
                    "held_out_reference_batch": str(held_out),
                    "biology_eval_cells_total": int(len(test_positions)),
                    "biology_eval_cells_known": 0,
                    "biology_label_coverage": 0.0,
                    "biology_macro_f1": np.nan,
                    "biology_balanced_accuracy": np.nan,
                    "biology_rare_macro_f1": np.nan,
                    "biology_rare_eval_cells": 0,
                    "biology_rare_label_count": 0,
                    "biology_min_class_recall": np.nan,
                    "biology_transfer_score": np.nan,
                    "biology_sufficiency_score": np.nan,
                }
            )
            continue

        eval_positions = test_positions[known]
        y_true = labels[eval_positions]
        x_test = _l2_normalize_sparse_rows(matrix[eval_positions])
        similarity = np.asarray(x_test @ centroid_matrix.T)
        y_pred = np.asarray(centroid_labels, dtype=object)[np.argmax(similarity, axis=1)]

        macro_f1 = float(f1_score(y_true, y_pred, average="macro", zero_division=0))
        balanced = float(balanced_accuracy_score(y_true, y_pred))

        eval_classes = sorted(np.unique(y_true).tolist())
        recalls = recall_score(
            y_true,
            y_pred,
            labels=eval_classes,
            average=None,
            zero_division=0,
        )
        min_class_recall = float(np.min(recalls)) if len(recalls) else np.nan

        rare_eval_labels = sorted(set(eval_classes) & rare_labels)
        rare_mask = np.isin(y_true, np.asarray(rare_eval_labels, dtype=object))
        if rare_eval_labels and np.any(rare_mask):
            rare_macro_f1 = float(
                f1_score(
                    y_true[rare_mask],
                    y_pred[rare_mask],
                    labels=rare_eval_labels,
                    average="macro",
                    zero_division=0,
                )
            )
            rare_eval_cells = int(np.sum(rare_mask))
        else:
            rare_macro_f1 = np.nan
            rare_eval_cells = 0

        coverage = float(len(eval_positions) / max(len(test_positions), 1))
        transfer_score = float(0.5 * macro_f1 + 0.5 * balanced)

        sufficiency_components = [macro_f1, balanced, min_class_recall]
        if np.isfinite(rare_macro_f1):
            sufficiency_components.append(rare_macro_f1)
        biology_sufficiency_score = float(np.mean(sufficiency_components))

        rows.append(
            {
                "held_out_reference_batch": str(held_out),
                "biology_eval_cells_total": int(len(test_positions)),
                "biology_eval_cells_known": int(len(eval_positions)),
                "biology_label_coverage": coverage,
                "biology_macro_f1": macro_f1,
                "biology_balanced_accuracy": balanced,
                "biology_rare_macro_f1": rare_macro_f1,
                "biology_rare_eval_cells": rare_eval_cells,
                "biology_rare_label_count": int(len(rare_eval_labels)),
                "biology_min_class_recall": min_class_recall,
                "biology_transfer_score": transfer_score,
                "biology_sufficiency_score": biology_sufficiency_score,
            }
        )
    return pd.DataFrame(rows)


def _composite_budget_biology_labels(obs: pd.DataFrame, biology_keys: list[str]) -> np.ndarray:
    if len(biology_keys) == 1:
        return obs[biology_keys[0]].astype(str).to_numpy()
    return (
        obs[biology_keys]
        .astype(str)
        .agg(" | ".join, axis=1)
        .to_numpy(dtype=object)
    )


def _within_biology_donor_instability(
    adata,
    sampled_positions: np.ndarray,
    genes: list[str],
    config: dict[str, Any],
) -> float:
    """Mean within-biology donor centroid cosine distance; lower is better."""
    batch_key = config["batch_key"]
    biology_keys = list(dict.fromkeys(config.get("biology_keys", [config["label_key"]])))
    min_cells = int(config.get("budget_min_cells_per_donor_biology", 5))
    matrix = _log_normalized_subset(adata, sampled_positions, genes, config)
    obs = adata.obs.iloc[sampled_positions].copy()
    obs["_budget_biology_stratum"] = _composite_budget_biology_labels(obs, biology_keys)
    donors = obs[batch_key].astype(str).to_numpy()
    strata = obs["_budget_biology_stratum"].astype(str).to_numpy()

    stratum_scores: list[float] = []
    for stratum in sorted(np.unique(strata)):
        stratum_mask = strata == stratum
        donor_centroids = []
        for donor in sorted(np.unique(donors[stratum_mask])):
            positions = np.flatnonzero(stratum_mask & (donors == donor))
            if positions.size < min_cells:
                continue
            centroid = np.asarray(matrix[positions].mean(axis=0)).ravel().astype(np.float32)
            norm = float(np.linalg.norm(centroid))
            if norm <= 1e-12:
                continue
            donor_centroids.append(centroid / norm)

        if len(donor_centroids) < 2:
            continue
        centroid_matrix = np.vstack(donor_centroids)
        similarity = np.clip(centroid_matrix @ centroid_matrix.T, -1.0, 1.0)
        upper = np.triu_indices(len(donor_centroids), k=1)
        distances = 1.0 - similarity[upper]
        if distances.size:
            stratum_scores.append(float(np.mean(distances)))

    return float(np.mean(stratum_scores)) if stratum_scores else np.nan


def _robust_location_scale(values: pd.Series) -> tuple[float, float]:
    finite = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return np.nan, np.nan
    median = float(finite.median())
    mad = float(np.median(np.abs(finite.to_numpy(dtype=float) - median)))
    robust_sigma = float(1.4826 * mad)
    return median, robust_sigma


def _adaptive_plateau_tolerance(
    frame: pd.DataFrame,
    *,
    std_column: str,
    n_column: str,
    minimum: float,
    maximum: float,
) -> float:
    """Noise-adaptive absolute plateau width, bounded conservatively."""
    if std_column in frame and n_column in frame:
        std = pd.to_numeric(frame[std_column], errors="coerce").to_numpy(dtype=float)
        n = pd.to_numeric(frame[n_column], errors="coerce").to_numpy(dtype=float)
        sem = np.divide(
            std,
            np.sqrt(n),
            out=np.full_like(std, np.nan, dtype=float),
            where=np.isfinite(std) & np.isfinite(n) & (n > 0),
        )
        finite = sem[np.isfinite(sem)]
        if finite.size:
            return float(np.clip(np.median(finite), minimum, maximum))
    return float(minimum)


def _recommend_joint_reference_budget(
    summary: pd.DataFrame,
    *,
    stability_tolerance: float,
    biology_tolerance: float,
    donor_guardrail_quantile: float,
    biology_plateau_min_tolerance: float = 0.01,
    biology_plateau_max_tolerance: float = 0.015,
    rare_biology_plateau_min_tolerance: float = 0.015,
    rare_biology_plateau_max_tolerance: float = 0.03,
    stability_mad_multiplier: float = 2.5,
    donor_instability_mad_multiplier: float = 2.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """v0.9 Reference-only minimum-sufficient budget rule.

    Selection priority is intentionally asymmetric:
      1. biology plateau is the primary sufficiency gate;
      2. within-biology donor instability excludes clear batch-heavy panels;
      3. HVG stability is only a low-quality/outlier guardrail, not a demand
         to sit near the global stability maximum.

    This prevents the U-shaped stability curve from automatically favouring
    very small or very large panels.
    """
    frame = summary.copy().sort_values("n_hvg").reset_index(drop=True)
    if frame.empty:
        raise RuntimeError("No budget metrics are available")

    # ------------------------------------------------------------
    # A. Stability = quality floor, NOT near-peak optimization.
    # ------------------------------------------------------------
    stability = pd.to_numeric(frame["stability_score_mean"], errors="coerce")
    if stability.notna().sum() == 0:
        raise RuntimeError("Reference-only HVG stability could not be estimated")
    peak_stability = float(stability.max())
    stability_median, stability_sigma = _robust_location_scale(stability)
    robust_stability_floor = (
        stability_median - float(stability_mad_multiplier) * stability_sigma
        if np.isfinite(stability_sigma)
        else np.nan
    )
    legacy_near_peak_floor = peak_stability - float(stability_tolerance)
    if np.isfinite(robust_stability_floor):
        stability_threshold = float(min(legacy_near_peak_floor, robust_stability_floor))
    else:
        stability_threshold = float(legacy_near_peak_floor)
    stability_threshold = max(0.0, stability_threshold)
    frame["stability_sufficient"] = stability >= stability_threshold
    frame["stability_quality_floor"] = stability_threshold

    # ------------------------------------------------------------
    # B. Biology = primary plateau. v0.9 uses the richer score when
    #    available: overall F1 + balanced accuracy + rare-cell F1 +
    #    minimum class recall. Old transfer score remains compatible.
    # ------------------------------------------------------------
    biology_column = (
        "biology_sufficiency_score_mean"
        if "biology_sufficiency_score_mean" in frame
        else "biology_transfer_score_mean"
    )
    biology = pd.to_numeric(frame[biology_column], errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    )
    if biology.notna().sum() == 0:
        raise RuntimeError(
            "Reference-only biology sufficiency could not be estimated. "
            "Ensure Reference donors share at least one label_key cell type."
        )

    peak_biology = float(biology.max())
    biology_peak_idx = int(biology.idxmax())

    adaptive_biology_tol = _adaptive_plateau_tolerance(
        frame,
        std_column=(
            "biology_sufficiency_score_std"
            if "biology_sufficiency_score_std" in frame
            else "biology_transfer_score_std"
        ),
        n_column="biology_n_reference_holdouts",
        minimum=float(biology_plateau_min_tolerance),
        maximum=float(biology_plateau_max_tolerance),
    )
    # Legacy biology_tolerance may make the rule stricter, but v0.9 will not
    # let the old 0.03 default widen the plateau beyond the new conservative cap.
    adaptive_biology_tol = float(min(adaptive_biology_tol, float(biology_tolerance)))
    biology_threshold = peak_biology - adaptive_biology_tol
    frame["biology_plateau_pass"] = biology >= biology_threshold

    # Rare-cell gate: only active when the dataset provides estimable rare labels.
    rare_available = False
    rare_threshold = np.nan
    rare_tol = np.nan
    if "biology_rare_macro_f1_mean" in frame:
        rare = pd.to_numeric(
            frame["biology_rare_macro_f1_mean"], errors="coerce"
        ).replace([np.inf, -np.inf], np.nan)
        if rare.notna().any():
            rare_available = True
            peak_rare = float(rare.max())
            rare_tol = _adaptive_plateau_tolerance(
                frame,
                std_column="biology_rare_macro_f1_std",
                n_column="biology_rare_n_reference_holdouts",
                minimum=float(rare_biology_plateau_min_tolerance),
                maximum=float(rare_biology_plateau_max_tolerance),
            )
            rare_threshold = peak_rare - rare_tol
            frame["rare_biology_sufficient"] = rare.notna() & (rare >= rare_threshold)
        else:
            frame["rare_biology_sufficient"] = True
    else:
        frame["rare_biology_sufficient"] = True

    frame["biology_sufficient"] = (
        frame["biology_plateau_pass"]
        & frame["rare_biology_sufficient"]
    )

    # ------------------------------------------------------------
    # C. Donor instability = robust upper-outlier guardrail.
    #    q75 is retained only as a compatibility safeguard; the threshold
    #    is never made stricter than q75.
    # ------------------------------------------------------------
    instability = pd.to_numeric(
        frame["within_biology_donor_instability_mean"], errors="coerce"
    ).replace([np.inf, -np.inf], np.nan)
    if instability.notna().any():
        donor_median, donor_sigma = _robust_location_scale(instability)
        robust_donor_threshold = (
            donor_median + float(donor_instability_mad_multiplier) * donor_sigma
            if np.isfinite(donor_sigma)
            else np.nan
        )
        quantile_threshold = float(
            instability.quantile(float(donor_guardrail_quantile))
        )
        if np.isfinite(robust_donor_threshold):
            donor_threshold = float(max(robust_donor_threshold, quantile_threshold))
        else:
            donor_threshold = quantile_threshold
        frame["donor_guardrail_pass"] = (
            instability.isna() | (instability <= donor_threshold)
        )
        donor_guardrail_available = True
    else:
        donor_threshold = np.nan
        frame["donor_guardrail_pass"] = True
        donor_guardrail_available = False

    frame["joint_sufficient"] = (
        frame["biology_sufficient"]
        & frame["donor_guardrail_pass"]
        & frame["stability_sufficient"]
    )

    eligible = frame.loc[frame["joint_sufficient"]].copy()
    selection_status = "biology_plateau_quality_pass"
    if not eligible.empty:
        best_n = int(eligible["n_hvg"].min())
    else:
        # Biology remains primary in fallback. First preserve quality guardrails,
        # then choose the best biology score and the smallest tied budget.
        quality = frame.loc[
            frame["donor_guardrail_pass"] & frame["stability_sufficient"]
        ].copy()
        if quality.empty:
            quality = frame.copy()
            selection_status = "biology_best_global_fallback"
        else:
            selection_status = "biology_best_among_quality_pass_fallback"

        best_biology = float(pd.to_numeric(quality[biology_column], errors="coerce").max())
        candidates = quality.loc[
            pd.to_numeric(quality[biology_column], errors="coerce")
            >= best_biology - adaptive_biology_tol
        ].copy()
        if rare_available:
            rare_candidates = candidates.loc[candidates["rare_biology_sufficient"]]
            if not rare_candidates.empty:
                candidates = rare_candidates
        best_n = int(candidates["n_hvg"].min())

    recommended = frame.loc[frame["n_hvg"].eq(best_n)].iloc[0]
    peak_stability_row = frame.loc[stability.idxmax()]
    peak_biology_row = frame.loc[biology_peak_idx]

    recommendation = pd.DataFrame(
        [
            {
                "recommended_minimum_sufficient_hvg": int(best_n),
                "observed_peak_hvg": int(peak_stability_row["n_hvg"]),
                "observed_peak_stability_hvg": int(peak_stability_row["n_hvg"]),
                "observed_peak_biology_hvg": int(peak_biology_row["n_hvg"]),
                "recommended_stability_score": float(recommended["stability_score_mean"]),
                "observed_peak_stability_score": peak_stability,
                "recommended_biology_transfer_score": float(
                    recommended.get("biology_transfer_score_mean", recommended[biology_column])
                ),
                "recommended_biology_sufficiency_score": float(
                    recommended[biology_column]
                ),
                "observed_peak_biology_transfer_score": float(
                    frame.get("biology_transfer_score_mean", biology).max()
                ),
                "observed_peak_biology_sufficiency_score": peak_biology,
                "recommended_within_biology_donor_instability": float(
                    recommended["within_biology_donor_instability_mean"]
                )
                if pd.notna(recommended["within_biology_donor_instability_mean"])
                else np.nan,
                "stability_quality_floor": float(stability_threshold),
                "stability_robust_median": float(stability_median),
                "stability_robust_sigma": float(stability_sigma),
                "legacy_stability_near_optimal_tolerance": float(stability_tolerance),
                "biology_plateau_tolerance_used": float(adaptive_biology_tol),
                "biology_plateau_threshold": float(biology_threshold),
                "rare_biology_gate_available": bool(rare_available),
                "rare_biology_plateau_tolerance_used": float(rare_tol)
                if np.isfinite(rare_tol) else np.nan,
                "rare_biology_plateau_threshold": float(rare_threshold)
                if np.isfinite(rare_threshold) else np.nan,
                "legacy_biology_near_optimal_tolerance": float(biology_tolerance),
                "donor_instability_guardrail_quantile": float(donor_guardrail_quantile),
                "donor_instability_guardrail_threshold": donor_threshold,
                "donor_guardrail_available": bool(donor_guardrail_available),
                "selection_status": selection_status,
                "selection_basis": (
                    "reference_biology_plateau_primary_rare_cell_guard_"
                    "donor_instability_outlier_guard_stability_quality_floor"
                ),
                "query_expression_used": False,
                "query_labels_used": False,
            }
        ]
    )
    return frame, recommendation


def _make_final_adata(
    adata,
    config: dict[str, Any],
    reference: np.ndarray,
    query: np.ndarray,
    decision: pd.DataFrame,
    candidate_genes: list[str],
    final_genes: list[str],
    audit: dict[str, Any],
):
    """Return all selected cells and all count-source genes with HVG annotations.

    ``highly_variable`` marks the final HVG panel. Raw counts are always copied
    to ``layers['counts']`` so the returned object can be passed to count-based
    downstream methods such as scVI/scANVI without re-running HVG selection.
    """
    import anndata as ad

    from .io import subset_count_matrix

    keep = reference | query
    count_genes = count_var_names(adata, config.get("counts_layer"))
    count_genes = pd.Index(count_genes.astype(str), name="gene")
    current = pd.Index(adata.var_names.astype(str))

    if count_genes.isin(current).all():
        final = adata[keep, count_genes].copy()
    else:
        counts_all = subset_count_matrix(
            adata,
            keep,
            count_genes.astype(str).tolist(),
            config.get("counts_layer"),
        ).copy()
        original_var = adata.var.copy()
        original_var.index = current
        var = original_var.reindex(count_genes).copy()
        final = ad.AnnData(X=counts_all.copy(), obs=adata.obs.loc[keep].copy(), var=var)

    counts = subset_count_matrix(
        adata,
        keep,
        count_genes.astype(str).tolist(),
        config.get("counts_layer"),
    ).copy()
    final.layers["counts"] = counts

    genes = pd.Index(final.var_names.astype(str))
    candidate_set = set(map(str, candidate_genes))
    final_set = set(map(str, final_genes))
    final.var["hvgdecision_candidate"] = genes.isin(candidate_set)
    final.var["highly_variable"] = genes.isin(final_set)
    final.var["hvgdecision_final"] = final.var["highly_variable"].astype(bool)
    final.var["hvgdecision_risk_flagged"] = False
    final.var["hvgdecision_marker_protected"] = False
    final.var["hvgdecision_harmful"] = False
    final.var["hvgdecision_removed"] = False
    final.var["hvgdecision_input_rank"] = np.nan
    final.var["hvgdecision_risk_score"] = np.nan
    final.var["hvgdecision_reason"] = "not in recommended base HVG panel"

    if not decision.empty:
        d = decision.copy()
        d["gene"] = d["gene"].astype(str)
        d = d.drop_duplicates("gene", keep="first").set_index("gene")
        shared = genes.intersection(d.index)
        if len(shared):
            final.var.loc[shared, "hvgdecision_input_rank"] = pd.to_numeric(
                d.loc[shared, "input_rank"], errors="coerce"
            ).to_numpy()
            if "risk_score" in d:
                final.var.loc[shared, "hvgdecision_risk_score"] = pd.to_numeric(
                    d.loc[shared, "risk_score"], errors="coerce"
                ).to_numpy()
            if "risk_flagged" in d:
                final.var.loc[shared, "hvgdecision_risk_flagged"] = _coerce_bool_series(
                    d.loc[shared, "risk_flagged"]
                ).to_numpy()
            if "hard_replicated_marker_protection" in d:
                final.var.loc[shared, "hvgdecision_marker_protected"] = _coerce_bool_series(
                    d.loc[shared, "hard_replicated_marker_protection"]
                ).to_numpy()
            if "harmful_gene" in d:
                final.var.loc[shared, "hvgdecision_harmful"] = _coerce_bool_series(
                    d.loc[shared, "harmful_gene"]
                ).to_numpy()
            if "final_action" in d:
                final.var.loc[shared, "hvgdecision_removed"] = d.loc[
                    shared, "final_action"
                ].astype(str).eq("remove").to_numpy()
            if "decision_reason" in d:
                final.var.loc[shared, "hvgdecision_reason"] = d.loc[
                    shared, "decision_reason"
                ].fillna("").astype(str).to_numpy()

    final.uns["hvgdecision"] = dict(audit)
    return final


@dataclass
class CountSourceResult:
    """Result of the raw-count source audit."""

    source: str | None
    location: str
    audit: pd.DataFrame
    valid: bool
    error: str = ""
    matrix: Any = None
    gene_names: Any = None

    def __repr__(self) -> str:
        return (
            "CountSourceResult("
            f"valid={self.valid}, source={self.source!r}, location={self.location!r})"
        )


@dataclass
class HVGRefinementResult:
    """Auditable decision returned by :meth:`HVGStudy.refine_hvg`."""

    decision_table: pd.DataFrame
    composition_audit: pd.DataFrame
    audit: dict[str, Any]
    output_dir: Path
    adata: Any

    @property
    def recommended_n_hvg(self) -> int:
        return int(self.audit["recommended_n_hvg"])

    @property
    def removed_genes(self) -> list[str]:
        part = self.decision_table.loc[self.decision_table["final_action"].eq("remove")]
        return part["gene"].astype(str).tolist()

    @property
    def retained_genes(self) -> list[str]:
        part = self.decision_table.loc[self.decision_table["in_final_panel"].fillna(False)]
        return (
            part.sort_values("input_rank", kind="mergesort")["gene"].astype(str).tolist()
        )

    @property
    def harmful_genes(self) -> list[str]:
        if "harmful_gene" not in self.decision_table:
            return self.removed_genes
        part = self.decision_table.loc[_coerce_bool_series(self.decision_table["harmful_gene"])]
        return part.sort_values("input_rank", kind="mergesort")["gene"].astype(str).tolist()

    @property
    def final_hvg_genes(self) -> list[str]:
        return self.retained_genes

    @property
    def final_n_hvg(self) -> int:
        return len(self.retained_genes)

    def __repr__(self) -> str:
        return (
            "HVGRefinementResult("
            f"recommended_n_hvg={self.recommended_n_hvg}, "
            f"removed={len(self.removed_genes)}, "
            f"mode={self.audit.get('mode')!r}, "
            f"control={self.audit.get('actual_control_mode')!r}, "
            f"output_dir='{self.output_dir}')"
        )


@dataclass
class BudgetSearchResult:
    """Auditable result returned by :meth:`HVGStudy.find_best_hvg`.

    ``best_n_hvg`` is the budget-search recommendation *before* targeted risk
    deletion.  ``final_n_hvg`` is the size of the panel actually returned in
    ``adata`` after the Reference-only risk gate.  Zero deletion is a valid
    outcome, in which case both values are identical.
    """

    recommendation: pd.DataFrame
    results: pd.DataFrame
    gene_table: pd.DataFrame
    composition_audit: pd.DataFrame
    decision_table: pd.DataFrame
    refinement_audit: dict[str, Any]
    output_dir: Path
    adata: Any
    budget_summary: pd.DataFrame | None = None

    @property
    def best_n_hvg(self) -> int:
        return int(self.recommendation.iloc[0]["recommended_minimum_sufficient_hvg"])

    @property
    def observed_peak_n_hvg(self) -> int:
        return int(self.recommendation.iloc[0]["observed_peak_hvg"])

    @property
    def observed_peak_biology_n_hvg(self) -> int:
        value = self.recommendation.iloc[0].get(
            "observed_peak_biology_hvg", self.observed_peak_n_hvg
        )
        return int(value)

    @property
    def selection_status(self) -> str:
        return str(self.recommendation.iloc[0].get("selection_status", "unknown"))

    @property
    def best_genes(self) -> list[str]:
        """Base recommended HVGs before risk deletion."""
        part = self.gene_table.loc[self.gene_table["n_hvg"].eq(self.best_n_hvg)]
        return part.sort_values("gene_rank", kind="mergesort")["gene"].astype(str).tolist()

    @property
    def removed_genes(self) -> list[str]:
        part = self.decision_table.loc[self.decision_table["final_action"].eq("remove")]
        return part.sort_values("input_rank", kind="mergesort")["gene"].astype(str).tolist()

    @property
    def retained_genes(self) -> list[str]:
        part = self.decision_table.loc[self.decision_table["in_final_panel"].fillna(False)]
        return part.sort_values("input_rank", kind="mergesort")["gene"].astype(str).tolist()

    @property
    def final_genes(self) -> list[str]:
        return self.retained_genes

    @property
    def final_n_hvg(self) -> int:
        return len(self.retained_genes)

    @property
    def recommended_n_hvg(self) -> int:
        return self.best_n_hvg

    @property
    def harmful_genes(self) -> list[str]:
        if "harmful_gene" not in self.decision_table:
            return self.removed_genes
        part = self.decision_table.loc[_coerce_bool_series(self.decision_table["harmful_gene"])]
        return part.sort_values("input_rank", kind="mergesort")["gene"].astype(str).tolist()

    @property
    def final_hvg_genes(self) -> list[str]:
        return self.final_genes

    def __repr__(self) -> str:
        n_budgets = int(self.gene_table["n_hvg"].nunique()) if "n_hvg" in self.gene_table else 0
        n_holdouts = (
            int(self.results["held_out_reference_batch"].nunique())
            if "held_out_reference_batch" in self.results
            else 0
        )
        return (
            "BudgetSearchResult("
            f"best_n_hvg={self.best_n_hvg}, "
            f"removed={len(self.removed_genes)}, "
            f"final_n_hvg={self.final_n_hvg}, "
            f"observed_peak_n_hvg={self.observed_peak_n_hvg}, "
            f"budgets={n_budgets}, reference_holdouts={n_holdouts}, "
            f"output_dir='{self.output_dir}')"
        )


class HVGStudy:
    """An in-memory, explicitly split Reference/Query experiment."""

    def __init__(self, adata, config: dict[str, Any]):
        self.adata = adata
        self.config = config

    @property
    def reference_mask(self) -> np.ndarray:
        return role_masks(self.adata, self.config)[0]

    @property
    def query_mask(self) -> np.ndarray:
        return role_masks(self.adata, self.config)[1]

    def inspect(self) -> dict[str, Any]:
        """Audit counts, keys, split sizes, batches and shared cell types."""
        return inspect_dataset(self.config)

    @property
    def mode(self) -> str:
        """Canonical selected mode: ``within_domain`` or ``cross_domain``."""
        return normalize_mode(self.config.get("mode"))

    def run(
        self,
        *,
        hvg_method: str = "seurat_v3",
        n_hvg: int = 2000,
        cross_domain_delete_budget: int = 5,
        return_details: bool = False,
        refinement_options: dict[str, Any] | None = None,
        **within_domain_options: Any,
    ) -> Any:
        """Run the workflow selected by ``study.mode``.

        ``within_domain`` performs Reference-only HVG-budget recommendation and
        final harmful-gene gating. ``cross_domain`` builds a Query HVG panel,
        scores Reference/Query distribution shift without Query labels, and
        removes the requested number of top Cross-domain Rule V3 genes.
        """
        mode = self.mode
        if mode == "within_domain":
            options = dict(within_domain_options)
            options.setdefault("hvg_method", hvg_method)
            options.setdefault("refinement_options", refinement_options)
            options.setdefault("return_details", return_details)
            return self.find_best_hvg(**options)

        if int(n_hvg) < 2:
            raise ValueError("n_hvg must be >= 2 in cross_domain mode")
        config = dict(self.config)
        query = self.query_mask
        query_config = {**config, "hvg_method": hvg_method, "hvg_source": "query"}
        ranking = _hvg_ranking(self.adata, query, query_config, int(n_hvg)).copy()
        ranking["n_hvg"] = int(len(ranking))
        ranking["method"] = str(hvg_method)
        options = {
            "mode": "cross_domain",
            "cross_domain_delete_budget": int(cross_domain_delete_budget),
            "cross_reference_hvg_method": str(hvg_method),
            "allow_internal_independent_refit": True,
            **(refinement_options or {}),
        }
        return self.refine_hvg(
            ranking[["method", "n_hvg", "gene_rank", "gene"]],
            method_name=str(hvg_method),
            initial_n_hvg=int(len(ranking)),
            control_mode="independent",
            run_refinement=True,
            refinement_options=options,
            return_details=return_details,
        )

    def find_best_hvg(
        self,
        budgets: Sequence[int] = DEFAULT_HVG_BUDGETS,
        *,
        hvg_method: str = "seurat_v3",
        selection_mode: str = "auto",
        near_optimal_tolerance: float = 0.02,
        biology_near_optimal_tolerance: float = 0.03,
        donor_instability_guardrail_quantile: float = 0.75,
        budget_max_cells_per_batch_label: int = 200,
        biology_plateau_min_tolerance: float = 0.01,
        biology_plateau_max_tolerance: float = 0.015,
        rare_biology_plateau_min_tolerance: float = 0.015,
        rare_biology_plateau_max_tolerance: float = 0.03,
        stability_mad_multiplier: float = 2.5,
        donor_instability_mad_multiplier: float = 2.5,
        rare_reference_quantile: float = 0.25,
        run_refinement: bool = True,
        refinement_options: dict[str, Any] | None = None,
        return_details: bool = False,
    ) -> Any:
        """Recommend an HVG budget and remove harmful HVGs using Reference only.

        No Harmony, BBKNN, Scanorama, scVI, or scANVI is run inside
        ``find_best_hvg``. Budget recommendation is Reference-only. In v0.9,
        donor-held-out biology preservation is the primary sufficiency gate:
        macro-F1, balanced accuracy, rare-cell macro-F1, and minimum class
        recall are combined into a biology sufficiency score. Within-biology
        donor instability excludes batch-heavy panels, while HVG stability is
        only a robust low-quality guardrail rather than a demand to sit near
        the global stability maximum. The smallest panel on the supported
        biology plateau is selected. Query expression and Query labels never
        participate in budget selection or risk estimation.

        The recommended base panel is then frozen and passed through the
        marker-protected harmful-gene gate.  Zero harmful genes is a valid
        result.  The returned AnnData keeps all genes from the selected count
        source; ``adata.var['highly_variable']`` marks the final HVGs that should
        be supplied to downstream integration software.
        """
        adata = self.adata
        config = dict(self.config)
        if normalize_mode(config.get("mode")) != "within_domain":
            raise RuntimeError(
                "find_best_hvg() is the within_domain budget-search workflow. "
                "For cross_domain mode use study.run(...) or study.refine_hvg(...)."
            )
        reference, query, split = role_masks(adata, config)
        source = "reference"
        source_mask = reference

        requested_selection_mode = str(selection_mode).lower()
        if requested_selection_mode == "nested":
            requested_selection_mode = "truncate"
        if requested_selection_mode not in {"auto", "truncate", "independent"}:
            raise ValueError("selection_mode must be 'auto', 'truncate', or 'independent'")
        if requested_selection_mode == "auto":
            actual_selection_mode = (
                "independent" if _is_seurat_v3_like(hvg_method) else "truncate"
            )
        else:
            actual_selection_mode = requested_selection_mode

        available_n_vars = len(count_var_names(adata, config.get("counts_layer")))
        valid_budgets = sorted(
            {int(value) for value in budgets if 2 <= int(value) <= available_n_vars - 1}
        )
        if len(valid_budgets) < 2:
            raise ValueError(
                "At least two valid HVG budgets are required after filtering to "
                f"2..{available_n_vars - 1}; observed={list(budgets)}"
            )
        if float(near_optimal_tolerance) < 0:
            raise ValueError("near_optimal_tolerance must be >= 0")
        if float(biology_near_optimal_tolerance) < 0:
            raise ValueError("biology_near_optimal_tolerance must be >= 0")
        if not 0.0 <= float(donor_instability_guardrail_quantile) <= 1.0:
            raise ValueError("donor_instability_guardrail_quantile must be in [0, 1]")
        if int(budget_max_cells_per_batch_label) < 1:
            raise ValueError("budget_max_cells_per_batch_label must be >= 1")
        if not 0.0 <= float(rare_reference_quantile) <= 1.0:
            raise ValueError("rare_reference_quantile must be in [0, 1]")
        if not 0.0 <= float(biology_plateau_min_tolerance) <= float(biology_plateau_max_tolerance):
            raise ValueError(
                "biology plateau tolerances must satisfy 0 <= min <= max"
            )
        if not 0.0 <= float(rare_biology_plateau_min_tolerance) <= float(rare_biology_plateau_max_tolerance):
            raise ValueError(
                "rare biology plateau tolerances must satisfy 0 <= min <= max"
            )
        if float(stability_mad_multiplier) < 0:
            raise ValueError("stability_mad_multiplier must be >= 0")
        if float(donor_instability_mad_multiplier) < 0:
            raise ValueError("donor_instability_mad_multiplier must be >= 0")

        reference_batches = sorted(
            adata.obs.loc[reference, config["batch_key"]].astype(str).unique().tolist()
        )
        if len(reference_batches) < 2:
            raise ValueError(
                "find_best_hvg requires at least two Reference batches/donors so "
                "HVG reproducibility can be assessed without using Query."
            )

        output = Path(config["output_dir"]).expanduser().resolve() / "00_budget_search"
        output.mkdir(parents=True, exist_ok=True)
        config.update(
            {
                "hvg_source": "reference",
                "hvg_method": hvg_method,
                "budget_rare_reference_quantile": float(rare_reference_quantile),
            }
        )
        budget_sample_positions = _balanced_reference_sample_indices(
            adata,
            reference,
            config,
            max_cells_per_batch_label=int(budget_max_cells_per_batch_label),
        )
        if budget_sample_positions.size < 3:
            raise RuntimeError("Too few sampled Reference cells for budget evaluation")

        # ------------------------------------------------------------
        # Stage 1A: full-Reference candidate panels.
        # Seurat-v3-like selectors are independently refit for every budget;
        # strictly nested selectors may reuse one maximum ranking by truncation.
        # ------------------------------------------------------------
        ranking_tables: dict[int, pd.DataFrame] = {}
        if actual_selection_mode == "truncate":
            maximum_ranking = _hvg_ranking(adata, source_mask, config, max(valid_budgets))
            for budget in valid_budgets:
                part = maximum_ranking.head(budget).copy()
                part["gene_rank"] = np.arange(1, len(part) + 1)
                ranking_tables[budget] = part
        else:
            for budget in tqdm(valid_budgets, desc="Full-Reference independent HVG fits"):
                ranking_tables[budget] = _hvg_ranking(adata, source_mask, config, budget)

        feature_method = {
            "seurat_v3": "Scanpy_SeuratV3_batch_aware",
            "seurat": "Scanpy_Seurat_batch_aware",
            "cell_ranger": "Scanpy_CellRanger_batch_aware",
        }.get(hvg_method, f"Scanpy_{hvg_method}_batch_aware")
        gene_table = pd.concat(
            [
                ranking_tables[budget][["gene", "gene_rank"]].assign(
                    dataset=config["dataset_name"],
                    feature_method=feature_method,
                    n_hvg=budget,
                )
                for budget in valid_budgets
            ],
            ignore_index=True,
        )
        gene_table = gene_table[["dataset", "feature_method", "n_hvg", "gene_rank", "gene"]]
        write_csv(gene_table, output / "hvg_budget_gene_long_table.csv")

        composition_audit = _budget_composition_audit(ranking_tables)
        composition_audit["requested_selection_mode"] = requested_selection_mode
        composition_audit["actual_selection_mode"] = actual_selection_mode
        write_csv(composition_audit, output / "hvg_budget_composition_audit.csv")

        # ------------------------------------------------------------
        # Stage 1B: Reference-only leave-one-batch-out HVG stability.
        # For each held-out Reference donor, fit HVGs on the remaining
        # Reference donors and ask how reproducible the full-Reference panel is.
        # Query is not touched.
        # ------------------------------------------------------------
        batch_values = adata.obs[config["batch_key"]].astype(str).to_numpy()
        stability_rows: list[dict[str, Any]] = []
        biology_rows: list[dict[str, Any]] = []
        for held_out in tqdm(reference_batches, desc="Reference-only HVG stability + biology"):
            training_mask = reference & (batch_values != held_out)
            if int(training_mask.sum()) < 3:
                raise ValueError(
                    f"Too few Reference cells remain after holding out batch {held_out!r}"
                )

            validation_rankings: dict[int, pd.DataFrame] = {}
            if actual_selection_mode == "truncate":
                maximum = _hvg_ranking(adata, training_mask, config, max(valid_budgets))
                for budget in valid_budgets:
                    part = maximum.head(budget).copy()
                    part["gene_rank"] = np.arange(1, len(part) + 1)
                    validation_rankings[budget] = part
            else:
                for budget in valid_budgets:
                    validation_rankings[budget] = _hvg_ranking(
                        adata, training_mask, config, budget
                    )

            for budget in valid_budgets:
                base_genes = (
                    ranking_tables[budget]
                    .sort_values("gene_rank", kind="mergesort")["gene"]
                    .astype(str)
                    .tolist()
                )
                validation_genes = (
                    validation_rankings[budget]
                    .sort_values("gene_rank", kind="mergesort")["gene"]
                    .astype(str)
                    .tolist()
                )
                overlap_n, overlap_fraction, jaccard, weighted = _weighted_panel_recovery(
                    base_genes, validation_genes
                )
                stability_rows.append(
                    {
                        "dataset": config["dataset_name"],
                        "feature_method": feature_method,
                        "n_hvg": int(budget),
                        "held_out_reference_batch": str(held_out),
                        "training_reference_cells": int(training_mask.sum()),
                        "overlap_n": int(overlap_n),
                        "overlap_fraction": float(overlap_fraction),
                        "jaccard": float(jaccard),
                        "rank_weighted_recovery": float(weighted),
                        "stability_score": float(0.5 * overlap_fraction + 0.5 * weighted),
                        "query_expression_used": False,
                        "query_labels_used": False,
                    }
                )

                # Nested Reference-only biology validation: the HVG panel used
                # to predict the held-out Reference donor was fitted WITHOUT
                # that donor (validation_genes), so biology sufficiency does
                # not benefit from held-out feature-selection leakage.
                biology = _reference_biology_transfer(
                    adata,
                    budget_sample_positions,
                    validation_genes,
                    config,
                    [str(held_out)],
                )
                if not biology.empty:
                    row = biology.iloc[0].to_dict()
                    row.update(
                        {
                            "dataset": config["dataset_name"],
                            "feature_method": feature_method,
                            "n_hvg": int(budget),
                            "feature_panel_fit_excludes_held_out_batch": True,
                            "query_expression_used": False,
                            "query_labels_used": False,
                        }
                    )
                    biology_rows.append(row)

        results = pd.DataFrame(stability_rows)
        biology_by_batch = pd.DataFrame(biology_rows)

        donor_instability_rows: list[dict[str, Any]] = []
        for budget in tqdm(valid_budgets, desc="Reference-only donor guardrail"):
            genes = (
                ranking_tables[budget]
                .sort_values("gene_rank", kind="mergesort")["gene"]
                .astype(str)
                .tolist()
            )
            donor_instability_rows.append(
                {
                    "n_hvg": int(budget),
                    "within_biology_donor_instability": _within_biology_donor_instability(
                        adata,
                        budget_sample_positions,
                        genes,
                        config,
                    ),
                }
            )

        if biology_by_batch.empty:
            raise RuntimeError(
                "No Reference-only donor-held-out biology transfer could be estimated"
            )
        write_csv(
            biology_by_batch,
            output / "reference_biology_transfer_by_batch.csv",
        )
        results = results.merge(
            biology_by_batch[
                [
                    "n_hvg",
                    "held_out_reference_batch",
                    "biology_eval_cells_total",
                    "biology_eval_cells_known",
                    "biology_label_coverage",
                    "biology_macro_f1",
                    "biology_balanced_accuracy",
                    "biology_rare_macro_f1",
                    "biology_rare_eval_cells",
                    "biology_rare_label_count",
                    "biology_min_class_recall",
                    "biology_transfer_score",
                    "biology_sufficiency_score",
                    "feature_panel_fit_excludes_held_out_batch",
                ]
            ],
            on=["n_hvg", "held_out_reference_batch"],
            how="left",
            validate="one_to_one",
        )
        write_csv(results, output / "reference_budget_metrics_by_batch.csv")
        donor_instability = pd.DataFrame(donor_instability_rows)
        write_csv(
            donor_instability,
            output / "reference_within_biology_donor_instability.csv",
        )
        write_csv(results, output / "reference_hvg_stability_by_batch.csv")
        stability_summary = results.groupby(
            ["dataset", "feature_method", "n_hvg"], as_index=False
        ).agg(
            n_reference_holdouts=("held_out_reference_batch", "nunique"),
            overlap_fraction_mean=("overlap_fraction", "mean"),
            overlap_fraction_min=("overlap_fraction", "min"),
            jaccard_mean=("jaccard", "mean"),
            rank_weighted_recovery_mean=("rank_weighted_recovery", "mean"),
            stability_score_mean=("stability_score", "mean"),
            stability_score_min=("stability_score", "min"),
            stability_score_std=("stability_score", "std"),
        )
        write_csv(stability_summary, output / "reference_hvg_stability_summary.csv")

        if stability_summary.empty or stability_summary["stability_score_mean"].isna().all():
            raise RuntimeError("Reference-only HVG stability could not be estimated")

        biology_summary = biology_by_batch.groupby("n_hvg", as_index=False).agg(
            biology_n_reference_holdouts=("held_out_reference_batch", "nunique"),
            biology_label_coverage_mean=("biology_label_coverage", "mean"),
            biology_macro_f1_mean=("biology_macro_f1", "mean"),
            biology_macro_f1_min=("biology_macro_f1", "min"),
            biology_balanced_accuracy_mean=("biology_balanced_accuracy", "mean"),
            biology_balanced_accuracy_min=("biology_balanced_accuracy", "min"),
            biology_rare_n_reference_holdouts=("biology_rare_macro_f1", "count"),
            biology_rare_macro_f1_mean=("biology_rare_macro_f1", "mean"),
            biology_rare_macro_f1_min=("biology_rare_macro_f1", "min"),
            biology_rare_macro_f1_std=("biology_rare_macro_f1", "std"),
            biology_min_class_recall_mean=("biology_min_class_recall", "mean"),
            biology_min_class_recall_min=("biology_min_class_recall", "min"),
            biology_transfer_score_mean=("biology_transfer_score", "mean"),
            biology_transfer_score_min=("biology_transfer_score", "min"),
            biology_transfer_score_std=("biology_transfer_score", "std"),
            biology_sufficiency_score_mean=("biology_sufficiency_score", "mean"),
            biology_sufficiency_score_min=("biology_sufficiency_score", "min"),
            biology_sufficiency_score_std=("biology_sufficiency_score", "std"),
        )

        summary = (
            stability_summary
            .merge(biology_summary, on="n_hvg", how="left", validate="one_to_one")
            .merge(donor_instability, on="n_hvg", how="left", validate="one_to_one")
            .rename(
                columns={
                    "within_biology_donor_instability":
                        "within_biology_donor_instability_mean"
                }
            )
        )

        summary, recommendation = _recommend_joint_reference_budget(
            summary,
            stability_tolerance=float(near_optimal_tolerance),
            biology_tolerance=float(biology_near_optimal_tolerance),
            donor_guardrail_quantile=float(donor_instability_guardrail_quantile),
            biology_plateau_min_tolerance=float(biology_plateau_min_tolerance),
            biology_plateau_max_tolerance=float(biology_plateau_max_tolerance),
            rare_biology_plateau_min_tolerance=float(rare_biology_plateau_min_tolerance),
            rare_biology_plateau_max_tolerance=float(rare_biology_plateau_max_tolerance),
            stability_mad_multiplier=float(stability_mad_multiplier),
            donor_instability_mad_multiplier=float(donor_instability_mad_multiplier),
        )
        best_n_hvg = int(
            recommendation.iloc[0]["recommended_minimum_sufficient_hvg"]
        )
        peak_budget = int(recommendation.iloc[0]["observed_peak_hvg"])
        write_csv(summary, output / "reference_hvg_budget_joint_summary.csv")
        write_csv(recommendation, output / "hvg_budget_recommendation.csv")
        write_json(
            {
                "version": "0.9.0",
                "dataset": config["dataset_name"],
                "split": split,
                "counts_source": config.get(
                    "counts_source_label", count_source_label(config.get("counts_layer"))
                ),
                "hvg_method": hvg_method,
                "requested_selection_mode": requested_selection_mode,
                "actual_selection_mode": actual_selection_mode,
                "budgets_evaluated": valid_budgets,
                "reference_validation_batches": reference_batches,
                "budget_selection_basis": (
                    "reference_biology_plateau_primary_rare_cell_guard_"
                    "donor_instability_outlier_guard_stability_quality_floor"
                ),
                "biology_near_optimal_tolerance_legacy_cap": float(
                    biology_near_optimal_tolerance
                ),
                "biology_plateau_min_tolerance": float(biology_plateau_min_tolerance),
                "biology_plateau_max_tolerance": float(biology_plateau_max_tolerance),
                "rare_biology_plateau_min_tolerance": float(
                    rare_biology_plateau_min_tolerance
                ),
                "rare_biology_plateau_max_tolerance": float(
                    rare_biology_plateau_max_tolerance
                ),
                "stability_mad_multiplier": float(stability_mad_multiplier),
                "donor_instability_mad_multiplier": float(
                    donor_instability_mad_multiplier
                ),
                "rare_reference_quantile": float(rare_reference_quantile),
                "donor_instability_guardrail_quantile": float(
                    donor_instability_guardrail_quantile
                ),
                "budget_max_cells_per_batch_label": int(
                    budget_max_cells_per_batch_label
                ),
                "integration_run_during_selection": False,
                "query_expression_used_for_selection": False,
                "query_labels_used_for_selection": False,
            },
            output / "design_audit.json",
        )

        best_ranking = (
            ranking_tables[best_n_hvg]
            .sort_values("gene_rank", kind="mergesort")[["gene", "gene_rank"]]
            .copy()
        )
        best_genes = best_ranking["gene"].astype(str).tolist()

        # ------------------------------------------------------------
        # Stage 2: Reference-only harmful-gene refinement.
        # ------------------------------------------------------------
        decision_output = Path(config["output_dir"]).expanduser().resolve() / "01_hvg_decision"
        decision_output.mkdir(parents=True, exist_ok=True)
        base_path = decision_output / "_working_base_hvg.csv"
        best_ranking.assign(method=feature_method).to_csv(base_path, index=False)
        refine_config = {
            **config,
            "output_dir": str(decision_output),
            "hvg_source": "reference",
            "hvg_method": hvg_method,
            "initial_hvg": int(best_n_hvg),
            "base_gene_table": str(base_path),
            "base_feature_method": feature_method,
            "base_gene_column": "gene",
            "base_rank_column": "gene_rank",
            "base_method_column": "method",
            "same_n_control_mode": actual_selection_mode,
            "allow_internal_independent_refit": True,
            "write_diagnostics": False,
            **(refinement_options or {}),
        }
        refine_audit = refine_panel(refine_config)
        refine_output = decision_output / "01_refine"
        risk = pd.read_csv(refine_output / "gene_risk_evidence.csv")
        risk["gene"] = risk["gene"].astype(str)
        statistical_risk = (
            _coerce_bool_series(risk["passes_maxT"])
            & _coerce_bool_series(risk["passes_stability"])
            & _coerce_bool_series(risk["passes_effect_floors"])
        )
        risk["risk_flagged_before_protection"] = statistical_risk
        harmful = set(
            risk.loc[_coerce_bool_series(risk["selected_risk_gene"]), "gene"].astype(str)
        )

        if run_refinement:
            final_genes = [gene for gene in best_genes if gene not in harmful]
            adoption_rule = "targeted_harmful_gene_deletion"
        else:
            final_genes = list(best_genes)
            adoption_rule = "base_budget_only_no_harmful_gene_deletion"
        if len(final_genes) < 2:
            raise RuntimeError("Final HVG panel would contain fewer than two genes")
        final_set = set(final_genes)

        decision = best_ranking.rename(columns={"gene_rank": "input_rank"}).copy()
        decision["gene"] = decision["gene"].astype(str)
        decision["available_in_counts"] = True
        decision = decision.merge(risk, on="gene", how="left", validate="one_to_one")
        decision["risk_flagged"] = _coerce_bool_series(
            decision["risk_flagged_before_protection"]
        )
        decision["harmful_gene"] = _coerce_bool_series(decision["selected_risk_gene"])
        decision["in_final_panel"] = decision["gene"].isin(final_set)
        decision["final_action"] = np.where(
            decision["in_final_panel"], "retain", "remove"
        )
        protected = _coerce_bool_series(
            decision.get(
                "hard_replicated_marker_protection",
                pd.Series(False, index=decision.index),
            )
        )
        reasons = []
        for index, row in decision.iterrows():
            if bool(row["in_final_panel"]):
                if bool(row["risk_flagged"]) and bool(protected.loc[index]):
                    reasons.append("risk signal observed but retained by biology/marker protection")
                elif bool(row["harmful_gene"]) and not run_refinement:
                    reasons.append("harmful gene retained because run_refinement=False")
                else:
                    reasons.append("retained in final HVG panel")
            elif bool(row["harmful_gene"]):
                reasons.append(
                    "removed by Cross-domain Rule V3"
                    if mode == "cross_domain"
                    else "removed as Reference-only harmful HVG after biology protection"
                )
            else:
                reasons.append("not adopted in final panel")
        decision["decision_reason"] = reasons
        decision["input_n_hvg"] = int(best_n_hvg)
        decision["final_n_hvg"] = int(len(final_genes))
        decision["hvg_method"] = hvg_method
        decision["requested_selection_mode"] = requested_selection_mode
        decision["actual_selection_mode"] = actual_selection_mode
        decision["run_refinement"] = bool(run_refinement)
        decision = decision.sort_values("input_rank", kind="mergesort")
        write_csv(decision, decision_output / "hvg_decision_table.csv")
        write_csv(
            pd.DataFrame(
                {"gene_rank": np.arange(1, len(final_genes) + 1), "gene": final_genes}
            ),
            decision_output / "final_hvg_genes.csv",
        )
        write_csv(
            decision.loc[decision["harmful_gene"], ["gene", "input_rank", "risk_score", "decision_reason"]],
            decision_output / "harmful_genes.csv",
        )

        public_audit = {
            "version": "0.9.0",
            "mode": "within_domain",
            "dataset": config["dataset_name"],
            "counts_source": config.get(
                "counts_source_label", count_source_label(config.get("counts_layer"))
            ),
            "split": split,
            "route": "python_internal_hvg",
            "hvg_method": hvg_method,
            "requested_selection_mode": requested_selection_mode,
            "actual_selection_mode": actual_selection_mode,
            "budget_selection_basis": (
                "reference_joint_hvg_stability_cross_donor_biology_"
                "within_biology_donor_guardrail"
            ),
            "recommended_budget_metrics": recommendation.iloc[0].to_dict(),
            "integration_run_during_selection": False,
            "recommended_base_n_hvg": int(best_n_hvg),
            "observed_peak_n_hvg": int(peak_budget),
            "risk_flagged_before_protection": int(decision["risk_flagged"].sum()),
            "harmful_genes": int(len(harmful)),
            "removed_genes": int(len(best_genes) - len(final_genes)),
            "final_n_hvg": int(len(final_genes)),
            "run_refinement": bool(run_refinement),
            "adoption_rule": adoption_rule,
            "zero_deletion_allowed": True,
            "query_labels_used_for_selection_or_risk": False,
            "query_expression_used_for_selection_or_risk": False,
            "biology_keys": list(config.get("biology_keys", [config.get("label_key")])),
            "protected_genes": list(config.get("protected_genes", [])),
            "candidate_hvg_genes": list(best_genes),
            "harmful_gene_names": [gene for gene in best_genes if gene in harmful],
            "final_hvg_genes": list(final_genes),
            "final_object_gene_space": "all genes from selected raw-count source",
            "final_hvg_column": "highly_variable",
            "harmful_gene_column": "hvgdecision_harmful",
            "raw_counts_layer": "counts",
            "budget_search_output_dir": str(output),
            "decision_output_dir": str(decision_output),
            "risk_gate": {
                "permutation_design": refine_audit.get("permutation_design"),
                "permutation_alpha": refine_audit.get("permutation_alpha"),
                "bootstrap_min_pass_fraction": refine_audit.get("bootstrap_min_pass_fraction"),
                "leakage_z_floor": refine_audit.get("leakage_z_floor"),
                "risk_z_floor": refine_audit.get("risk_z_floor"),
                "biology_z_ceiling": refine_audit.get("biology_z_ceiling"),
            },
        }
        write_json(public_audit, decision_output / "decision.json")
        base_path.unlink(missing_ok=True)

        final_adata = _make_final_adata(
            adata,
            config,
            reference,
            query,
            decision,
            best_genes,
            final_genes,
            public_audit,
        )
        details = BudgetSearchResult(
            recommendation=recommendation,
            results=results,
            gene_table=gene_table,
            composition_audit=composition_audit,
            decision_table=decision,
            refinement_audit=public_audit,
            output_dir=decision_output,
            adata=final_adata,
            budget_summary=summary,
        )
        self.last_result = details
        return details if return_details else final_adata

    def refine_hvg(
        self,
        hvg: Any,
        *,
        method_name: str = "external",
        selector: str | None = None,
        gene_column: str = "gene",
        rank_column: str = "gene_rank",
        budget_column: str = "n_hvg",
        method_column: str = "method",
        initial_n_hvg: int | None = None,
        control_mode: str = "auto",
        independent_same_n: Any | None = None,
        run_refinement: bool | None = None,
        refinement_options: dict[str, Any] | None = None,
        return_details: bool = False,
    ) -> Any:
        """Refine a frozen HVG list or R-exported HVG table without reselecting it.

        The supplied gene order is preserved.  When multiple ``n_hvg`` budgets
        are present, their composition is audited automatically.  ``auto`` uses
        independent same-N logic for Seurat-v3-like/non-nested selectors and
        simple truncation for strict-prefix selectors.  ``run_refinement`` may
        override whether the final panel removes the flagged risk genes
        themselves (True) or only adopts the recommended panel size by taking
        the original top-N genes (False).
        """
        if selector is not None:
            method_name = str(selector)
        adata = self.adata
        config = dict(self.config)
        mode = normalize_mode(config.get("mode"))
        reference, query, split = role_masks(adata, config)
        source = "query" if mode == "cross_domain" else "reference"

        imported = _read_external_hvg_input(
            hvg,
            gene_column=gene_column,
            rank_column=rank_column,
            budget_column=budget_column,
            method_column=method_column,
            method_name=method_name,
        )
        methods = imported["method"].astype(str).unique().tolist()
        if str(method_name) in methods:
            selected_method = str(method_name)
        elif len(methods) == 1:
            selected_method = methods[0]
        else:
            raise ValueError(
                f"External HVG table contains multiple methods {methods}; "
                "set method_name to one of them"
            )
        imported = imported.loc[imported["method"].eq(selected_method)].copy()
        budgets = sorted(imported["n_hvg"].astype(int).unique().tolist())
        if not budgets:
            raise ValueError("External HVG input contains no genes")
        base_budget = int(initial_n_hvg) if initial_n_hvg is not None else int(max(budgets))
        if base_budget not in budgets:
            if len(budgets) == 1 and len(imported) >= base_budget:
                imported = (
                    imported.sort_values(["gene_rank", "_input_order"], kind="mergesort")
                    .head(base_budget)
                    .copy()
                )
                imported["n_hvg"] = base_budget
                imported["gene_rank"] = np.arange(1, len(imported) + 1)
                budgets = [base_budget]
            else:
                raise ValueError(
                    f"initial_n_hvg={base_budget} is not present in external budgets {budgets}"
                )

        ranking_tables: dict[int, pd.DataFrame] = {}
        for budget in budgets:
            part = imported.loc[imported["n_hvg"].eq(budget)].copy()
            part = part.sort_values(["gene_rank", "_input_order"], kind="mergesort")
            ranking_tables[int(budget)] = part[["gene", "gene_rank"]].reset_index(drop=True)
        composition = _budget_composition_audit(ranking_tables)
        if len(composition) > 1:
            lower = composition.loc[composition["n_hvg"].lt(max(budgets))]
            detected_strict_nested = bool(
                not lower.empty and lower["strict_prefix_truncation"].all()
            )
            detected_non_nested = bool(
                not lower.empty and (~lower["strict_prefix_truncation"]).any()
            )
        else:
            detected_strict_nested = not _is_seurat_v3_like(selected_method)
            detected_non_nested = _is_seurat_v3_like(selected_method)
        composition["method"] = selected_method
        composition["detected_strict_nested"] = detected_strict_nested
        composition["detected_non_nested"] = detected_non_nested

        requested_control_mode = str(control_mode).lower()
        if requested_control_mode == "nested":
            requested_control_mode = "truncate"
        if requested_control_mode not in {"auto", "truncate", "independent"}:
            raise ValueError("control_mode must be 'auto', 'truncate', or 'independent'")
        if requested_control_mode == "auto":
            actual_control_mode = (
                "independent"
                if independent_same_n is not None
                or detected_non_nested
                or _is_seurat_v3_like(selected_method)
                else "truncate"
            )
        else:
            actual_control_mode = requested_control_mode
        if run_refinement is None:
            actual_run_refinement = (
                True if mode == "cross_domain" else actual_control_mode == "independent"
            )
        else:
            actual_run_refinement = bool(run_refinement)

        output = Path(config["output_dir"]).expanduser().resolve() / "01_hvg_decision"
        output.mkdir(parents=True, exist_ok=True)
        composition["requested_control_mode"] = requested_control_mode
        composition["actual_control_mode"] = actual_control_mode
        composition["run_refinement"] = actual_run_refinement
        write_csv(composition, output / "hvg_composition_audit.csv")

        available = set(count_var_names(adata, config.get("counts_layer")))
        base_all = imported.loc[imported["n_hvg"].eq(base_budget)].copy()
        base_all = base_all.sort_values(["gene_rank", "_input_order"], kind="mergesort")
        base_all = base_all.drop_duplicates("gene", keep="first").copy()
        base_all["available_in_counts"] = base_all["gene"].isin(available)
        available_base = base_all.loc[base_all["available_in_counts"]].copy()
        if len(available_base) < 2:
            raise ValueError("Fewer than two supplied HVGs are present in the raw-count matrix")
        available_base["gene_rank"] = np.arange(1, len(available_base) + 1)
        base_path = output / "_working_base_hvg.csv"
        available_base[["gene", "gene_rank"]].assign(method=selected_method).to_csv(
            base_path, index=False
        )

        independent_path: Path | None = None
        if independent_same_n is not None:
            independent_frame = _read_external_hvg_input(
                independent_same_n,
                gene_column=gene_column,
                rank_column=rank_column,
                budget_column=budget_column,
                method_column=method_column,
                method_name=selected_method,
            )
            if selected_method in independent_frame["method"].astype(str).unique():
                independent_frame = independent_frame.loc[
                    independent_frame["method"].astype(str).eq(selected_method)
                ].copy()
            independent_path = output / "_working_independent_same_n.csv"
            independent_frame[["method", "n_hvg", "gene_rank", "gene"]].to_csv(
                independent_path, index=False
            )

        refine_config = {
            **config,
            "output_dir": str(output),
            "hvg_source": source,
            "hvg_method": selected_method,
            "initial_hvg": int(len(available_base)),
            "base_gene_table": str(base_path),
            "base_feature_method": selected_method,
            "base_gene_column": "gene",
            "base_rank_column": "gene_rank",
            "base_method_column": "method",
            "same_n_control_mode": actual_control_mode,
            "selection_is_non_nested": detected_non_nested,
            "write_diagnostics": False,
            **(refinement_options or {}),
        }
        if independent_path is not None:
            refine_config["independent_same_n_table"] = str(independent_path)
            refine_config["independent_gene_column"] = "gene"
            refine_config["independent_rank_column"] = "gene_rank"
            refine_config["independent_budget_column"] = "n_hvg"

        refine_audit = refine_panel(refine_config)
        refine_output = output / "01_refine"
        risk = pd.read_csv(refine_output / "gene_risk_evidence.csv")
        risk["gene"] = risk["gene"].astype(str)
        if "risk_flagged_before_protection" in risk:
            risk["risk_flagged_before_protection"] = _coerce_bool_series(
                risk["risk_flagged_before_protection"]
            )
        else:
            risk["risk_flagged_before_protection"] = (
                _coerce_bool_series(risk["passes_maxT"])
                & _coerce_bool_series(risk["passes_stability"])
                & _coerce_bool_series(risk["passes_effect_floors"])
            )
        harmful = set(
            risk.loc[_coerce_bool_series(risk["selected_risk_gene"]), "gene"].astype(str)
        )
        base_genes = available_base["gene"].astype(str).tolist()
        risk_refined = [gene for gene in base_genes if gene not in harmful]
        recommended_n = len(risk_refined)
        if recommended_n < 2:
            raise RuntimeError("Recommended HVG panel would contain fewer than two genes")
        if actual_run_refinement:
            final_genes = risk_refined
            adoption_rule = "targeted_risk_deletion"
        else:
            final_genes = base_genes[:recommended_n]
            adoption_rule = "count_only_top_n_truncation"
        final_set = set(final_genes)

        decision = base_all.rename(columns={"gene_rank": "input_rank"})[
            ["gene", "input_rank", "available_in_counts"]
        ].copy()
        decision = decision.merge(risk, on="gene", how="left", validate="one_to_one")
        decision["risk_flagged"] = _coerce_bool_series(
            decision["risk_flagged_before_protection"]
        )
        decision["harmful_gene"] = _coerce_bool_series(decision["selected_risk_gene"])
        decision["in_final_panel"] = decision["gene"].isin(final_set)
        decision["final_action"] = np.where(
            ~decision["available_in_counts"],
            "not_found",
            np.where(decision["in_final_panel"], "retain", "remove"),
        )

        protected = decision.get(
            "hard_replicated_marker_protection", pd.Series(False, index=decision.index)
        )
        protected = _coerce_bool_series(protected)
        reasons = []
        for index, row in decision.iterrows():
            if not bool(row["available_in_counts"]):
                reasons.append("gene not found in selected raw-count source")
            elif bool(row["in_final_panel"]):
                if bool(row["risk_flagged"]) and bool(protected.loc[index]):
                    reasons.append("risk signal observed but retained by biology/marker protection")
                elif bool(row["harmful_gene"]) and not actual_run_refinement:
                    reasons.append("harmful gene retained because user selected count-only mode")
                else:
                    reasons.append("retained in final HVG panel")
            elif actual_run_refinement and bool(row["harmful_gene"]):
                reasons.append("removed as marker-protected Reference-only harmful HVG")
            elif not actual_run_refinement:
                reasons.append("removed by recommended top-N truncation")
            else:
                reasons.append("not adopted in final panel")
        decision["decision_reason"] = reasons
        decision["method"] = selected_method
        decision["input_n_hvg"] = int(base_budget)
        decision["recommended_n_hvg"] = int(recommended_n)
        decision["requested_control_mode"] = requested_control_mode
        decision["actual_control_mode"] = refine_audit.get(
            "actual_control_mode", actual_control_mode
        )
        decision["run_refinement"] = actual_run_refinement
        decision = decision.sort_values("input_rank", kind="mergesort")
        write_csv(decision, output / "hvg_decision_table.csv")
        write_csv(
            pd.DataFrame({"gene_rank": np.arange(1, recommended_n + 1), "gene": final_genes}),
            output / "final_hvg_genes.csv",
        )

        public_audit = {
            "version": "0.9.0",
            "mode": mode,
            "dataset": config["dataset_name"],
            "counts_source": config.get(
                "counts_source_label", count_source_label(config.get("counts_layer"))
            ),
            "split": split,
            "method": selected_method,
            "input_n_hvg": int(base_budget),
            "input_available_hvg": int(len(base_genes)),
            "not_found_genes": int((~base_all["available_in_counts"]).sum()),
            "risk_flagged_before_protection": int(decision["risk_flagged"].sum()),
            "harmful_genes": int(len(harmful)),
            "recommended_n_hvg": int(recommended_n),
            "requested_control_mode": requested_control_mode,
            "actual_control_mode": refine_audit.get("actual_control_mode", actual_control_mode),
            "direct_control_available": bool(
                refine_audit.get("direct_control_available", False)
            ),
            "direct_control_definition": refine_audit.get("direct_control_definition", ""),
            "detected_strict_nested": detected_strict_nested,
            "detected_non_nested": detected_non_nested,
            "run_refinement": actual_run_refinement,
            "adoption_rule": adoption_rule,
            "zero_deletion_allowed": True,
            "biology_keys": list(config.get("biology_keys", [config.get("label_key")])),
            "protected_genes": list(config.get("protected_genes", [])),
            "query_labels_used_for_selection_or_risk": False,
            "query_expression_used_for_selection_or_risk": bool(
                refine_audit.get("query_expression_used_for_selection_or_risk", False)
            ),
            "candidate_hvg_genes": list(base_genes),
            "harmful_gene_names": [gene for gene in base_genes if gene in harmful],
            "final_hvg_genes": list(final_genes),
            "final_object_gene_space": "all genes from selected raw-count source",
            "final_hvg_column": "highly_variable",
            "harmful_gene_column": "hvgdecision_harmful",
            "raw_counts_layer": "counts",
            "output_dir": str(output),
        }
        write_json(public_audit, output / "decision.json")
        write_csv(
            decision.loc[decision["harmful_gene"], ["gene", "input_rank", "risk_score", "decision_reason"]],
            output / "harmful_genes.csv",
        )
        for working_path in (base_path, independent_path):
            if working_path is not None:
                working_path.unlink(missing_ok=True)
        final_adata = _make_final_adata(
            adata,
            config,
            reference,
            query,
            decision,
            base_genes,
            final_genes,
            public_audit,
        )
        details = HVGRefinementResult(
            decision_table=decision,
            composition_audit=composition,
            audit=public_audit,
            output_dir=output,
            adata=final_adata,
        )
        self.last_result = details
        return details if return_details else final_adata


def find_raw_counts(
    adata, counts_layer: Any = "auto", *, source: Any = _UNSET
) -> CountSourceResult:
    """Find and audit the raw integer-count matrix in an AnnData object.

    ``counts_layer`` may be ``"auto"``, any layer name, ``"raw"``/``"raw.X"``,
    or ``None``/``"X"``. Auto mode audits every available location and selects
    the first valid candidate by priority. An explicitly requested source is
    audited on its own and is never silently replaced by another source.
    """
    if source is not _UNSET:
        if counts_layer != "auto":
            raise ValueError("Pass either counts_layer or source, not both")
        counts_layer = source
    if isinstance(counts_layer, CountSourceResult):
        return counts_layer
    if not isinstance(counts_layer, (str, type(None))) or looks_like_count_table_path(
        counts_layer
    ):
        try:
            matrix, gene_names, location, alignment = external_count_matrix(
                adata, counts_layer
            )
            row = audit_counts(matrix, location)
            row["alignment"] = alignment
            row["selected"] = bool(row["valid"])
            error = "" if row["valid"] else str(row["error"])
            return CountSourceResult(
                source="external",
                location=location,
                audit=pd.DataFrame([row]),
                valid=bool(row["valid"]),
                error=error,
                matrix=matrix,
                gene_names=gene_names,
            )
        except Exception as exception:
            location = (
                str(counts_layer)
                if isinstance(counts_layer, (str, Path))
                else f"user object ({type(counts_layer).__name__})"
            )
            row = {
                "source": location,
                "valid": False,
                "selected": False,
                "error": repr(exception),
            }
            return CountSourceResult(
                source="external",
                location=location,
                audit=pd.DataFrame([row]),
                valid=False,
                error=repr(exception),
            )
    requested = counts_layer
    source, audit_rows = resolve_counts_source(adata, requested)
    audit = pd.DataFrame(audit_rows)
    valid = bool(
        not audit.empty
        and "selected" in audit
        and audit["selected"].fillna(False).astype(bool).any()
    )
    if valid:
        selected_row = audit.loc[audit["selected"].fillna(False).astype(bool)].iloc[0]
        location = str(selected_row["source"])
        error = ""
    else:
        location = "NOT_FOUND"
        error = " | ".join(
            f"{row.get('source', 'unknown')}: {row.get('error', 'invalid')}"
            for row in audit_rows
        )
    return CountSourceResult(
        source=source,
        location=location,
        audit=audit,
        valid=valid,
        error=error,
    )


def setup_reference_query(
    adata,
    *,
    mode: str,
    batch_key: str,
    label_key: str,
    reference: str | Sequence[str],
    query: str | Sequence[str],
    split_key: str | None = None,
    biology_keys: Sequence[str] | None = None,
    protected_genes: Sequence[str] | None = None,
    counts_layer: Any = "auto",
    counts: Any = _UNSET,
    dataset_name: str = "dataset",
    output_dir: str | Path = "HVGDecision_results",
) -> HVGStudy:
    """Define a Reference/Query experiment directly from a Scanpy AnnData.

    ``split_key`` defaults to ``batch_key``. ``counts_layer='auto'`` audits all
    candidate locations. When the user explicitly supplies a source, that
    exact source is validated and is never silently replaced.
    """
    mode = normalize_mode(mode)
    biology_keys = list(dict.fromkeys(str(key) for key in (biology_keys or [label_key])))
    if label_key not in biology_keys:
        biology_keys.insert(0, label_key)
    protected_genes = list(dict.fromkeys(str(gene) for gene in (protected_genes or [])))
    for key in (batch_key, label_key, split_key or batch_key, *biology_keys):
        if key not in adata.obs:
            raise KeyError(f"adata.obs is missing required column {key!r}")
    if adata.n_obs < 3:
        raise ValueError("AnnData has fewer than three cells")

    split_key = split_key or batch_key
    reference_values = _values(reference, "reference")
    query_values = _values(query, "query")
    overlap = sorted(set(reference_values) & set(query_values))
    if overlap:
        raise ValueError(f"Reference and Query values overlap: {overlap}")
    observed = set(adata.obs[split_key].astype(str))
    missing_reference = sorted(set(reference_values) - observed)
    missing_query = sorted(set(query_values) - observed)
    if missing_reference or missing_query:
        raise ValueError(
            f"Split values are absent from adata.obs[{split_key!r}]: "
            f"missing_reference={missing_reference}, missing_query={missing_query}"
        )
    selected_roles = adata.obs[split_key].astype(str).isin(reference_values + query_values)
    if adata.obs.loc[selected_roles, batch_key].isna().any():
        raise ValueError("batch_key cannot contain missing values in Reference/Query cells")
    reference_rows = adata.obs[split_key].astype(str).isin(reference_values)
    if adata.obs.loc[reference_rows, label_key].isna().any():
        raise ValueError("label_key cannot contain missing values in Reference cells")
    missing_biology = [
        key for key in biology_keys if adata.obs.loc[reference_rows, key].isna().any()
    ]
    if missing_biology:
        raise ValueError(
            f"biology_keys cannot contain missing values in Reference cells: {missing_biology}"
        )

    if counts is not _UNSET:
        if counts_layer != "auto":
            raise ValueError("Pass either counts or counts_layer, not both")
        counts_layer = counts
    count_result = find_raw_counts(adata, counts_layer)
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    write_csv(count_result.audit, root / "raw_count_source_audit.csv")
    if not count_result.valid:
        raise ValueError(
            "The requested/available expression source is not valid raw counts. "
            f"Inspect {root / 'raw_count_source_audit.csv'}. Details: "
            f"{count_result.error}"
        )
    if count_result.matrix is not None:
        import anndata as ad

        external_genes = pd.Index(count_result.gene_names, name="gene")
        original_var = adata.var.copy()
        original_var.index = original_var.index.astype(str)
        if external_genes.isin(original_var.index).all():
            external_var = original_var.loc[external_genes].copy()
        else:
            external_var = pd.DataFrame(index=external_genes)
        adata = ad.AnnData(
            X=count_result.matrix,
            obs=adata.obs.copy(),
            var=external_var,
        )
        resolved_layer = None
    else:
        resolved_layer = count_result.source
    selected_var_names = count_var_names(adata, resolved_layer)
    if not selected_var_names.is_unique:
        raise ValueError(f"Gene identifiers in {count_result.location} must be unique")
    if len(selected_var_names) < 3:
        raise ValueError("The selected raw-count source has fewer than three genes")

    config: dict[str, Any] = {
        "_adata": adata,
        "mode": mode,
        "dataset_name": str(dataset_name),
        "output_dir": str(Path(output_dir).expanduser().resolve()),
        "counts_layer": resolved_layer,
        "counts_source_label": count_result.location,
        "batch_key": batch_key,
        "label_key": label_key,
        "biology_keys": biology_keys,
        "protected_genes": protected_genes,
        "split_key": split_key,
        "reference_values": reference_values,
        "query_values": query_values,
        "hvg_source": "reference",
    }
    reference_mask, query_mask, split = role_masks(adata, config)
    role = np.full(adata.n_obs, "excluded", dtype=object)
    role[reference_mask] = "reference"
    role[query_mask] = "query"
    adata.obs["_hvgdecision_role"] = pd.Categorical(
        role, categories=["reference", "query", "excluded"]
    )

    root = Path(config["output_dir"])
    write_json(
        {
            "dataset": config["dataset_name"],
            "mode": mode,
            "n_cells": int(adata.n_obs),
            "n_genes": int(len(selected_var_names)),
            "reference_cells": int(reference_mask.sum()),
            "query_cells": int(query_mask.sum()),
            "excluded_cells": int((~(reference_mask | query_mask)).sum()),
            "split": split,
            "batch_key": batch_key,
            "label_key": label_key,
            "biology_keys": biology_keys,
            "protected_genes": protected_genes,
            "counts_source": count_result.location,
        },
        root / "reference_query_setup.json",
    )
    return HVGStudy(adata, config)
