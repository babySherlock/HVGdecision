"""Frozen statistical engine; see docs/ENGINE_PROVENANCE.json."""
import numpy as np
import pandas as pd
from scipy import sparse
import warnings

RISK_ENGINE_VERSION = "scenario_specific_v2.0"
DONOR_REPLICATE_RISK_CONFIG = {
    "min_group_cells": 15,
    "technical_quantile": 0.975,
    "risk_quantile": 0.985,
    "minimum_absolute_risk": 0.30,
    "maximum_biological_support_quantile": 0.45,
    "minimum_direction_agreement": 0.75,
    "minimum_single_donor_dominance": 0.55,
}

def _safe_quantile(values, q, default=0.0):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.quantile(values, q)) if values.size else float(default)

def _percentile_rank(values):
    values = np.asarray(values, dtype=float)
    ranked = pd.Series(values).rank(method="average", pct=True).to_numpy(dtype=float, copy=True)
    ranked[~np.isfinite(ranked)] = 0.0
    return ranked

def _group_means(matrix, codes, n_groups):
    codes = np.asarray(codes, dtype=int)
    n_cells = matrix.shape[0]
    if len(codes) != n_cells:
        raise ValueError("Grouping vector length does not match the count matrix")
    valid = codes >= 0
    indicator = sparse.csr_matrix(
        (
            np.ones(int(valid.sum()), dtype=np.float32),
            (codes[valid], np.flatnonzero(valid)),
        ),
        shape=(n_groups, n_cells),
    )
    group_n = np.bincount(codes[valid], minlength=n_groups).astype(float)
    sums = indicator @ matrix
    if sparse.issparse(sums):
        sums = sums.toarray()
    else:
        sums = np.asarray(sums)
    means = sums / np.maximum(group_n[:, None], 1.0)
    return means.astype(np.float64, copy=False), group_n

def _log_normalized_hvg_matrix(counts, all_genes, hvg_genes):
    all_genes = pd.Index(all_genes.astype(str))
    hvg_genes = list(map(str, hvg_genes))
    if len(hvg_genes) != len(set(hvg_genes)):
        raise ValueError("HVG list contains duplicate genes")
    positions = all_genes.get_indexer(hvg_genes)
    missing = [hvg_genes[i] for i, position in enumerate(positions) if position < 0]
    if missing:
        raise KeyError(f"{len(missing)} HVGs are absent from the count matrix: {missing[:20]}")

    library_size = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    library_size = np.maximum(library_size, 1.0)
    scale = (1.0e4 / library_size).astype(np.float32)
    selected = counts[:, positions]

    if sparse.issparse(selected):
        selected = selected.tocsr().astype(np.float32, copy=True)
        selected = selected.multiply(scale[:, None]).tocsr()
        selected.data = np.log1p(selected.data)
        selected.eliminate_zeros()
        return selected

    selected = np.array(selected, dtype=np.float32, copy=True)
    selected *= scale[:, None]
    np.log1p(selected, out=selected)
    return selected

def _mean_and_variance(matrix):
    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=0)).ravel().astype(float)
        mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel().astype(float)
    else:
        matrix = np.asarray(matrix, dtype=np.float64)
        mean = matrix.mean(axis=0)
        mean_sq = np.square(matrix).mean(axis=0)
    variance = np.maximum(mean_sq - np.square(mean), 1.0e-12)
    return mean, variance

def _compute_replicate_components(
    adata,
    counts,
    genes,
    hvg_genes,
    batch_key,
    label_key,
    *,
    min_group_cells=15,
):
    for key in [batch_key, label_key]:
        if key not in adata.obs.columns:
            raise KeyError(f"Missing required Reference metadata column: {key!r}")

    obs = adata.obs[[batch_key, label_key]].copy()
    obs[batch_key] = obs[batch_key].astype(str)
    obs[label_key] = obs[label_key].astype(str)
    if obs[batch_key].isna().any() or obs[label_key].isna().any():
        raise ValueError("Reference batch/label metadata contains missing values")

    batch_levels = sorted(obs[batch_key].unique().tolist())
    label_levels = sorted(obs[label_key].unique().tolist())
    if len(batch_levels) < 2:
        raise ValueError(
            f"The scenario-specific risk engine needs at least 2 Reference groups in {batch_key!r}"
        )
    if len(label_levels) < 2:
        raise ValueError(
            f"The scenario-specific risk engine needs at least 2 cell types in {label_key!r}"
        )

    batch_code = pd.Categorical(obs[batch_key], categories=batch_levels).codes
    label_code = pd.Categorical(obs[label_key], categories=label_levels).codes
    n_batch = len(batch_levels)
    n_label = len(label_levels)
    combo_code = label_code * n_batch + batch_code
    n_combo = n_label * n_batch

    matrix = _log_normalized_hvg_matrix(counts, pd.Index(genes), hvg_genes)
    overall_mean, total_variance = _mean_and_variance(matrix)
    label_mean, label_n = _group_means(matrix, label_code, n_label)
    combo_mean, combo_n = _group_means(matrix, combo_code, n_combo)

    between_label = (
        label_n[:, None] * np.square(label_mean - overall_mean[None, :])
    ).sum(axis=0) / max(float(label_n.sum()), 1.0)
    celltype_eta2 = np.clip(between_label / total_variance, 0.0, 1.0)

    combo_label = np.repeat(np.arange(n_label), n_batch)
    combo_deviation = combo_mean - label_mean[combo_label]
    combo_valid = combo_n >= int(min_group_cells)
    contribution = combo_n[:, None] * np.square(combo_deviation)
    contribution[~combo_valid, :] = 0.0
    effective_n = max(float(combo_n[combo_valid].sum()), 1.0)
    nuisance_eta2 = np.clip(
        contribution.sum(axis=0) / effective_n / total_variance,
        0.0,
        1.0,
    )
    contribution_total = contribution.sum(axis=0)
    max_group_contribution = contribution.max(axis=0) / np.maximum(
        contribution_total, 1.0e-12
    )

    top_label_index = np.argmax(label_mean, axis=0)
    top_label = np.asarray(label_levels, dtype=object)[top_label_index]
    top_prevalence = label_n[top_label_index] / max(float(label_n.sum()), 1.0)

    n_gene = len(hvg_genes)
    gene_index = np.arange(n_gene)
    effects = np.full((n_batch, n_gene), np.nan, dtype=float)
    for batch_index in range(n_batch):
        rows = np.arange(n_label) * n_batch + batch_index
        valid_label = combo_n[rows] >= int(min_group_cells)
        if int(valid_label.sum()) < 2:
            continue
        valid_rows = rows[valid_label]
        total_n = float(combo_n[valid_rows].sum())
        total_sum = (
            combo_n[valid_rows, None] * combo_mean[valid_rows]
        ).sum(axis=0)
        top_rows = top_label_index * n_batch + batch_index
        top_n = combo_n[top_rows]
        top_mean = combo_mean[top_rows, gene_index]
        other_n = total_n - top_n
        usable = (top_n >= int(min_group_cells)) & (other_n >= int(min_group_cells))
        other_mean = (total_sum - top_n * top_mean) / np.maximum(other_n, 1.0)
        effects[batch_index, usable] = top_mean[usable] - other_mean[usable]

    finite_effect = np.isfinite(effects)
    replicate_count = finite_effect.sum(axis=0).astype(int)
    direction_agreement = np.divide(
        np.nansum(effects > 0.0, axis=0),
        np.maximum(replicate_count, 1),
    )
    absolute_effect = np.abs(effects)
    effect_sum = np.nansum(absolute_effect, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        effect_max = np.nanmax(absolute_effect, axis=0)
        effect_sd = np.nanstd(effects, axis=0)
        effect_mean_abs = np.nanmean(absolute_effect, axis=0)
    effect_max[~np.isfinite(effect_max)] = 0.0
    effect_heterogeneity = np.clip(
        effect_sd / np.maximum(effect_mean_abs, 1.0e-8), 0.0, 5.0
    )
    donor_or_source_dominance = effect_max / np.maximum(effect_sum, 1.0e-8)

    lodo_effect = np.full_like(effects, np.nan)
    for excluded in range(n_batch):
        retained = effects.copy()
        retained[excluded, :] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            lodo_effect[excluded, :] = np.nanmean(retained, axis=0)
    lodo_valid = np.isfinite(lodo_effect)
    lodo_count = lodo_valid.sum(axis=0).astype(int)
    lodo_positive_fraction = np.divide(
        np.nansum(lodo_effect > 0.0, axis=0),
        np.maximum(lodo_count, 1),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        lodo_sd = np.nanstd(lodo_effect, axis=0)
        lodo_mean_abs = np.nanmean(np.abs(lodo_effect), axis=0)
    lodo_instability = np.clip(
        lodo_sd / np.maximum(lodo_mean_abs, 1.0e-8), 0.0, 5.0
    )

    minimum_replicates = max(2, int(np.ceil(n_batch / 2)))
    replicated_marker_protection = (
        (celltype_eta2 >= _safe_quantile(celltype_eta2, 0.90))
        & (direction_agreement >= 0.80)
        & (replicate_count >= minimum_replicates)
        & (lodo_positive_fraction >= 0.80)
    )
    rare_marker_protection = (
        (top_prevalence <= 0.05)
        & (celltype_eta2 >= _safe_quantile(celltype_eta2, 0.75))
        & (direction_agreement >= 0.75)
        & (replicate_count >= minimum_replicates)
    )

    return pd.DataFrame(
        {
            "gene": list(map(str, hvg_genes)),
            "hvg_rank": np.arange(1, len(hvg_genes) + 1, dtype=int),
            "top_supported_celltype": top_label,
            "top_celltype_prevalence": top_prevalence,
            "celltype_eta2": celltype_eta2,
            "nuisance_eta2_within_celltype": nuisance_eta2,
            "max_group_contribution_fraction": max_group_contribution,
            "replicate_count": replicate_count,
            "celltype_effect_direction_agreement": direction_agreement,
            "single_group_effect_dominance": donor_or_source_dominance,
            "celltype_effect_heterogeneity": effect_heterogeneity,
            "lodo_positive_fraction": lodo_positive_fraction,
            "lodo_instability": lodo_instability,
            "replicated_marker_protection": replicated_marker_protection,
            "rare_marker_protection": rare_marker_protection,
        }
    )

def _finish_decision_table(table, protected_genes, mode):
    protected_set = set(map(str, protected_genes or []))
    table = table.copy()
    table["user_protected"] = table["gene"].isin(protected_set)
    table["hard_protected"] = (
        table["replicated_marker_protection"].astype(bool)
        | table["rare_marker_protection"].astype(bool)
        | table["user_protected"].astype(bool)
    )
    table["risk_flagged"] = (
        table["risk_gate_passed"].astype(bool)
        & ~table["hard_protected"].astype(bool)
    )
    table["final_action"] = np.where(table["risk_flagged"], "remove", "keep")

    reasons = np.full(len(table), "kept_below_mode_specific_risk_gate", dtype=object)
    reasons[table["replicated_marker_protection"].to_numpy(dtype=bool)] = (
        "kept_replicated_celltype_marker"
    )
    reasons[table["rare_marker_protection"].to_numpy(dtype=bool)] = (
        "kept_replicated_rare_cell_marker"
    )
    reasons[table["user_protected"].to_numpy(dtype=bool)] = "kept_user_protected"
    reasons[table["risk_flagged"].to_numpy(dtype=bool)] = (
        "removed_mode_specific_nonreproducible_hvg"
    )
    table["decision_reason"] = reasons
    table["algorithm_mode"] = mode
    table["risk_engine_version"] = RISK_ENGINE_VERSION
    table = table.sort_values(
        ["risk_flagged", "risk_score", "hvg_rank"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    harmful = table.loc[table["risk_flagged"], "gene"].astype(str).tolist()
    return table, harmful

def find_harmful_within_protocol_donor_replicate(
    adata,
    counts,
    genes,
    hvg_genes,
    batch_key,
    label_key,
    protected_genes=None,
    seed=20260829,
    risk_config=None,
    **kwargs,
):
    """Reference-only donor-repeatability risk for a tightly matched protocol/run."""
    del seed, kwargs
    config = {**DONOR_REPLICATE_RISK_CONFIG, **(risk_config or {})}
    table = _compute_replicate_components(
        adata,
        counts,
        genes,
        hvg_genes,
        batch_key,
        label_key,
        min_group_cells=config["min_group_cells"],
    )

    nuisance_rank = _percentile_rank(table["nuisance_eta2_within_celltype"])
    dominance_rank = _percentile_rank(table["single_group_effect_dominance"])
    heterogeneity_rank = _percentile_rank(table["celltype_effect_heterogeneity"])
    lodo_rank = _percentile_rank(table["lodo_instability"])
    celltype_rank = _percentile_rank(table["celltype_eta2"])
    direction = table["celltype_effect_direction_agreement"].to_numpy(dtype=float)

    table["donor_nonreproducibility_score"] = (
        0.35 * nuisance_rank
        + 0.25 * dominance_rank
        + 0.20 * heterogeneity_rank
        + 0.20 * lodo_rank
    )
    table["biological_support_score"] = 0.55 * celltype_rank + 0.45 * direction
    table["risk_score"] = (
        table["donor_nonreproducibility_score"]
        - table["biological_support_score"]
    )

    technical_cut = _safe_quantile(
        table["donor_nonreproducibility_score"], config["technical_quantile"]
    )
    risk_cut = max(
        config["minimum_absolute_risk"],
        _safe_quantile(table["risk_score"], config["risk_quantile"]),
    )
    biological_cut = _safe_quantile(
        table["biological_support_score"],
        config["maximum_biological_support_quantile"],
    )

    table["technical_threshold"] = technical_cut
    table["risk_threshold"] = risk_cut
    table["biological_support_ceiling"] = biological_cut
    table["risk_gate_passed"] = (
        (table["donor_nonreproducibility_score"] >= technical_cut)
        & (table["risk_score"] >= risk_cut)
        & (table["biological_support_score"] <= biological_cut)
        & (table["replicate_count"] >= 3)
        & (
            (table["single_group_effect_dominance"] >= config["minimum_single_donor_dominance"])
            | (table["lodo_positive_fraction"] < 0.80)
        )
        & (
            table["celltype_effect_direction_agreement"]
            < config["minimum_direction_agreement"]
        )
    )

    table, harmful = _finish_decision_table(
        table,
        protected_genes,
        mode="within_protocol_donor_replicate",
    )
    lodo_audit = table[
        [
            "gene",
            "replicate_count",
            "single_group_effect_dominance",
            "celltype_effect_direction_agreement",
            "lodo_positive_fraction",
            "lodo_instability",
            "risk_flagged",
            "decision_reason",
        ]
    ].copy()
    return {
        "harmful_genes": harmful,
        "table": table,
        "lodo_audit": lodo_audit,
    }
