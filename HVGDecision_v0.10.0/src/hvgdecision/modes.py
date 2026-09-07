"""Mode-specific HVGDecision refinement rules.

The package exposes two explicit modes:

``within_domain``
    Multi-donor / multi-batch refinement inside a shared experimental domain.
    Query expression and Query labels are not used for the risk decision.

``cross_domain``
    Cross-dataset / cross-technology refinement. Query expression is used only
    for label-free distribution-shift evidence; Query true labels are never
    used by feature selection.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm.auto import tqdm

from .config import require, resolve_output
from .io import (
    count_var_names,
    gene_hash,
    load_adata,
    resolve_config_count_source,
    role_masks,
    subset_count_matrix,
    write_csv,
    write_json,
)
from .statistics import benjamini_hochberg, biology_components, compose_risk, donor_risk_components
from .refine import (
    _external_independent_panel,
    _hvg_ranking,
    _marker_protection,
    _precomputed_ranking,
    _sample_reference,
    _is_seurat_v3_like,
)

VERSION = "0.9.1"

MODE_ALIASES = {
    "within": "within_domain",
    "within_domain": "within_domain",
    "within-domain": "within_domain",
    "multi_donor": "within_domain",
    "multi-donor": "within_domain",
    "cross": "cross_domain",
    "cross_domain": "cross_domain",
    "cross-domain": "cross_domain",
    "cross_dataset": "cross_domain",
    "cross-dataset": "cross_domain",
    "cross_technology": "cross_domain",
    "cross-technology": "cross_domain",
}


def normalize_mode(value: str | None) -> str:
    if value is None or not str(value).strip():
        raise ValueError(
            "HVGDecision mode must be selected explicitly: "
            "mode='within_domain' or mode='cross_domain'."
        )
    key = str(value).strip().lower().replace(" ", "_")
    if key not in MODE_ALIASES:
        raise ValueError(
            f"Unknown mode {value!r}. Choose 'within_domain' or 'cross_domain'."
        )
    return MODE_ALIASES[key]


def _pct_rank(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    out = pd.Series(x).rank(method="average", pct=True).to_numpy(dtype=float)
    out[~np.isfinite(out)] = 0.0
    return out


def _base_ranking(adata, source_mask: np.ndarray, config: dict[str, Any]) -> tuple[pd.DataFrame, str]:
    if config.get("base_gene_table"):
        ranking = _precomputed_ranking(
            config,
            set(count_var_names(adata, config.get("counts_layer", "counts"))),
        )
        return ranking, "precomputed_table"
    return (
        _hvg_ranking(adata, source_mask, config, int(config["initial_hvg"])),
        "scanpy_fit",
    )


def _reference_risk_evidence(
    adata,
    reference: np.ndarray,
    genes: list[str],
    config: dict[str, Any],
    *,
    output: Path,
    prefix: str = "",
    write_diagnostics: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Compute the final Within-domain V1 evidence table on Reference only."""
    matrix, donors, labels, sampled_obs, sampling = _sample_reference(
        adata, reference, genes, config
    )
    if len(np.unique(donors)) < 2:
        raise ValueError("Reference risk estimation requires at least two batches/donors")

    minimum = int(config.get("min_group_cells", 10))
    leakage, instability = donor_risk_components(matrix, donors, labels, minimum)
    biology, consistency = biology_components(matrix, donors, labels, minimum)
    raw_risk, z_values = compose_risk(leakage, instability, biology, consistency)

    # ------------------------------------------------------------
    # Conditional donor-label permutation within biology strata.
    # Per-gene leakage p values, followed by BH-FDR.
    # ------------------------------------------------------------
    n_permutations = int(config.get("n_permutations", 100))
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")
    rng = np.random.default_rng(int(config.get("selection_seed", 20260829)) + 1)
    exceed = np.zeros(len(genes), dtype=np.int64)
    permutation_rows: list[dict[str, Any]] = []
    for permutation in tqdm(
        range(n_permutations),
        desc=f"{prefix or 'Within-domain'} conditional permutation",
    ):
        permuted = donors.copy()
        for label in sorted(np.unique(labels)):
            positions = np.flatnonzero(labels == label)
            if positions.size > 1:
                permuted[positions] = rng.permutation(permuted[positions])
        null_leakage, _ = donor_risk_components(matrix, permuted, labels, minimum)
        exceed += (null_leakage >= leakage).astype(np.int64)
        if write_diagnostics:
            permutation_rows.extend(
                {
                    "permutation": permutation + 1,
                    "gene": gene,
                    "null_donor_leakage": float(value),
                }
                for gene, value in zip(genes, null_leakage)
            )
    permutation_p = (1.0 + exceed) / (n_permutations + 1.0)
    permutation_fdr = benjamini_hochberg(permutation_p)

    # ------------------------------------------------------------
    # Bootstrap recurrence of the effect gate.
    # ------------------------------------------------------------
    n_bootstraps = int(config.get("n_bootstraps", 20))
    bootstrap_fraction = float(config.get("bootstrap_fraction", 0.80))
    if n_bootstraps < 1:
        raise ValueError("n_bootstraps must be >= 1")
    if not 0 < bootstrap_fraction <= 1:
        raise ValueError("bootstrap_fraction must be in (0, 1]")
    leakage_floor = float(config.get("leakage_z_floor", 1.0))
    risk_floor = float(config.get("risk_z_floor", 1.0))
    biology_ceiling = float(config.get("biology_z_ceiling", 0.0))

    strata: dict[tuple[str, str], np.ndarray] = {}
    for donor in sorted(np.unique(donors)):
        for label in sorted(np.unique(labels)):
            positions = np.flatnonzero((donors == donor) & (labels == label))
            if positions.size:
                strata[(str(donor), str(label))] = positions

    bootstrap_passes = np.zeros(len(genes), dtype=np.int64)
    bootstrap_rows: list[dict[str, Any]] = []
    bootstrap_rng = np.random.default_rng(int(config.get("selection_seed", 20260829)) + 2)
    for bootstrap in tqdm(range(n_bootstraps), desc=f"{prefix or 'Within-domain'} bootstrap"):
        sampled: list[int] = []
        for positions in strata.values():
            take = max(1, int(np.ceil(bootstrap_fraction * len(positions))))
            take = min(len(positions), take)
            sampled.extend(bootstrap_rng.choice(positions, take, replace=False).tolist())
        sampled_idx = np.asarray(sorted(sampled), dtype=int)
        boot_leakage, boot_instability = donor_risk_components(
            matrix[sampled_idx], donors[sampled_idx], labels[sampled_idx], minimum
        )
        boot_biology, boot_consistency = biology_components(
            matrix[sampled_idx], donors[sampled_idx], labels[sampled_idx], minimum
        )
        boot_raw, boot_z = compose_risk(
            boot_leakage, boot_instability, boot_biology, boot_consistency
        )
        passes = (
            (boot_z["donor_leakage_z"] >= leakage_floor)
            & (boot_z["risk_z"] >= risk_floor)
            & (boot_z["biology_z"] <= biology_ceiling)
        )
        bootstrap_passes += passes.astype(np.int64)
        if write_diagnostics:
            bootstrap_rows.extend(
                {
                    "bootstrap": bootstrap + 1,
                    "gene": gene,
                    "risk_score_raw": float(score),
                    "risk_z": float(rz),
                    "passes_effect_gate": bool(hit),
                }
                for gene, score, rz, hit in zip(
                    genes, boot_raw, boot_z["risk_z"], passes
                )
            )
    bootstrap_pass_fraction = bootstrap_passes / float(n_bootstraps)

    marker_long, marker_summary = _marker_protection(
        matrix, donors, sampled_obs, genes, config
    )

    risk = pd.DataFrame(
        {
            "dataset": config.get("dataset_name", "dataset"),
            "gene": genes,
            "base_hvg_rank": np.arange(1, len(genes) + 1),
            "donor_leakage": leakage,
            "interaction_instability": instability,
            "biology_eta2": biology,
            "protected_biology": biology,
            "cross_donor_consistency": consistency,
            "risk_score_raw": raw_risk,
            "risk_score": raw_risk,
            "permutation_p": permutation_p,
            "permutation_fdr": permutation_fdr,
            "stability_pass_bootstraps": bootstrap_passes,
            "stability_total_bootstraps": n_bootstraps,
            "bootstrap_pass_fraction": bootstrap_pass_fraction,
            **z_values,
        }
    ).merge(marker_summary, on="gene", how="left", validate="one_to_one")

    alpha = float(config.get("permutation_alpha", 0.05))
    bootstrap_min = float(config.get("bootstrap_min_pass_fraction", 0.80))
    significant = risk["permutation_fdr"] <= alpha
    effect_gate = (
        (risk["donor_leakage_z"] >= leakage_floor)
        & (risk["risk_z"] >= risk_floor)
        & (risk["biology_z"] <= biology_ceiling)
    )
    stable = risk["bootstrap_pass_fraction"] >= bootstrap_min
    protected = risk["hard_replicated_marker_protection"].fillna(False).astype(bool)

    risk["passes_permutation_fdr"] = significant
    risk["passes_effect_floors"] = effect_gate
    risk["passes_stability"] = stable
    # Backward compatibility: old public tables called the significance gate maxT.
    risk["passes_maxT"] = significant
    risk["risk_flagged_before_protection"] = significant & effect_gate & stable
    risk["selected_risk_gene"] = risk["risk_flagged_before_protection"] & ~protected
    risk["harmful"] = risk["selected_risk_gene"]
    risk["harmful_gene"] = risk["selected_risk_gene"]

    risk = risk.sort_values(
        ["selected_risk_gene", "risk_score_raw", "gene"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    if write_diagnostics:
        write_csv(sampling, output / f"{prefix}reference_sampling_audit.csv")
        write_csv(pd.DataFrame(permutation_rows), output / f"{prefix}conditional_permutation_long.csv")
        write_csv(pd.DataFrame(bootstrap_rows), output / f"{prefix}stability_bootstrap_long.csv")
        write_csv(marker_long, output / f"{prefix}marker_evidence_by_donor.csv")
        write_csv(marker_summary, output / f"{prefix}replicated_marker_protection.csv")

    audit = {
        "risk_equation": "Z(donor_leakage) + 0.75*Z(interaction_instability) - Z(biology_eta2)",
        "permutation_design": "donor labels permuted within biology strata; per-gene leakage p; BH-FDR",
        "n_permutations": n_permutations,
        "permutation_alpha": alpha,
        "bootstrap_design": "resample within donor x biology strata; recurrence of leakage/risk/biology effect gate",
        "n_bootstraps": n_bootstraps,
        "bootstrap_min_pass_fraction": bootstrap_min,
        "leakage_z_floor": leakage_floor,
        "risk_z_floor": risk_floor,
        "biology_z_ceiling": biology_ceiling,
    }
    return risk, audit


def _control_metadata(
    adata,
    source_mask: np.ndarray,
    genes: list[str],
    retained: list[str],
    config: dict[str, Any],
) -> dict[str, Any]:
    requested = str(config.get("same_n_control_mode", "auto")).lower()
    if requested == "nested":
        requested = "truncate"
    if requested not in {"auto", "independent", "truncate"}:
        raise ValueError("same_n_control_mode must be auto, independent, or truncate")
    if requested == "auto":
        if config.get("independent_same_n_table"):
            actual = "independent"
        elif bool(config.get("selection_is_non_nested", False)):
            actual = "independent"
        elif (
            not config.get("base_gene_table")
            or bool(config.get("allow_internal_independent_refit", False))
        ) and _is_seurat_v3_like(config.get("hvg_method", "seurat_v3")):
            actual = "independent"
        else:
            actual = "truncate"
    else:
        actual = requested

    direct: list[str] | None = None
    if actual == "truncate":
        direct = genes[: len(retained)]
        definition = "top-N truncation of the original ranked panel"
    else:
        direct = _external_independent_panel(config, len(retained))
        if direct is not None:
            definition = "user-supplied independent same-N panel"
        elif not config.get("base_gene_table") or bool(config.get("allow_internal_independent_refit", False)):
            direct = _hvg_ranking(adata, source_mask, config, len(retained))["gene"].astype(str).tolist()
            definition = "independent Scanpy HVG refit at matched N"
        else:
            definition = "independent same-N requested but unavailable"

    return {
        "requested_control_mode": requested,
        "actual_control_mode": actual,
        "direct_control_available": direct is not None,
        "direct_control_definition": definition,
        "direct_control_genes": direct,
    }


def refine_within_domain(config: dict[str, Any]) -> dict[str, Any]:
    require(config, "output_dir", "batch_key", "label_key", "initial_hvg")
    output = resolve_output(config) / "01_refine"
    output.mkdir(parents=True, exist_ok=True)
    adata = load_adata(config)
    count_audit = resolve_config_count_source(adata, config)
    write_csv(count_audit, resolve_output(config) / "raw_count_source_audit.csv")
    reference, query, split = role_masks(adata, config)

    source = str(config.get("hvg_source", "reference")).lower()
    if source not in {"reference", "query", "all"}:
        raise ValueError("hvg_source must be reference, query, or all")
    source_mask = {"reference": reference, "query": query, "all": reference | query}[source]
    ranking, ranking_source = _base_ranking(adata, source_mask, config)
    genes = ranking["gene"].astype(str).tolist()
    write_csv(ranking, output / "base_hvg_ranking.csv")

    risk, risk_audit = _reference_risk_evidence(
        adata,
        reference,
        genes,
        config,
        output=output,
        write_diagnostics=bool(config.get("write_diagnostics", False)),
    )
    write_csv(risk, output / "gene_risk_evidence.csv")
    selected = risk.loc[risk["selected_risk_gene"].astype(bool)].copy()
    write_csv(selected, output / "selected_risk_genes.csv")

    removed = selected["gene"].astype(str).tolist()
    removed_set = set(removed)
    retained = [gene for gene in genes if gene not in removed_set]
    if len(retained) < 2:
        raise RuntimeError("Within-domain gate retained fewer than two genes")

    control = _control_metadata(adata, source_mask, genes, retained, config)
    if bool(config.get("write_controls", False)):
        panels = {"original_hvg": genes, "refined": retained, "rank_tail_same_n": genes[: len(retained)]}
        if control["direct_control_genes"] is not None:
            panels["direct_same_n"] = control["direct_control_genes"]
        long = pd.concat(
            [
                pd.DataFrame({"panel": name, "gene_rank": np.arange(1, len(panel) + 1), "gene": panel})
                for name, panel in panels.items()
            ],
            ignore_index=True,
        )
        manifest = pd.DataFrame(
            [{"panel": name, "n_genes": len(panel), "sha256": gene_hash(panel)} for name, panel in panels.items()]
        )
        write_csv(long, output / "panel_gene_long_table.csv")
        write_csv(manifest, output / "panel_manifest.csv")

    audit = {
        "dataset": config.get("dataset_name", "dataset"),
        "version": VERSION,
        "mode": "within_domain",
        "split": split,
        "hvg_source": source,
        "ranking_source": ranking_source,
        "risk_source": "reference",
        "query_expression_used_for_selection_or_risk": source in {"query", "all"},
        "query_labels_used_for_selection_or_risk": False,
        "initial_hvg": len(genes),
        "removed_genes": len(removed),
        "final_panel_size": len(retained),
        "zero_deletion_allowed": True,
        "harmful_gene_names": removed,
        "final_hvg_genes": retained,
        **risk_audit,
        **{k: v for k, v in control.items() if k != "direct_control_genes"},
    }
    write_json(audit, output / "refinement_audit.json")
    return audit


def _lognorm_selected(full_counts, all_genes: pd.Index, selected_genes: list[str]):
    positions = all_genes.get_indexer(pd.Index(selected_genes, dtype=str))
    if np.any(positions < 0):
        missing = [selected_genes[i] for i, p in enumerate(positions) if p < 0]
        raise KeyError(f"Selected genes missing from count matrix: {missing[:20]}")
    selected = full_counts[:, positions]
    library = np.asarray(full_counts.sum(axis=1)).ravel().astype(np.float64)
    scale = (1e4 / np.maximum(library, 1.0)).astype(np.float32)
    if sparse.issparse(selected):
        x = selected.tocsr().astype(np.float32, copy=True)
        x = x.multiply(scale[:, None]).tocsr()
        x.data = np.log1p(x.data)
        x.eliminate_zeros()
        return x
    x = np.asarray(selected, dtype=np.float32).copy()
    x *= scale[:, None]
    np.log1p(x, out=x)
    return x


def _mean_var(x) -> tuple[np.ndarray, np.ndarray]:
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel().astype(float)
        mean2 = np.asarray(x.multiply(x).mean(axis=0)).ravel().astype(float)
    else:
        x = np.asarray(x, dtype=np.float64)
        mean = x.mean(axis=0)
        mean2 = np.square(x).mean(axis=0)
    variance = np.maximum(mean2 - np.square(mean), 1e-8)
    return mean, variance


def _detection_rate(full_counts, all_genes: pd.Index, selected_genes: list[str]) -> np.ndarray:
    positions = all_genes.get_indexer(pd.Index(selected_genes, dtype=str))
    if np.any(positions < 0):
        raise KeyError("Selected genes missing while calculating detection rates")
    x = full_counts[:, positions]
    if sparse.issparse(x):
        detected = x.copy().tocsr()
        detected.data = np.ones_like(detected.data, dtype=np.float32)
        return np.asarray(detected.mean(axis=0)).ravel().astype(float)
    return np.mean(np.asarray(x) > 0, axis=0).astype(float)


def _reference_biology_protection(reference_risk: pd.DataFrame) -> pd.DataFrame:
    result = reference_risk.copy()
    biology_raw = pd.to_numeric(result["biology_eta2"], errors="coerce").fillna(0).to_numpy(float)
    biology_pct = _pct_rank(biology_raw)
    replication = pd.to_numeric(
        result.get("marker_replication_fraction", pd.Series(0.0, index=result.index)),
        errors="coerce",
    ).fillna(0).clip(0, 1).to_numpy(float)
    hard = result.get(
        "hard_replicated_marker_protection", pd.Series(False, index=result.index)
    ).fillna(False).astype(bool).to_numpy()
    protection = np.clip(0.60 * biology_pct + 0.25 * replication + 0.15 * hard.astype(float), 0, 1)
    result["reference_biology_raw"] = biology_raw
    result["reference_biology_percentile"] = biology_pct
    result["marker_replication_fraction_v3"] = replication
    result["crossdomain_hard_protected"] = hard
    result["biology_protection_score_v3"] = protection
    return result


def refine_cross_domain(config: dict[str, Any]) -> dict[str, Any]:
    require(config, "output_dir", "batch_key", "label_key", "initial_hvg")
    output = resolve_output(config) / "01_refine"
    output.mkdir(parents=True, exist_ok=True)
    adata = load_adata(config)
    count_audit = resolve_config_count_source(adata, config)
    write_csv(count_audit, resolve_output(config) / "raw_count_source_audit.csv")
    reference, query, split = role_masks(adata, config)
    if not reference.any() or not query.any():
        raise ValueError("cross_domain mode requires non-empty Reference and Query")

    # The final cross-domain method refines the Query HVG panel.
    query_config = dict(config)
    query_config["hvg_source"] = "query"
    query_ranking, ranking_source = _base_ranking(adata, query, query_config)
    query_genes = query_ranking["gene"].astype(str).tolist()
    query_rank = {gene: rank for rank, gene in enumerate(query_genes, start=1)}
    write_csv(query_ranking, output / "base_hvg_ranking.csv")
    write_csv(query_ranking, output / "QUERY_hvg_ranking.csv")

    # Reference HVGs are fitted independently at the same nominal budget.
    reference_config = dict(config)
    reference_config.pop("base_gene_table", None)
    reference_config.pop("base_feature_method", None)
    reference_config["hvg_source"] = "reference"
    reference_config["hvg_method"] = str(
        config.get("cross_reference_hvg_method", "seurat_v3")
    )
    reference_ranking = _hvg_ranking(adata, reference, reference_config, len(query_genes))
    reference_genes = reference_ranking["gene"].astype(str).tolist()
    write_csv(reference_ranking, output / "REFERENCE_hvg_ranking.csv")

    reference_risk, reference_audit = _reference_risk_evidence(
        adata,
        reference,
        reference_genes,
        reference_config,
        output=output,
        prefix="REFERENCE_",
        write_diagnostics=bool(config.get("write_diagnostics", False)),
    )
    reference_risk = _reference_biology_protection(reference_risk)
    reference_risk["reference_rule_percentile"] = _pct_rank(
        pd.to_numeric(reference_risk["risk_score_raw"], errors="coerce").fillna(0).to_numpy(float)
    )
    write_csv(reference_risk, output / "REFERENCE_RULE_V1_full_audit_for_crossdomain.csv")

    transportable = reference_risk.loc[reference_risk["gene"].isin(query_rank)].copy()
    if "base_hvg_rank" in transportable.columns:
        transportable = transportable.rename(columns={"base_hvg_rank": "reference_hvg_rank"})
    if transportable.empty:
        raise RuntimeError("No Reference HVGs overlap the Query HVG panel")
    transportable["query_hvg_rank"] = transportable["gene"].map(query_rank).astype(int)
    candidate_genes = transportable["gene"].astype(str).tolist()

    layer = config.get("counts_layer", "counts")
    all_genes = count_var_names(adata, layer).astype(str)
    reference_counts = subset_count_matrix(adata, reference, all_genes.tolist(), layer)
    query_counts = subset_count_matrix(adata, query, all_genes.tolist(), layer)
    if not sparse.issparse(reference_counts):
        reference_counts = np.asarray(reference_counts)
    if not sparse.issparse(query_counts):
        query_counts = np.asarray(query_counts)

    ref_log = _lognorm_selected(reference_counts, all_genes, candidate_genes)
    qry_log = _lognorm_selected(query_counts, all_genes, candidate_genes)
    ref_mean, ref_var = _mean_var(ref_log)
    qry_mean, qry_var = _mean_var(qry_log)
    pooled_sd = np.sqrt(0.5 * (ref_var + qry_var) + 0.05**2)
    mean_shift_effect = np.abs(qry_mean - ref_mean) / np.maximum(pooled_sd, 1e-6)
    ref_detect = _detection_rate(reference_counts, all_genes, candidate_genes)
    qry_detect = _detection_rate(query_counts, all_genes, candidate_genes)
    detection_shift = np.abs(qry_detect - ref_detect)
    mean_shift_pct = _pct_rank(mean_shift_effect)
    detection_shift_pct = _pct_rank(detection_shift)
    dataset_shift_score = 0.70 * mean_shift_pct + 0.30 * detection_shift_pct

    transportable["reference_logmean"] = ref_mean
    transportable["query_logmean"] = qry_mean
    transportable["mean_shift_effect"] = mean_shift_effect
    transportable["mean_shift_percentile"] = mean_shift_pct
    transportable["reference_detection_rate"] = ref_detect
    transportable["query_detection_rate"] = qry_detect
    transportable["detection_rate_shift"] = detection_shift
    transportable["detection_shift_percentile"] = detection_shift_pct
    transportable["dataset_shift_score"] = dataset_shift_score

    r = transportable["reference_rule_percentile"].to_numpy(float)
    s = transportable["dataset_shift_score"].to_numpy(float)
    b = transportable["biology_protection_score_v3"].to_numpy(float)
    transportable["technical_consensus"] = np.minimum(r, s)
    transportable["rule_shift_agreement"] = 1.0 - np.abs(r - s)
    transportable["cross_risk_score"] = (
        transportable["technical_consensus"].to_numpy(float)
        * (0.75 + 0.25 * transportable["rule_shift_agreement"].to_numpy(float))
        * (1.0 - 0.75 * b)
    )

    eligible = transportable.loc[
        ~transportable["crossdomain_hard_protected"].fillna(False).astype(bool)
    ].sort_values(
        [
            "cross_risk_score",
            "technical_consensus",
            "rule_shift_agreement",
            "dataset_shift_score",
            "reference_rule_percentile",
            "query_hvg_rank",
        ],
        ascending=[False, False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    eligible["cross_risk_rank"] = np.arange(1, len(eligible) + 1)
    eligible["cross_risk_percentile"] = eligible["cross_risk_score"].rank(
        method="average", pct=True
    )

    delete_budget = int(config.get("cross_domain_delete_budget", 5))
    if delete_budget < 0:
        raise ValueError("cross_domain_delete_budget must be >= 0")
    if delete_budget > len(eligible):
        raise ValueError(
            f"cross_domain_delete_budget={delete_budget} exceeds {len(eligible)} eligible genes"
        )
    removed = eligible.head(delete_budget)["gene"].astype(str).tolist()
    removed_set = set(removed)
    retained = [gene for gene in query_genes if gene not in removed_set]
    if len(retained) < 2:
        raise RuntimeError("Cross-domain refinement retained fewer than two genes")

    # Full Query-panel evidence table for a stable public output schema.
    evidence = pd.DataFrame(
        {
            "gene": query_genes,
            "base_hvg_rank": np.arange(1, len(query_genes) + 1),
            "query_hvg_rank": np.arange(1, len(query_genes) + 1),
        }
    ).merge(transportable, on=["gene", "query_hvg_rank"], how="left")
    rank_map = dict(zip(eligible["gene"].astype(str), eligible["cross_risk_rank"].astype(int)))
    evidence["cross_risk_rank"] = evidence["gene"].map(rank_map)
    evidence["risk_score_raw"] = pd.to_numeric(
        evidence.get("cross_risk_score"), errors="coerce"
    )
    evidence["risk_score"] = evidence["risk_score_raw"]
    evidence["selected_risk_gene"] = evidence["gene"].isin(removed_set)
    evidence["harmful"] = evidence["selected_risk_gene"]
    evidence["harmful_gene"] = evidence["selected_risk_gene"]
    evidence["passes_maxT"] = evidence["selected_risk_gene"]  # legacy output alias
    evidence["passes_stability"] = True
    evidence["passes_effect_floors"] = True
    evidence["risk_flagged_before_protection"] = evidence["selected_risk_gene"]
    if "hard_replicated_marker_protection" not in evidence:
        evidence["hard_replicated_marker_protection"] = False
    evidence["hard_replicated_marker_protection"] = evidence[
        "hard_replicated_marker_protection"
    ].fillna(False).astype(bool)
    evidence["mode"] = "cross_domain"

    write_csv(transportable, output / "CROSSDOMAIN_reference_query_shift_audit.csv")
    write_csv(eligible, output / "CROSSDOMAIN_RULE_V3_ranking.csv")
    write_csv(evidence, output / "gene_risk_evidence.csv")
    write_csv(
        evidence.loc[evidence["selected_risk_gene"].astype(bool)].copy(),
        output / "selected_risk_genes.csv",
    )

    membership = eligible.head(max(delete_budget, 20 if len(eligible) >= 20 else delete_budget)).copy()
    for budget in (5, 10, 20):
        membership[f"in_k{budget:03d}"] = membership["cross_risk_rank"] <= budget
    write_csv(membership, output / "CROSSDOMAIN_RULE_V3_budget_gene_membership.csv")
    write_csv(
        pd.DataFrame(
            [
                {
                    "mode": "cross_domain",
                    "selection_type": "fixed_budget",
                    "k": delete_budget,
                    "input_n_hvg": len(query_genes),
                    "final_n_hvg": len(retained),
                    "removed_genes": "|".join(removed),
                }
            ]
        ),
        output / "CROSSDOMAIN_RULE_V3_selection_manifest.csv",
    )
    write_csv(
        eligible.head(delete_budget).copy(),
        output / "CROSSDOMAIN_RULE_V3_removal_audit.csv",
    )

    control = _control_metadata(adata, query, query_genes, retained, query_config)
    audit = {
        "dataset": config.get("dataset_name", "dataset"),
        "version": VERSION,
        "mode": "cross_domain",
        "split": split,
        "hvg_source": "query",
        "ranking_source": ranking_source,
        "risk_source": "reference_plus_label_free_reference_query_shift",
        "query_expression_used_for_selection_or_risk": True,
        "query_labels_used_for_selection_or_risk": False,
        "reference_hvg": len(reference_genes),
        "initial_hvg": len(query_genes),
        "cross_domain_delete_budget": delete_budget,
        "primary_supported_budgets": [5, 10, 20],
        "removed_genes": len(removed),
        "final_panel_size": len(retained),
        "zero_deletion_allowed": True,
        "harmful_gene_names": removed,
        "final_hvg_genes": retained,
        "cross_domain_rule": (
            "min(R,S) * (0.75 + 0.25*(1-|R-S|)) * (1 - 0.75*B)"
        ),
        "dataset_shift_rule": "0.70*percentile(mean_shift_effect) + 0.30*percentile(detection_shift)",
        "biology_protection_rule": "0.60*biology_percentile + 0.25*marker_replication + 0.15*hard_protection",
        **reference_audit,
        **{k: v for k, v in control.items() if k != "direct_control_genes"},
    }
    write_json(audit, output / "refinement_audit.json")
    return audit
