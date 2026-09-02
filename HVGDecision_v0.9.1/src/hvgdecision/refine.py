"""HVG selection and marker-protected Reference-only risk refinement."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import f as f_distribution
from tqdm.auto import tqdm

from .config import require, resolve_output
from .io import (
    count_var_names,
    gene_hash,
    load_adata,
    role_masks,
    subset_count_matrix,
    write_csv,
    write_json,
)
from .statistics import (
    benjamini_hochberg,
    biology_components,
    compose_risk,
    donor_risk_components,
    eta_squared,
)




def _normalized_method_name(value: str) -> str:
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")


def _is_seurat_v3_like(value: str) -> bool:
    name = _normalized_method_name(value)
    return name in {"seurat_v3", "seuratv3", "scanpy_seurat_v3", "scanpy_seuratv3"} or (
        "seurat" in name and "v3" in name
    )


def _external_independent_panel(config: dict[str, Any], n_genes: int) -> list[str] | None:
    path_value = config.get("independent_same_n_table")
    if not path_value:
        return None
    path = Path(path_value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    gene_column = config.get("independent_gene_column", "gene")
    rank_column = config.get("independent_rank_column", "gene_rank")
    budget_column = config.get("independent_budget_column", "n_hvg")
    if gene_column not in frame:
        raise KeyError(f"Independent same-N table is missing {gene_column!r}")
    if budget_column in frame:
        part = frame.loc[pd.to_numeric(frame[budget_column], errors="coerce").eq(n_genes)].copy()
    else:
        part = frame.copy() if len(frame) == n_genes else frame.iloc[0:0].copy()
    if part.empty:
        return None
    if rank_column in part:
        part[rank_column] = pd.to_numeric(part[rank_column], errors="raise")
        part = part.sort_values([rank_column, gene_column], kind="mergesort")
    genes = part[gene_column].astype(str).drop_duplicates().tolist()
    if len(genes) != n_genes:
        return None
    return genes

def _precomputed_ranking(config: dict[str, Any], available_genes: set[str]) -> pd.DataFrame:
    path = Path(config["base_gene_table"]).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    gene_column = config.get("base_gene_column", "gene")
    rank_column = config.get("base_rank_column", "gene_rank")
    method_column = config.get("base_method_column", "method")
    dataset_column = config.get("base_dataset_column", "dataset")
    missing = [column for column in (gene_column, rank_column) if column not in frame]
    if missing:
        raise KeyError(f"Precomputed feature table is missing columns: {missing}")
    if config.get("base_feature_method") is not None:
        if method_column not in frame:
            raise KeyError(f"base_feature_method was set but {method_column!r} is missing")
        frame = frame.loc[
            frame[method_column].astype(str).eq(str(config["base_feature_method"]))
        ].copy()
    if dataset_column in frame and config.get("dataset_name") is not None:
        matching = frame[dataset_column].astype(str).eq(str(config["dataset_name"]))
        if matching.any():
            frame = frame.loc[matching].copy()
    frame = frame.rename(columns={gene_column: "gene", rank_column: "gene_rank"})
    frame["gene"] = frame["gene"].astype(str)
    frame["gene_rank"] = pd.to_numeric(frame["gene_rank"], errors="raise")
    frame = frame.loc[frame["gene"].isin(available_genes)].copy()
    frame = frame.sort_values(["gene_rank", "gene"], kind="mergesort")
    frame = frame.drop_duplicates("gene", keep="first")
    requested = int(config["initial_hvg"])
    if len(frame) < requested:
        raise RuntimeError(
            f"Precomputed ranking has only {len(frame)} available unique genes; "
            f"{requested} are required"
        )
    frame = frame.head(requested).copy()
    frame["gene_rank"] = np.arange(1, len(frame) + 1)
    return frame.reset_index(drop=True)


def _hvg_ranking(adata, mask: np.ndarray, config: dict[str, Any], n_top: int) -> pd.DataFrame:
    import anndata as ad
    import scanpy as sc

    layer = config.get("counts_layer", "counts")
    genes = count_var_names(adata, layer).astype(str).tolist()
    counts = subset_count_matrix(adata, mask, genes, layer).copy()
    work = ad.AnnData(
        X=counts,
        obs=adata.obs.loc[mask].copy(),
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    flavor = str(config.get("hvg_method", "seurat_v3"))
    requested = min(int(n_top), work.n_vars - 1)
    kwargs: dict[str, Any] = {
        "flavor": flavor,
        "n_top_genes": requested,
        "inplace": True,
    }
    batch_key = config["batch_key"]
    if work.obs[batch_key].astype(str).nunique() >= 2:
        kwargs["batch_key"] = batch_key

    # Seurat-v3 expects raw counts.  Classic Seurat/Cell Ranger flavors expect
    # log-normalized expression, so normalize only the temporary fitting object.
    if flavor in {"seurat", "cell_ranger"}:
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)

    if _is_seurat_v3_like(flavor):
        spans = [float(config.get("seurat_v3_span", 0.3)), 0.5, 0.7, 1.0]
        failures = []
        for span in dict.fromkeys(spans):
            try:
                attempt = work.copy()
                sc.pp.highly_variable_genes(attempt, span=span, **kwargs)
                work = attempt
                break
            except (ValueError, np.linalg.LinAlgError) as error:
                failures.append(f"span={span}: {error!r}")
        else:
            raise RuntimeError(
                "Seurat-v3 HVG failed for all span values: " + " | ".join(failures)
            )
    else:
        sc.pp.highly_variable_genes(work, **kwargs)

    frame = work.var.reset_index(names="gene").copy()
    frame["gene"] = frame["gene"].astype(str)
    frame = frame.loc[frame["highly_variable"].fillna(False)].copy()
    if "highly_variable_rank" in frame and frame["highly_variable_rank"].notna().any():
        frame = frame.sort_values(["highly_variable_rank", "gene"], kind="mergesort")
    elif "variances_norm" in frame:
        frame = frame.sort_values(["variances_norm", "gene"], ascending=[False, True], kind="mergesort")
    elif "dispersions_norm" in frame:
        frame = frame.sort_values(["dispersions_norm", "gene"], ascending=[False, True], kind="mergesort")
    else:
        frame = frame.sort_values("gene", kind="mergesort")
    frame = frame.head(requested).copy()
    frame["gene_rank"] = np.arange(1, len(frame) + 1)
    if len(frame) != requested or frame["gene"].nunique() != requested:
        raise RuntimeError(
            f"Expected exactly {requested} unique HVGs, observed rows={len(frame)}, "
            f"unique={frame['gene'].nunique()}"
        )
    return frame.reset_index(drop=True)


def _composite_biology_labels(obs: pd.DataFrame, biology_keys: list[str]) -> np.ndarray:
    """Build categorical strata that preserve all declared biology axes."""
    if len(biology_keys) == 1:
        return obs[biology_keys[0]].astype(str).to_numpy()
    return (
        obs[biology_keys]
        .astype(str)
        .agg(" | ".join, axis=1)
        .to_numpy(dtype=object)
    )


def _sample_reference(adata, reference: np.ndarray, genes: list[str], config: dict[str, Any]):
    layer = config.get("counts_layer", "counts")
    batch_key = config["batch_key"]
    label_key = config["label_key"]
    biology_keys = list(dict.fromkeys(config.get("biology_keys", [label_key])))
    maximum = int(config.get("max_cells_per_batch_celltype", 200))
    rng = np.random.default_rng(int(config.get("selection_seed", 20260828)))
    meta = adata.obs.loc[reference, [batch_key, *biology_keys]].copy()
    meta["_biology_stratum"] = _composite_biology_labels(meta, biology_keys)
    reference_indices = np.flatnonzero(reference)
    selected = []
    rows = []
    for (batch, stratum), group in meta.groupby(
        [batch_key, "_biology_stratum"], observed=True, sort=True
    ):
        local = meta.index.get_indexer(group.index)
        available = reference_indices[local]
        take = min(len(available), maximum)
        chosen = np.sort(rng.choice(available, take, replace=False))
        selected.extend(chosen.tolist())
        rows.append(
            {
                "batch": str(batch),
                "biology_stratum": str(stratum),
                "available": len(available),
                "sampled": take,
            }
        )
    selected = np.asarray(sorted(set(selected)), dtype=int)
    counts = subset_count_matrix(adata, selected, genes, layer)
    counts = (
        counts.tocsr().astype(np.float32)
        if sparse.issparse(counts)
        else sparse.csr_matrix(counts, dtype=np.float32)
    )
    library = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=np.float32), where=library > 0)
    logged = (sparse.diags(scale) @ counts).tocsr()
    logged.data = np.log1p(logged.data)
    obs = adata.obs.iloc[selected].copy()
    labels = _composite_biology_labels(obs, biology_keys)
    return (
        logged.toarray().astype(np.float32, copy=False),
        obs[batch_key].astype(str).to_numpy(),
        labels,
        obs,
        pd.DataFrame(rows),
    )


def _marker_protection_single_axis(
    matrix: np.ndarray,
    donors: np.ndarray,
    labels: np.ndarray,
    genes: list[str],
    config: dict[str, Any],
    biology_key: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    fdr = float(config.get("marker_fdr", 0.01))
    effect_floor = float(config.get("marker_min_log_effect", 0.50))
    min_cells_per_side = int(config.get("marker_min_cells_per_side", 10))
    donor_tables = []
    for donor in sorted(np.unique(donors)):
        positions = np.flatnonzero(donors == donor)
        donor_matrix = matrix[positions]
        donor_labels = labels[positions]
        categories = sorted(np.unique(donor_labels))
        if len(categories) < 2:
            continue
        eta = eta_squared(donor_matrix, donor_labels)
        n_cells, n_groups = donor_matrix.shape[0], len(categories)
        numerator = eta / max(n_groups - 1, 1)
        denominator = (1.0 - eta) / max(n_cells - n_groups, 1)
        statistic = np.divide(
            numerator, denominator, out=np.full_like(numerator, np.inf), where=denominator > 1e-15
        )
        p_value = f_distribution.sf(statistic, max(n_groups - 1, 1), max(n_cells - n_groups, 1))
        q_value = benjamini_hochberg(np.nan_to_num(p_value, nan=1.0))
        means = np.vstack(
            [donor_matrix[donor_labels == category].mean(axis=0) for category in categories]
        )
        top_index = np.argmax(means, axis=0)
        top_label = np.asarray(categories, dtype=object)[top_index]
        effects = np.empty(len(genes), dtype=float)
        eligible = np.zeros(len(genes), dtype=bool)
        for gene_index, category_index in enumerate(top_index):
            top_mask = donor_labels == categories[category_index]
            rest_mask = ~top_mask
            rest = donor_matrix[rest_mask, gene_index]
            effects[gene_index] = means[category_index, gene_index] - float(rest.mean())
            eligible[gene_index] = (
                int(top_mask.sum()) >= min_cells_per_side
                and int(rest_mask.sum()) >= min_cells_per_side
            )
        donor_tables.append(
            pd.DataFrame(
                {
                    "gene": genes,
                    "donor": donor,
                    "biology_key": biology_key,
                    "marker_q_value": q_value,
                    "top_enriched_group": top_label,
                    "top_vs_rest_log_effect": effects,
                    "marker_donor_eligible": eligible,
                    "donor_marker_significant": (
                        eligible & (q_value <= fdr) & (effects >= effect_floor)
                    ),
                }
            )
        )
    long = pd.concat(donor_tables, ignore_index=True) if donor_tables else pd.DataFrame()
    minimum_donors = int(
        config.get("marker_min_eligible_donors", config.get("marker_min_donors", 2))
    )
    minimum_fraction = float(
        config.get("marker_replication_fraction", config.get("marker_min_donor_fraction", 0.80))
    )
    records = []
    for gene in genes:
        part = long.loc[long["gene"].eq(gene)] if not long.empty else pd.DataFrame()
        eligible_part = (
            part.loc[part["marker_donor_eligible"].astype(bool)]
            if "marker_donor_eligible" in part
            else part
        )
        significant = (
            eligible_part.loc[eligible_part["donor_marker_significant"].astype(bool)]
            if "donor_marker_significant" in eligible_part
            else eligible_part.iloc[0:0]
        )
        labels_seen = (
            significant.get("top_enriched_group", pd.Series(dtype=str)).astype(str).tolist()
        )
        mode_label, mode_count = Counter(labels_seen).most_common(1)[0] if labels_seen else ("", 0)
        n_eligible = len(eligible_part)
        required = max(
            minimum_donors,
            int(np.ceil(minimum_fraction * max(n_eligible, 1))),
        )
        records.append(
            {
                "gene": gene,
                "biology_key": biology_key,
                "n_reference_donors": len(part),
                "n_eligible_marker_donors": n_eligible,
                "n_significant_marker_donors": len(significant),
                "replicated_marker_group": mode_label,
                "n_donors_same_marker_group": mode_count,
                "required_marker_donors": required,
                "axis_marker_replication_fraction": (
                    float(mode_count / n_eligible) if n_eligible else 0.0
                ),
                "axis_marker_protection": bool(
                    n_eligible >= minimum_donors and mode_count >= required
                ),
            }
        )
    return long, pd.DataFrame(records)


def _marker_protection(
    matrix: np.ndarray,
    donors: np.ndarray,
    sampled_obs: pd.DataFrame,
    genes: list[str],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Protect replicated markers on any declared biological axis.

    ``biology_keys`` are biological variables whose signal must not be treated
    as batch contamination (for example cell type, treatment, or CD45 status).
    ``protected_genes`` are an explicit hard-prior escape hatch for known
    defining genes such as PTPRC when the reference design is sparse.
    """
    biology_keys = list(dict.fromkeys(config.get("biology_keys", [config["label_key"]])))
    long_parts = []
    summary_parts = []
    for key in biology_keys:
        long, summary = _marker_protection_single_axis(
            matrix, donors, sampled_obs[key].astype(str).to_numpy(), genes, config, key
        )
        if not long.empty:
            long_parts.append(long)
        summary_parts.append(summary)
    long = pd.concat(long_parts, ignore_index=True) if long_parts else pd.DataFrame()
    stacked = pd.concat(summary_parts, ignore_index=True) if summary_parts else pd.DataFrame()
    protected_set = set(str(gene) for gene in config.get("protected_genes", []))
    records = []
    for gene in genes:
        part = stacked.loc[stacked["gene"].eq(gene)] if not stacked.empty else pd.DataFrame()
        protected_axes = (
            part.loc[part["axis_marker_protection"].astype(bool), "biology_key"]
            .astype(str)
            .tolist()
            if "axis_marker_protection" in part
            else []
        )
        explicit = gene in protected_set
        replication_fraction = (
            float(pd.to_numeric(part.get("axis_marker_replication_fraction", 0), errors="coerce").fillna(0).max())
            if not part.empty
            else 0.0
        )
        records.append(
            {
                "gene": gene,
                "protected_by_biology_keys": ";".join(protected_axes),
                "explicit_protected_gene": explicit,
                "marker_replication_fraction": replication_fraction,
                "marker_protected": bool(protected_axes),
                "explicitly_protected": bool(explicit),
                "hard_protected": bool(protected_axes or explicit),
                "hard_replicated_marker_protection": bool(protected_axes or explicit),
            }
        )
    return long, pd.DataFrame(records)


def _stability(
    matrix: np.ndarray,
    donors: np.ndarray,
    labels: np.ndarray,
    threshold: float,
    config: dict[str, Any],
) -> tuple[np.ndarray, pd.DataFrame]:
    iterations = int(config.get("n_bootstraps", 50))
    fraction = float(config.get("bootstrap_fraction", 0.80))
    minimum = int(config.get("min_group_cells", 10))
    rng = np.random.default_rng(int(config.get("selection_seed", 20260828)) + 2)
    strata = {}
    for donor in sorted(np.unique(donors)):
        for label in sorted(np.unique(labels)):
            positions = np.flatnonzero((donors == donor) & (labels == label))
            if positions.size:
                strata[(donor, label)] = positions
    passes = np.zeros(matrix.shape[1], dtype=int)
    rows = []
    for iteration in tqdm(range(iterations), desc="Reference stability bootstrap"):
        sample = []
        for positions in strata.values():
            take = min(len(positions), max(minimum, int(np.ceil(fraction * len(positions)))))
            sample.extend(rng.choice(positions, take, replace=False).tolist())
        sample = np.asarray(sorted(sample), dtype=int)
        leakage, instability = donor_risk_components(
            matrix[sample], donors[sample], labels[sample], minimum
        )
        biology, consistency = biology_components(
            matrix[sample], donors[sample], labels[sample], minimum
        )
        score, _ = compose_risk(leakage, instability, biology, consistency)
        flag = score > threshold
        passes += flag.astype(int)
        rows.extend(
            {
                "bootstrap": iteration + 1,
                "gene_index": index,
                "risk_score": float(value),
                "passes": bool(hit),
            }
            for index, (value, hit) in enumerate(zip(score, flag))
        )
    return passes, pd.DataFrame(rows)


def _refine_panel_v08_legacy(config: dict[str, Any]) -> dict[str, Any]:
    require(config, "output_dir", "batch_key", "label_key", "initial_hvg")
    if not config.get("input_h5ad") and config.get("_adata") is None:
        raise KeyError("Either input_h5ad or an in-memory AnnData object is required")
    output = resolve_output(config) / "01_refine"
    output.mkdir(parents=True, exist_ok=True)
    adata = load_adata(config)
    reference, query, split = role_masks(adata, config)
    source = str(config.get("hvg_source", "reference")).lower()
    source_mask = {"reference": reference, "query": query, "all": reference | query}.get(source)
    if source_mask is None:
        raise ValueError("hvg_source must be reference, query, or all")
    if config.get("base_gene_table"):
        ranking = _precomputed_ranking(
            config,
            set(count_var_names(adata, config.get("counts_layer", "counts"))),
        )
        ranking_source = "precomputed_table"
    else:
        ranking = _hvg_ranking(adata, source_mask, config, int(config["initial_hvg"]))
        ranking_source = "scanpy_fit"
    genes = ranking["gene"].astype(str).tolist()
    write_csv(ranking, output / "base_hvg_ranking.csv")
    matrix, donors, labels, sampled_obs, sampling = _sample_reference(adata, reference, genes, config)
    write_diagnostics = bool(config.get("write_diagnostics", False))
    if write_diagnostics:
        write_csv(sampling, output / "reference_sampling_audit.csv")
    if len(np.unique(donors)) < 2:
        raise ValueError("Risk estimation requires at least two Reference batches/donors")
    minimum = int(config.get("min_group_cells", 10))
    leakage, instability = donor_risk_components(matrix, donors, labels, minimum)
    biology, consistency = biology_components(matrix, donors, labels, minimum)
    observed_score, z_values = compose_risk(leakage, instability, biology, consistency)
    n_permutations = int(config.get("n_permutations", 500))
    alpha = float(config.get("permutation_alpha", 0.05))
    rng = np.random.default_rng(int(config.get("selection_seed", 20260828)) + 1)
    null_maxima = np.empty(n_permutations, dtype=float)
    for index in tqdm(range(n_permutations), desc="Reference conditional maxT"):
        permuted = donors.copy()
        for label in sorted(np.unique(labels)):
            positions = np.flatnonzero(labels == label)
            permuted[positions] = rng.permutation(permuted[positions])
        null_leakage, null_instability = donor_risk_components(matrix, permuted, labels, minimum)
        null_score, _ = compose_risk(null_leakage, null_instability, biology, consistency)
        null_maxima[index] = float(np.nanmax(null_score))
    threshold = float(np.quantile(null_maxima, 1.0 - alpha))
    if write_diagnostics:
        write_csv(
            pd.DataFrame(
                {"permutation": np.arange(1, n_permutations + 1), "maximum_null_risk": null_maxima}
            ),
            output / "conditional_maxT_null.csv",
        )
    stability_passes, stability_long = _stability(matrix, donors, labels, threshold, config)
    stability_long["gene"] = [genes[index] for index in stability_long["gene_index"]]
    if write_diagnostics:
        write_csv(stability_long.drop(columns="gene_index"), output / "stability_bootstrap_long.csv")
    marker_long, marker_summary = _marker_protection(matrix, donors, sampled_obs, genes, config)
    if write_diagnostics:
        write_csv(marker_long, output / "marker_evidence_by_donor.csv")
        write_csv(marker_summary, output / "replicated_marker_protection.csv")
    risk = pd.DataFrame(
        {
            "dataset": config.get("dataset_name", "dataset"),
            "gene": genes,
            "base_hvg_rank": np.arange(1, len(genes) + 1),
            "donor_leakage": leakage,
            "interaction_instability": instability,
            "protected_biology": biology,
            "cross_donor_consistency": consistency,
            "risk_score": observed_score,
            "empirical_maxT_threshold": threshold,
            "empirical_fwer_p": [
                (1 + np.sum(null_maxima >= score)) / (n_permutations + 1)
                for score in observed_score
            ],
            "stability_pass_bootstraps": stability_passes,
            "stability_total_bootstraps": int(config.get("n_bootstraps", 50)),
            **z_values,
        }
    ).merge(marker_summary, on="gene", how="left", validate="one_to_one")
    risk["passes_maxT"] = risk["risk_score"] > threshold
    risk["passes_stability"] = risk["stability_pass_bootstraps"] >= np.ceil(
        float(config.get("bootstrap_min_pass_fraction", 0.80)) * int(config.get("n_bootstraps", 50))
    )
    risk["passes_effect_floors"] = (
        risk["donor_leakage_z"] >= float(config.get("leakage_z_floor", 1.0))
    ) & (risk["protected_biology_z"] <= float(config.get("biology_z_ceiling", 0.0)))
    risk["selected_risk_gene"] = (
        risk["passes_maxT"]
        & risk["passes_stability"]
        & risk["passes_effect_floors"]
        & ~risk["hard_replicated_marker_protection"].fillna(False)
    )
    risk = risk.sort_values(
        ["selected_risk_gene", "risk_score", "gene"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    write_csv(risk, output / "gene_risk_evidence.csv")
    removed = risk.loc[risk["selected_risk_gene"], "gene"].astype(str).tolist()
    retained = [gene for gene in genes if gene not in set(removed)]
    if len(retained) < 2:
        raise RuntimeError(
            f"Risk gate retained only {len(retained)} genes; refusing an unusable panel"
        )
    write_csv(risk.loc[risk["selected_risk_gene"]].copy(), output / "selected_risk_genes.csv")
    requested_control_mode = str(config.get("same_n_control_mode", "auto")).lower()
    if requested_control_mode == "nested":
        requested_control_mode = "truncate"
    if requested_control_mode not in {"auto", "independent", "truncate"}:
        raise ValueError("same_n_control_mode must be auto, independent, or truncate")
    if requested_control_mode == "auto":
        if config.get("independent_same_n_table"):
            actual_control_mode = "independent"
        elif bool(config.get("selection_is_non_nested", False)):
            actual_control_mode = "independent"
        elif (
            not config.get("base_gene_table")
            or bool(config.get("allow_internal_independent_refit", False))
        ) and _is_seurat_v3_like(config.get("hvg_method", "seurat_v3")):
            actual_control_mode = "independent"
        else:
            actual_control_mode = "truncate"
    else:
        actual_control_mode = requested_control_mode

    direct: list[str] | None
    if actual_control_mode == "truncate":
        direct = genes[: len(retained)]
        direct_control_definition = "top-N truncation of the original ranked panel"
    else:
        direct = _external_independent_panel(config, len(retained))
        if direct is not None:
            direct_control_definition = "user-supplied independent same-N panel"
        elif (
            not config.get("base_gene_table")
            or bool(config.get("allow_internal_independent_refit", False))
        ):
            direct = (
                _hvg_ranking(adata, source_mask, config, len(retained))["gene"]
                .astype(str)
                .tolist()
            )
            direct_control_definition = "independent Scanpy Seurat-v3 refit at matched N"
        else:
            direct_control_definition = (
                "independent same-N requested but unavailable for the external method; "
                "no surrogate panel was fabricated"
            )

    # Matched control panels were useful during method-development benchmarks,
    # but are not part of the compact public output in v0.8. They can still be
    # written explicitly for method diagnostics.
    if bool(config.get("write_controls", False)):
        panels: dict[str, list[str]] = {
            "original_hvg": genes,
            "refined": retained,
            "rank_tail_same_n": genes[: len(retained)],
        }
        if direct is not None:
            panels["direct_same_n"] = direct
        long = pd.concat(
            [
                pd.DataFrame(
                    {
                        "panel": name,
                        "gene_rank": np.arange(1, len(panel_genes) + 1),
                        "gene": panel_genes,
                    }
                )
                for name, panel_genes in panels.items()
            ],
            ignore_index=True,
        )
        manifest = pd.DataFrame(
            [
                {"panel": name, "n_genes": len(panel_genes), "sha256": gene_hash(panel_genes)}
                for name, panel_genes in panels.items()
            ]
        )
        write_csv(long, output / "panel_gene_long_table.csv")
        write_csv(manifest, output / "panel_manifest.csv")
    audit = {
        "dataset": config.get("dataset_name", "dataset"),
        "version": "0.8.0",
        "split": split,
        "hvg_source": source,
        "ranking_source": ranking_source,
        "risk_source": "reference",
        "biology_keys": list(config.get("biology_keys", [config.get("label_key")])),
        "protected_genes": list(config.get("protected_genes", [])),
        "risk_conditioning": "conditional within joint biology strata",
        "marker_protection": "replicated markers on any biology_key plus explicit protected_genes",
        "query_expression_used_for_hvg_selection": source in {"query", "all"},
        "query_labels_used_for_selection_or_risk": False,
        "initial_hvg": len(genes),
        "removed_genes": len(removed),
        "final_panel_size": len(retained),
        "requested_control_mode": requested_control_mode,
        "actual_control_mode": actual_control_mode,
        "direct_control_available": direct is not None,
        "direct_control_definition": direct_control_definition,
        "zero_deletion_allowed": True,
        "empirical_maxT_threshold": threshold,
        "reference_batches": sorted(np.unique(donors).tolist()),
        "stability_design": "stratified_cell_bootstrap"
        if len(np.unique(donors)) == 2
        else "stratified_cell_bootstrap; independent donor-split validation recommended",
    }
    write_json(audit, output / "audit.json")
    report = f"""# HVGDecision panel-refinement report

- Dataset: {audit["dataset"]}
- HVG source: {source}
- Risk-estimation source: Reference only
- Initial panel: {len(genes)} genes
- Empirical maxT threshold: {threshold:.4f}
- Automatically removed: {len(removed)} genes
- Final refined panel: {len(retained)} genes
- Query labels used before final evaluation: no

Zero deletion was permitted. Removal required maxT significance, bootstrap
stability, effect-size gates, and absence of replicated-marker / explicit
biological protection. Declared biology_keys are conditioned out of donor-risk
assessment rather than being treated as batch contamination.
The final HVG panel is intended for an independent downstream integration workflow.
"""
    (output / "decision_report.md").write_text(report, encoding="utf-8")
    return audit


def refine_panel(config: dict[str, Any]) -> dict[str, Any]:
    """Dispatch refinement to the explicitly selected HVGDecision mode.

    ``mode`` is mandatory in v0.9. Use ``within_domain`` for multi-donor /
    multi-batch refinement inside one experimental domain and ``cross_domain``
    for cross-dataset or cross-technology Reference→Query refinement.
    """
    from .modes import normalize_mode, refine_cross_domain, refine_within_domain

    mode = normalize_mode(config.get("mode"))
    routed = dict(config)
    routed["mode"] = mode
    if mode == "within_domain":
        return refine_within_domain(routed)
    return refine_cross_domain(routed)
