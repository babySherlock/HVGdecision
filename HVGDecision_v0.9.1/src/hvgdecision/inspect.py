"""Input inspection command."""

from __future__ import annotations

from typing import Any

from .config import require, resolve_output
from .io import (
    count_matrix,
    load_adata,
    resolve_config_count_source,
    role_masks,
    validate_counts,
    write_csv,
    write_json,
)
from .modes import normalize_mode


def inspect_dataset(config: dict[str, Any]) -> dict[str, Any]:
    require(config, "batch_key", "label_key", "output_dir")
    if not config.get("input_h5ad") and config.get("_adata") is None:
        raise KeyError("Either input_h5ad or an in-memory AnnData object is required")
    output = resolve_output(config) / "00_inspect"
    output.mkdir(parents=True, exist_ok=True)
    adata = load_adata(config)
    source_audit = resolve_config_count_source(adata, config)
    write_csv(source_audit, resolve_output(config) / "raw_count_source_audit.csv")
    batch_key = config["batch_key"]
    label_key = config["label_key"]
    missing = [key for key in (batch_key, label_key) if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing adata.obs columns: {missing}")
    reference, query, split = role_masks(adata, config)
    count_audit = validate_counts(
        count_matrix(adata, config.get("counts_layer", "counts")), "input"
    )
    reference_labels = set(
        adata.obs.loc[reference, label_key].dropna().astype(str)
    )
    query_labels = set(
        adata.obs.loc[query, label_key].dropna().astype(str)
    )
    shared_labels = sorted(reference_labels & query_labels)
    batch_label = (
        adata.obs.assign(_role="other")
        .assign(
            _role=lambda frame: frame["_role"].mask(reference, "reference").mask(query, "query")
        )
        .groupby(["_role", batch_key, label_key], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    write_csv(batch_label, output / "batch_celltype_counts.csv")
    mode = normalize_mode(config.get("mode")) if config.get("mode") else None
    audit = {
        "dataset": config.get("dataset_name", "dataset"),
        "mode": mode,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_batches": int(adata.obs[batch_key].astype(str).nunique()),
        "n_celltypes": int(adata.obs[label_key].astype(str).nunique()),
        "reference_cells": int(reference.sum()),
        "query_cells": int(query.sum()),
        "reference_celltypes": len(reference_labels),
        "query_celltypes": len(query_labels),
        "shared_celltypes": len(shared_labels),
        "shared_celltype_names": shared_labels,
        "counts": count_audit,
        "split": split,
        "hvg_source": config.get(
            "hvg_source", "query" if mode == "cross_domain" else "reference"
        ),
        "query_expression_used_for_hvg_selection": bool(mode == "cross_domain"),
        "query_labels_used_for_risk_estimation": False,
    }
    write_json(audit, output / "audit.json")
    return audit
