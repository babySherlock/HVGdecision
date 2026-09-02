"""AnnData and tabular I/O helpers."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse


def write_csv(frame: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def write_json(value: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    return path


def gene_hash(genes: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(map(str, genes)).encode()).hexdigest()


def matrix_values(matrix: Any, maximum: int = 100_000, seed: int = 20260828) -> np.ndarray:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if values.size > maximum:
        rng = np.random.default_rng(seed)
        values = rng.choice(values, maximum, replace=False)
    return values


def validate_counts(matrix: Any, name: str, tolerance: float = 1e-6) -> dict[str, Any]:
    audit = audit_counts(matrix, name, tolerance=tolerance)
    if not audit["valid"]:
        raise ValueError(str(audit["error"]))
    return {
        key: value
        for key, value in audit.items()
        if key not in {"source", "valid", "error"}
    }


def audit_counts(matrix: Any, name: str, tolerance: float = 1e-6) -> dict[str, Any]:
    """Describe a candidate matrix without hiding why it is not raw counts."""
    values = matrix_values(matrix)
    shape = tuple(int(value) for value in matrix.shape)
    storage = "sparse" if sparse.issparse(matrix) else "dense"
    matrix_dtype = getattr(matrix, "dtype", None)
    if matrix_dtype is None:
        matrix_dtype = np.asarray(matrix).dtype
    dtype = str(matrix_dtype)
    base = {
        "source": name,
        "shape": f"{shape[0]} x {shape[1]}",
        "n_cells": shape[0],
        "n_genes": shape[1],
        "storage": storage,
        "dtype": dtype,
        "sampled_values": int(values.size),
    }
    if values.size == 0:
        return {
            **base,
            "minimum": np.nan,
            "maximum": np.nan,
            "integer_like_fraction": np.nan,
            "valid": False,
            "error": f"{name}: count matrix is empty",
        }
    minimum = float(values.min())
    maximum = float(values.max())
    integer_fraction = float(
        np.mean(np.isclose(values, np.round(values), atol=tolerance, rtol=tolerance))
    )
    statistics = {
        **base,
        "minimum": minimum,
        "maximum": maximum,
        "integer_like_fraction": integer_fraction,
    }
    if float(values.min()) < 0:
        return {
            **statistics,
            "valid": False,
            "error": f"{name}: count matrix contains negative values",
        }
    if integer_fraction < 0.999:
        return {
            **statistics,
            "valid": False,
            "error": (
                f"{name}: counts are not integer-like "
                f"(fraction={integer_fraction:.6f}); do not run "
                "Seurat-v3/scVI on normalized expression"
            ),
        }
    return {
        **statistics,
        "valid": True,
        "error": "",
    }


def is_raw_source(layer: str | None) -> bool:
    return str(layer).lower() in {"raw", "raw.x", "adata.raw.x"}


def count_source_label(layer: str | None) -> str:
    if is_raw_source(layer):
        return "adata.raw.X"
    if layer in (None, "", "X"):
        return "adata.X"
    return f"adata.layers[{layer!r}]"


def looks_like_count_table_path(value: Any) -> bool:
    if not isinstance(value, (str, Path)):
        return False
    text = str(value).lower()
    return text.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz"))


def external_count_matrix(adata, source: Any):
    """Coerce user-supplied counts to cells x genes with explicit gene names.

    Supported inputs are a CSV/TSV path, pandas DataFrame, AnnData/Raw-like
    object, a numeric matrix aligned to ``adata.var_names``, a
    ``(matrix, gene_names)`` pair, or a mapping with ``matrix`` and
    ``gene_names`` entries.
    """
    location = f"user object ({type(source).__name__})"
    payload = source
    supplied_gene_names = None
    supplied_obs_names = None

    if source is adata.X:
        supplied_gene_names = adata.var_names
        supplied_obs_names = adata.obs_names
        location = "user supplied adata.X"
    elif adata.raw is not None and source is adata.raw.X:
        supplied_gene_names = adata.raw.var_names
        supplied_obs_names = adata.obs_names
        location = "user supplied adata.raw.X"
    elif any(source is matrix for matrix in adata.layers.values()):
        matched_layer = next(
            name for name, matrix in adata.layers.items() if source is matrix
        )
        supplied_gene_names = adata.var_names
        supplied_obs_names = adata.obs_names
        location = f"user supplied adata.layers[{matched_layer!r}]"
    elif looks_like_count_table_path(source):
        path = Path(source).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        separator = "\t" if any(
            str(path).lower().endswith(suffix)
            for suffix in (".tsv", ".tsv.gz", ".txt", ".txt.gz")
        ) else ","
        payload = pd.read_csv(path, sep=separator, index_col=0)
        location = str(path)
    elif isinstance(source, dict):
        if "matrix" not in source:
            raise KeyError("External count mapping requires a 'matrix' entry")
        payload = source["matrix"]
        supplied_gene_names = source.get("gene_names", source.get("var_names"))
        supplied_obs_names = source.get("obs_names", source.get("cell_names"))
        location = str(source.get("name", location))
    elif isinstance(source, tuple) and len(source) == 2:
        payload, supplied_gene_names = source
        location = "user tuple (matrix, gene_names)"

    target_obs = pd.Index(adata.obs_names.astype(str))
    target_vars = pd.Index(adata.var_names.astype(str))
    alignment = ""

    if isinstance(payload, pd.DataFrame):
        frame = payload.copy()
        frame.index = frame.index.astype(str)
        frame.columns = frame.columns.astype(str)
        if target_obs.isin(frame.index).all():
            frame = frame.loc[target_obs]
            gene_names = pd.Index(frame.columns, dtype=str)
            alignment = "DataFrame rows matched and reordered by adata.obs_names"
        elif target_obs.isin(frame.columns).all():
            frame = frame.loc[:, target_obs].T
            gene_names = pd.Index(frame.columns, dtype=str)
            alignment = "DataFrame columns matched adata.obs_names; transposed to cells x genes"
        elif frame.shape[0] == adata.n_obs and frame.shape[1] != adata.n_obs:
            gene_names = pd.Index(frame.columns, dtype=str)
            alignment = "DataFrame row order assumed to match adata.obs_names"
        elif frame.shape[1] == adata.n_obs and frame.shape[0] != adata.n_obs:
            frame = frame.T
            gene_names = pd.Index(frame.columns, dtype=str)
            alignment = "DataFrame transposed by dimension; column order assumed to match adata.obs_names"
        else:
            raise ValueError(
                "Cannot align the external count table to adata.obs_names. "
                f"table_shape={frame.shape}, adata_shape={adata.shape}. Supply cell "
                "names on rows/columns or pass {'matrix': ..., 'obs_names': ..., "
                "'gene_names': ...}."
            )
        matrix = frame.to_numpy()
        if not np.issubdtype(matrix.dtype, np.number):
            matrix = frame.apply(pd.to_numeric, errors="raise").to_numpy()
    elif hasattr(payload, "X") and hasattr(payload, "var_names"):
        object_obs = pd.Index(getattr(payload, "obs_names", target_obs).astype(str))
        if not target_obs.isin(object_obs).all():
            missing = target_obs[~target_obs.isin(object_obs)].tolist()
            raise ValueError(
                f"External AnnData/Raw object is missing {len(missing)} cells; "
                f"first 20={missing[:20]}"
            )
        positions = object_obs.get_indexer(target_obs)
        matrix = payload.X[positions, :]
        gene_names = pd.Index(payload.var_names.astype(str))
        alignment = "AnnData/Raw-like object reordered by adata.obs_names"
    else:
        matrix = payload
        if not hasattr(matrix, "shape") or len(matrix.shape) != 2:
            raise TypeError(
                "External counts must be a two-dimensional numeric matrix, DataFrame, "
                "CSV/TSV path, AnnData/Raw-like object, (matrix, gene_names), or mapping"
            )
        if supplied_obs_names is not None:
            external_obs = pd.Index(list(map(str, supplied_obs_names)))
            if len(external_obs) == matrix.shape[0] and target_obs.isin(external_obs).all():
                matrix = matrix[external_obs.get_indexer(target_obs), :]
                alignment = "External matrix rows reordered by supplied obs_names"
            elif len(external_obs) == matrix.shape[1] and target_obs.isin(external_obs).all():
                matrix = matrix[:, external_obs.get_indexer(target_obs)].T
                alignment = "External matrix columns reordered by supplied obs_names and transposed"
            else:
                raise ValueError("supplied obs_names cannot be aligned to adata.obs_names")
        elif matrix.shape[0] == adata.n_obs:
            alignment = "External matrix row order assumed to match adata.obs_names"
        elif matrix.shape[1] == adata.n_obs:
            matrix = matrix.T
            alignment = "External matrix transposed; column order assumed to match adata.obs_names"
        else:
            raise ValueError(
                f"External matrix shape={matrix.shape} has no dimension equal to "
                f"adata.n_obs={adata.n_obs}"
            )
        if supplied_gene_names is None:
            if matrix.shape[1] != adata.n_vars:
                raise ValueError(
                    "External matrix has a different gene dimension from adata. "
                    "Pass (matrix, gene_names) or a mapping containing gene_names."
                )
            gene_names = target_vars
        else:
            gene_names = pd.Index(list(map(str, supplied_gene_names)))

    if matrix.shape[0] != adata.n_obs:
        raise ValueError(
            f"Aligned external counts have {matrix.shape[0]} cells; expected {adata.n_obs}"
        )
    if len(gene_names) != matrix.shape[1]:
        raise ValueError(
            f"gene_names length={len(gene_names)} does not match matrix genes={matrix.shape[1]}"
        )
    if not gene_names.is_unique:
        duplicated = gene_names[gene_names.duplicated()].unique().tolist()
        raise ValueError(f"External count gene names are not unique; first 20={duplicated[:20]}")
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    return matrix, gene_names, location, alignment


def count_matrix(adata, layer: str | None = "counts"):
    """Return raw counts from ``adata.X`` or a named layer."""
    if is_raw_source(layer):
        if adata.raw is None:
            raise KeyError("Requested adata.raw.X, but adata.raw is None")
        return adata.raw.X
    if layer in (None, "", "X"):
        return adata.X
    if layer not in adata.layers:
        raise KeyError(
            f"Missing raw-count layer {layer!r}; available={list(adata.layers)}. "
            "Pass counts_layer=None only when adata.X contains raw integer counts."
        )
    return adata.layers[layer]


def count_var_names(adata, layer: str | None) -> pd.Index:
    """Return the gene index associated with the selected count source."""
    if is_raw_source(layer):
        if adata.raw is None:
            raise KeyError("Requested adata.raw.X, but adata.raw is None")
        return pd.Index(adata.raw.var_names.astype(str))
    return pd.Index(adata.var_names.astype(str))


def subset_count_matrix(adata, obs_indexer, genes: list[str], layer: str | None):
    """Subset the selected raw-count source by cells and genes."""
    if is_raw_source(layer):
        if adata.raw is None:
            raise KeyError("Requested adata.raw.X, but adata.raw is None")
        available = set(adata.raw.var_names.astype(str))
        missing = [gene for gene in genes if gene not in available]
        if missing:
            raise KeyError(
                f"{len(missing)} selected genes are absent from adata.raw.var_names; "
                f"first 20={missing[:20]}"
            )
        return adata.raw[obs_indexer, genes].X
    return count_matrix(adata[obs_indexer, genes], layer)


def resolve_counts_source(
    adata, requested: str | None = "auto"
) -> tuple[str | None, list[dict[str, Any]]]:
    """Audit count candidates and select the first valid source by priority."""
    if requested != "auto":
        label = count_source_label(requested)
        try:
            audit = audit_counts(count_matrix(adata, requested), label)
        except KeyError as error:
            audit = {
                "source": label,
                "valid": False,
                "selected": False,
                "error": str(error),
            }
            return requested, [audit]
        audit["selected"] = bool(audit["valid"])
        return requested, [audit]

    aliases = ["counts", "raw_counts", "raw.counts", "rawcounts", "umi_counts", "count"]
    by_lower = {str(name).lower(): str(name) for name in adata.layers.keys()}
    candidates: list[str | None] = []
    for alias in aliases:
        if alias in by_lower and by_lower[alias] not in candidates:
            candidates.append(by_lower[alias])
    if adata.raw is not None:
        candidates.append("raw")
    candidates.append(None)
    # Also inspect unusually named layers so the user can locate a raw matrix
    # even when its name is not one of the recognized automatic aliases.
    for layer_name in map(str, adata.layers.keys()):
        if layer_name not in candidates:
            candidates.append(layer_name)

    audits: list[dict[str, Any]] = []
    selected: str | None | object = _NO_COUNT_SOURCE
    for candidate in candidates:
        label = count_source_label(candidate)
        try:
            audit = audit_counts(count_matrix(adata, candidate), label)
        except KeyError as error:
            audits.append(
                {"source": label, "valid": False, "selected": False, "error": str(error)}
            )
            continue
        if audit["valid"] and selected is _NO_COUNT_SOURCE:
            selected = candidate
            audit["selected"] = True
        else:
            audit["selected"] = False
        audits.append(audit)
    if selected is not _NO_COUNT_SOURCE:
        return selected, audits
    return None, audits


_NO_COUNT_SOURCE = object()


def resolve_config_count_source(adata, config: dict[str, Any]) -> pd.DataFrame:
    """Resolve/validate an AnnData count source for config/CLI workflows.

    Python ``setup_reference_query`` already resolves ``counts_layer='auto'``.
    This helper gives YAML/CLI workflows the same behavior so ``auto`` is not
    mistaken for a literal layer named ``auto``.
    """
    requested = config.get("counts_layer", "auto")
    source, rows = resolve_counts_source(adata, requested)
    audit = pd.DataFrame(rows)
    selected = (
        not audit.empty
        and "selected" in audit
        and audit["selected"].fillna(False).astype(bool).any()
    )
    if not selected:
        details = " | ".join(
            f"{row.get('source', 'unknown')}: {row.get('error', 'invalid')}"
            for row in rows
        )
        raise ValueError(
            "No valid raw-count source was found for the config workflow. "
            f"requested={requested!r}. Details: {details}"
        )
    config["counts_layer"] = source
    config["counts_source_label"] = count_source_label(source)
    return audit


def load_adata(config: dict[str, Any]):
    import anndata as ad

    if config.get("_adata") is not None:
        adata = config["_adata"]
        if not isinstance(adata, ad.AnnData):
            raise TypeError("config['_adata'] must be an anndata.AnnData object")
    else:
        path = Path(config["input_h5ad"]).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        adata = ad.read_h5ad(path)
    layer = config.get("counts_layer", "counts")
    matrix = count_matrix(adata, layer)
    if sparse.issparse(matrix) and matrix.format != "csr" and not is_raw_source(layer):
        if layer in (None, "", "X"):
            adata.X = matrix.tocsr()
        else:
            adata.layers[layer] = matrix.tocsr()
    if not all(isinstance(value, str) for value in adata.var_names[: min(100, adata.n_vars)]):
        adata.var_names = adata.var_names.astype(str)
    source_var_names = count_var_names(adata, layer)
    if not source_var_names.is_unique:
        raise ValueError(
            f"Gene identifiers in {count_source_label(layer)} must be unique"
        )
    return adata


def role_masks(adata, config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    split_key = config.get("split_key")
    if split_key in (None, "", "null"):
        mask = np.ones(adata.n_obs, dtype=bool)
        return mask, mask.copy(), {"mode": "all_batches", "split_key": None}
    if split_key not in adata.obs:
        raise KeyError(f"split_key={split_key!r} is absent from adata.obs")
    values = adata.obs[split_key].astype(str)
    reference_values = {str(value) for value in config.get("reference_values", [])}
    query_values = {str(value) for value in config.get("query_values", [])}
    if not reference_values or not query_values:
        raise ValueError("reference_values and query_values are required when split_key is used")
    overlap = reference_values & query_values
    if overlap:
        raise ValueError(f"Reference and Query values overlap: {sorted(overlap)}")
    reference = values.isin(reference_values).to_numpy()
    query = values.isin(query_values).to_numpy()
    if not reference.any() or not query.any():
        raise ValueError(
            f"Empty split: reference_cells={reference.sum()}, query_cells={query.sum()}, "
            f"observed_values={sorted(values.unique())}"
        )
    return (
        reference,
        query,
        {
            "mode": "reference_to_query",
            "split_key": split_key,
            "reference_values": sorted(reference_values),
            "query_values": sorted(query_values),
        },
    )


def clean_obs_strings(adata) -> None:
    for column in adata.obs.columns:
        if isinstance(adata.obs[column].dtype, pd.StringDtype):
            adata.obs[column] = adata.obs[column].astype(object)
