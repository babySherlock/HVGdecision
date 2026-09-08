"""Reference -> Query transfer workflow for HVGDecision 0.10.0.

This module deliberately does *not* reimplement HVGDecision's donor-aware risk
engine. Risk discovery is delegated to the installed ``hvgdecision.refine`` so
that the exact v0.10.0 core remains authoritative.

Every independent Seurat-v3 fit starts at span=0.3 and retries 0.5, 0.7, 1.0
only after a numerical LOESS failure. Query selection always starts again from
0.3 and never inherits the span used by Reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_SPANS = (0.3, 0.5, 0.7, 1.0)


def _values(value: str | Sequence[str]) -> list[str]:
    if isinstance(value, str):
        out = [value]
    else:
        out = list(map(str, value))
    out = list(dict.fromkeys(out))
    if not out:
        raise ValueError("Reference/Query values cannot be empty")
    return out


def _mask(adata, batch_key: str, values: Sequence[str]) -> np.ndarray:
    if batch_key not in adata.obs:
        raise KeyError(f"adata.obs is missing {batch_key!r}")
    observed = adata.obs[batch_key].astype(str)
    missing = sorted(set(map(str, values)) - set(observed.unique()))
    if missing:
        raise ValueError(f"Unknown {batch_key} values: {missing}")
    return observed.isin(list(map(str, values))).to_numpy(dtype=bool)


def _resolve_counts(adata, counts: Any):
    """Resolve the same count source selected by HVGDecision when possible."""
    import hvgdecision as hd

    result = counts
    if result is None or result == "auto":
        result = hd.find_raw_counts(adata)

    # Current HVGDecision CountSourceResult can represent internal AnnData
    # sources without storing the matrix itself (matrix=None).
    matrix = getattr(result, "matrix", None)
    source = getattr(result, "source", None)
    location = str(getattr(result, "location", ""))
    gene_names = getattr(result, "gene_names", None)

    if matrix is not None:
        genes = pd.Index(
            gene_names if gene_names is not None else adata.var_names,
            dtype=str,
        )
        return matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix), genes, result

    # Prefer the package's own source resolver. This preserves v0.10.0 count
    # routing when available.
    try:
        from hvgdecision.io import count_matrix, count_var_names

        matrix = count_matrix(adata, source)
        genes = pd.Index(count_var_names(adata, source).astype(str))
        return matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix), genes, result
    except Exception:
        pass

    # Conservative fallback for common internal locations.
    if source is None or location in {"adata.X", "X"} or location.endswith(".X") and "raw" not in location:
        matrix = adata.X
        genes = pd.Index(adata.var_names.astype(str))
    elif source == "raw" or "raw.X" in location:
        if adata.raw is None:
            raise ValueError("Count source points to adata.raw.X, but adata.raw is None")
        matrix = adata.raw.X
        genes = pd.Index(adata.raw.var_names.astype(str))
    elif isinstance(source, str) and source in adata.layers:
        matrix = adata.layers[source]
        genes = pd.Index(adata.var_names.astype(str))
    else:
        raise ValueError(
            "Could not resolve the raw-count matrix from the supplied count source. "
            f"source={source!r}, location={location!r}"
        )

    return matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix), genes, result


def select_seurat_v3_hvgs(
    adata,
    *,
    matrix,
    gene_names: Sequence[str],
    mask: np.ndarray,
    batch_key: str,
    n_hvg: int = 2000,
    spans: Sequence[float] = DEFAULT_SPANS,
    role: str = "selection",
):
    """Independent batch-aware Seurat-v3 HVG fit with auditable span fallback."""
    import anndata as ad
    import scanpy as sc

    mask = np.asarray(mask, dtype=bool)
    genes = pd.Index(map(str, gene_names))
    if mask.sum() < 2:
        raise ValueError(f"{role}: fewer than two cells")

    x = matrix[mask, :]
    obs = adata.obs.loc[mask].copy()
    n_top = min(int(n_hvg), max(2, len(genes) - 1))
    attempts: list[dict[str, Any]] = []
    fit = None
    used_span = None

    for span in dict.fromkeys(float(x) for x in spans):
        attempt = ad.AnnData(
            X=x.copy(),
            obs=pd.DataFrame(
                {batch_key: obs[batch_key].astype(str).to_numpy()},
                index=obs.index.copy(),
            ),
            var=pd.DataFrame(index=genes.copy()),
        )
        kwargs = dict(
            flavor="seurat_v3",
            n_top_genes=n_top,
            span=span,
            subset=False,
            inplace=True,
        )
        if attempt.obs[batch_key].astype(str).nunique() >= 2:
            kwargs["batch_key"] = batch_key
        try:
            sc.pp.highly_variable_genes(attempt, **kwargs)
            attempts.append({"role": role, "span": span, "status": "success", "error": ""})
            fit = attempt
            used_span = span
            break
        except (ValueError, np.linalg.LinAlgError) as error:
            attempts.append(
                {
                    "role": role,
                    "span": span,
                    "status": "failed",
                    "error": f"{type(error).__name__}: {error}",
                }
            )

    audit = pd.DataFrame(attempts)
    if fit is None:
        raise RuntimeError(
            f"{role}: Seurat-v3 failed for all spans {tuple(spans)}\n"
            + audit.to_string(index=False)
        )

    hv = fit.var.copy()
    hv["gene"] = hv.index.astype(str)
    if "highly_variable_rank" in hv.columns:
        hv = hv.sort_values(
            ["highly_variable", "highly_variable_rank"],
            ascending=[False, True],
            kind="mergesort",
        )
    elif "variances_norm" in hv.columns:
        hv = hv.sort_values(
            ["highly_variable", "variances_norm"],
            ascending=[False, False],
            kind="mergesort",
        )
    else:
        hv = hv.sort_values("highly_variable", ascending=False, kind="mergesort")

    selected = (
        hv.loc[hv["highly_variable"].fillna(False), "gene"]
        .astype(str)
        .head(n_top)
        .tolist()
    )
    if len(selected) != n_top:
        raise RuntimeError(f"{role}: expected {n_top} HVGs, obtained {len(selected)}")
    return selected, float(used_span), audit, hv.reset_index(drop=True)


def _build_final_adata(adata, matrix, genes: pd.Index, final_genes: Sequence[str]):
    import anndata as ad

    final_genes = list(map(str, final_genes))
    idx = genes.get_indexer(final_genes)
    if np.any(idx < 0):
        missing = [final_genes[i] for i, x in enumerate(idx) if x < 0]
        raise KeyError(f"Final panel genes missing from count matrix: {missing[:20]}")
    x = matrix[:, idx]
    if sparse.issparse(x):
        x = x.tocsr()

    adata_genes = pd.Index(adata.var_names.astype(str))
    if genes.equals(adata_genes):
        out = adata[:, final_genes].copy()
        out.X = x.copy()
    else:
        out = ad.AnnData(
            X=x.copy(),
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=pd.Index(final_genes, name="gene")),
        )
    out.layers["counts"] = out.X.copy()
    return out


@dataclass
class ReferenceQueryResult:
    reference_result: Any
    reference_hvg_genes: list[str]
    query_hvg_genes: list[str]
    reference_risk_genes: list[str]
    removed_genes: list[str]
    final_hvg_genes: list[str]
    reference_span: float
    query_span: float
    span_audit: pd.DataFrame
    transfer_audit: pd.DataFrame
    adata: Any
    output_dir: Path

    @property
    def final_n_hvg(self) -> int:
        return len(self.final_hvg_genes)

    def __repr__(self) -> str:
        return (
            "ReferenceQueryResult("
            f"reference_risk={len(self.reference_risk_genes)}, "
            f"query_base={len(self.query_hvg_genes)}, "
            f"removed={len(self.removed_genes)}, "
            f"final={self.final_n_hvg}, "
            f"reference_span={self.reference_span}, query_span={self.query_span})"
        )


def refine_reference_query(
    adata,
    *,
    batch_key: str,
    label_key: str,
    reference: str | Sequence[str],
    query: str | Sequence[str],
    counts: Any = "auto",
    n_hvg: int = 2000,
    spans: Sequence[float] = DEFAULT_SPANS,
    output_dir: str | Path = "HVGDecision_reference_query",
    refine_kwargs: dict[str, Any] | None = None,
) -> ReferenceQueryResult:
    """Reference risk -> independently fitted Query HVG -> intersection removal.

    The installed HVGDecision 0.10.0 ``hd.refine`` performs Reference risk
    discovery. This wrapper only supplies the correctly frozen Reference HVGs
    and performs the independent Query HVG fit and transfer step.
    """
    import hvgdecision as hd

    outdir = Path(output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    ref_values = _values(reference)
    qry_values = _values(query)
    overlap = sorted(set(ref_values) & set(qry_values))
    if overlap:
        raise ValueError(f"Reference and Query overlap: {overlap}")

    matrix, genes, count_result = _resolve_counts(adata, counts)
    ref_mask = _mask(adata, batch_key, ref_values)
    qry_mask = _mask(adata, batch_key, qry_values)

    ref_hvgs, ref_span, ref_audit, _ = select_seurat_v3_hvgs(
        adata,
        matrix=matrix,
        gene_names=genes,
        mask=ref_mask,
        batch_key=batch_key,
        n_hvg=n_hvg,
        spans=spans,
        role="reference",
    )
    ref_audit.to_csv(outdir / "REFERENCE_span_fallback_audit.csv", index=False)
    pd.DataFrame(
        {
            "gene": ref_hvgs,
            "reference_hvg_rank": np.arange(1, len(ref_hvgs) + 1),
            "span_used": ref_span,
        }
    ).to_csv(outdir / "REFERENCE_hvg_2000.csv", index=False)

    kwargs = dict(refine_kwargs or {})
    # These arguments are controlled by the transfer wrapper and must not be
    # accidentally overridden through refine_kwargs.
    for reserved in [
        "batch_key", "label_key", "counts", "reference", "hvg_genes",
        "n_hvg", "output_dir", "return_details",
    ]:
        kwargs.pop(reserved, None)

    reference_result = hd.refine(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        counts=count_result,
        reference=ref_values,
        hvg_genes=ref_hvgs,
        output_dir=str(outdir / "01_reference_risk_discovery"),
        return_details=True,
        **kwargs,
    )
    ref_risk = list(map(str, reference_result.removed_genes))
    ref_risk_set = set(ref_risk)
    pd.DataFrame({"gene": ref_risk}).to_csv(outdir / "REFERENCE_risk_genes.csv", index=False)

    # IMPORTANT: independent second fit, starting from 0.3 again.
    qry_hvgs, qry_span, qry_audit, _ = select_seurat_v3_hvgs(
        adata,
        matrix=matrix,
        gene_names=genes,
        mask=qry_mask,
        batch_key=batch_key,
        n_hvg=n_hvg,
        spans=spans,
        role="query",
    )
    qry_audit.to_csv(outdir / "QUERY_span_fallback_audit.csv", index=False)
    pd.DataFrame(
        {
            "gene": qry_hvgs,
            "query_hvg_rank": np.arange(1, len(qry_hvgs) + 1),
            "span_used": qry_span,
        }
    ).to_csv(outdir / "QUERY_hvg_2000.csv", index=False)

    removed = [g for g in qry_hvgs if g in ref_risk_set]
    removed_set = set(removed)
    final_hvgs = [g for g in qry_hvgs if g not in removed_set]
    if len(final_hvgs) < 2:
        raise RuntimeError("Reference risk filtering leaves fewer than two Query HVGs")

    transfer = pd.DataFrame(
        {
            "gene": qry_hvgs,
            "query_hvg_rank": np.arange(1, len(qry_hvgs) + 1),
            "reference_risk_gene": [g in ref_risk_set for g in qry_hvgs],
            "removed_from_query": [g in removed_set for g in qry_hvgs],
        }
    )
    transfer.to_csv(outdir / "REFERENCE_to_QUERY_transfer_audit.csv", index=False)
    pd.DataFrame({"gene": removed}).to_csv(
        outdir / "QUERY_removed_by_REFERENCE_risk.csv", index=False
    )
    pd.DataFrame(
        {"gene": final_hvgs, "final_rank": np.arange(1, len(final_hvgs) + 1)}
    ).to_csv(outdir / "QUERY_refined_hvg_panel.csv", index=False)

    span_audit = pd.concat([ref_audit, qry_audit], ignore_index=True)
    span_audit.to_csv(outdir / "seuratv3_span_fallback_audit.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "reference_cells": int(ref_mask.sum()),
                "reference_donors": len(ref_values),
                "reference_span_used": ref_span,
                "reference_hvg_n": len(ref_hvgs),
                "reference_risk_n": len(ref_risk),
                "query_cells": int(qry_mask.sum()),
                "query_donors": len(qry_values),
                "query_span_used": qry_span,
                "query_hvg_n": len(qry_hvgs),
                "risk_in_query_hvg_n": len(removed),
                "final_query_panel_n": len(final_hvgs),
                "removed_genes": "|".join(removed),
            }
        ]
    )
    summary.to_csv(outdir / "experiment_summary.csv", index=False)

    final_adata = _build_final_adata(adata, matrix, genes, final_hvgs)
    final_adata.uns["hvgdecision_reference_query"] = {
        "reference": ref_values,
        "query": qry_values,
        "reference_span": ref_span,
        "query_span": qry_span,
        "reference_risk_genes": ref_risk,
        "removed_genes": removed,
        "final_n_hvg": len(final_hvgs),
    }

    return ReferenceQueryResult(
        reference_result=reference_result,
        reference_hvg_genes=ref_hvgs,
        query_hvg_genes=qry_hvgs,
        reference_risk_genes=ref_risk,
        removed_genes=removed,
        final_hvg_genes=final_hvgs,
        reference_span=ref_span,
        query_span=qry_span,
        span_audit=span_audit,
        transfer_audit=transfer,
        adata=final_adata,
        output_dir=outdir,
    )
