"""Metadata-only replication adequacy and design-identifiability audit."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components


@dataclass
class DesignAuditResult:
    summary: dict
    coverage: pd.DataFrame
    by_celltype: pd.DataFrame
    by_donor: pd.DataFrame

    @property
    def audit(self) -> pd.DataFrame:
        return pd.DataFrame([self.summary])


def selection_mask(adata, batch_key: str, label_key: str, reference=None) -> np.ndarray:
    if batch_key == label_key:
        raise ValueError("batch_key and label_key must be distinct columns")

    for key in (batch_key, label_key):
        if key not in adata.obs:
            raise KeyError(f"Missing adata.obs column: {key!r}")

    if adata.n_obs == 0:
        raise ValueError("AnnData contains no cells")

    if reference is None:
        mask = np.ones(adata.n_obs, dtype=bool)
    else:
        if isinstance(reference, str):
            raise TypeError("reference must be a list of donor IDs, not a string")
        reference = list(map(str, reference))
        if not reference or len(reference) != len(set(reference)):
            raise ValueError("reference must contain distinct donor IDs")
        observed = set(adata.obs[batch_key].astype(str))
        missing = set(reference) - observed
        if missing:
            raise ValueError(f"Reference donors not found: {sorted(missing)}")
        mask = adata.obs[batch_key].astype(str).isin(reference).to_numpy()

    obs = adata.obs.loc[mask, [batch_key, label_key]]
    for key in (batch_key, label_key):
        values = obs[key]
        if values.isna().any():
            raise ValueError(f"{key!r} contains missing metadata in the discovery cohort")
        if values.astype(str).str.strip().eq("").any():
            raise ValueError(f"{key!r} contains empty metadata in the discovery cohort")

    return mask


def audit_metadata(adata, *, batch_key: str, label_key: str, reference=None, tau: float = 15.0):
    """Compute weighted replication adequacy and supported-design identifiability.

    These quantities are descriptive audits only; they do not enter the
    gene-level dual-evidence fusion equation.
    """

    mask = selection_mask(adata, batch_key, label_key, reference)
    obs = adata.obs.loc[mask, [batch_key, label_key]].astype(str)

    donors = sorted(obs[batch_key].unique())
    labels = sorted(obs[label_key].unique())
    cross = pd.crosstab(obs[batch_key], obs[label_key]).reindex(
        index=donors, columns=labels, fill_value=0
    )

    n = cross.to_numpy(dtype=float)
    d, k = n.shape
    support = n / (n + float(tau))
    effective_by_type = support.sum(axis=0)
    type_weight = np.minimum(1.0, n.sum(axis=0) / (d * float(tau)))
    effective_donor_support = float(
        type_weight.sum() / np.sum(type_weight / (effective_by_type + 1e-12))
    )

    donor_mean_support = support.mean(axis=1)
    donor_balance = float(
        np.clip(
            donor_mean_support.sum() ** 2
            / (d * np.square(donor_mean_support).sum()),
            0,
            1,
        )
    )
    replication_adequacy = effective_donor_support * np.sqrt(donor_balance)

    edges = n >= float(tau)
    active_d = edges.any(axis=1)
    active_k = edges.any(axis=0)
    supported = edges[np.ix_(active_d, active_k)]
    ds, ks = supported.shape
    expected_rank = max(ds + ks - 1, 0)

    if supported.any():
        graph = sparse.bmat(
            [
                [None, sparse.csr_matrix(supported)],
                [sparse.csr_matrix(supported.T), None],
            ],
            format="csr",
        )
        components = int(
            connected_components(graph, directed=False, return_labels=False)
        )
        design_rank = ds + ks - components
    else:
        components = 0
        design_rank = 0

    identifiability = float(design_rank / expected_rank) if expected_rank else 0.0
    identifiable = bool(
        ds >= 2
        and ks >= 2
        and components == 1
        and design_rank == expected_rank
    )

    summary = {
        "selected_mode": "dual_evidence",
        "reason": "dataset_structure_audit_only",
        "n_cells": int(len(obs)),
        "n_donors": int(d),
        "n_celltypes": int(k),
        "n_supported_donors": int(ds),
        "n_supported_celltypes": int(ks),
        "supported_cell_fraction": float(n[edges].sum() / n.sum()),
        "effective_donor_support": float(effective_donor_support),
        "donor_balance": float(donor_balance),
        "replication_adequacy": float(replication_adequacy),
        "replication_formula": "weighted_harmonic_v1",
        "replication_adequacy_role": "audit_metadata_only_not_used_in_fusion",
        "design_rank": int(design_rank),
        "expected_design_rank": int(expected_rank),
        "identifiability": float(identifiability),
        "connected_components": int(components),
        "identifiable": identifiable,
        "support_tau": float(tau),
        "component_selection_used": False,
        "dataset_dependent_component_weight_used": False,
    }

    coverage = (
        cross.rename_axis(index="donor", columns="celltype")
        .stack()
        .rename("n_cells")
        .reset_index()
    )
    coverage["support"] = support.ravel()
    coverage["supported_edge"] = edges.ravel()

    by_type = pd.DataFrame(
        {
            "celltype": labels,
            "n_cells": n.sum(axis=0).astype(int),
            "effective_support": effective_by_type,
            "weight": type_weight,
        }
    )
    by_donor = pd.DataFrame(
        {
            "donor": donors,
            "n_cells": n.sum(axis=1).astype(int),
            "mean_support": donor_mean_support,
        }
    )

    return DesignAuditResult(summary, coverage, by_type, by_donor)
