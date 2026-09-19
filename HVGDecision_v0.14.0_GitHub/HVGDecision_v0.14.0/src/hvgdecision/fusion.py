"""Dual-evidence signed-margin fusion for HVGDecision v0.14.0."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Optional

import numpy as np
import pandas as pd

from ._signed_margin_core import (
    CellGateConfig,
    DonorGateConfig,
    cell_soft_gate,
    donor_soft_gate,
)


@dataclass(frozen=True)
class UnifiedFusionConfig:
    """Configuration for centered log-mean-exp fusion.

    Parameters
    ----------
    tau
        Fusion temperature. Smaller values approach the larger effective
        component margin; larger values impose a stronger disagreement penalty.
    discovery_cutoff
        Final threshold on the unified signed margin. The default 0 preserves
        the signed interpretation: positive values are discovery-risk calls.
    soft_temperature
        Internal scale used only while reconstructing legacy gate-consistency
        diagnostics. It does not enter the final v0.14 fusion equation.
    boundary_epsilon
        Numerical sign-preserving epsilon for exact gate boundaries.
    scale_floor
        Minimum acceptable robust scale when standardizing gate distances.
    """

    tau: float = 0.10
    discovery_cutoff: float = 0.0
    soft_temperature: float = 1.0
    boundary_epsilon: float = 1e-9
    scale_floor: float = 1e-8
    donor: DonorGateConfig = DonorGateConfig()
    cell: CellGateConfig = CellGateConfig()

    def __post_init__(self):
        for name in (
            "tau",
            "discovery_cutoff",
            "soft_temperature",
            "boundary_epsilon",
            "scale_floor",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        for name in ("tau", "soft_temperature", "boundary_epsilon", "scale_floor"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be > 0")
        if not isinstance(self.donor, DonorGateConfig):
            raise TypeError("donor must be DonorGateConfig")
        if not isinstance(self.cell, CellGateConfig):
            raise TypeError("cell must be CellGateConfig")


# Friendly synonym used in documentation.
DualEvidenceFusionConfig = UnifiedFusionConfig


def centered_logmeanexp(donor_margin, cell_margin, *, tau: float = 0.10) -> np.ndarray:
    """Return the centered log-mean-exp of two signed-margin arrays.

    M_u = tau * log((exp(M_d/tau) + exp(M_c/tau)) / 2)
    """

    tau = float(tau)
    if not math.isfinite(tau) or tau <= 0:
        raise ValueError("tau must be finite and > 0")

    md = np.asarray(donor_margin, dtype=float)
    mc = np.asarray(cell_margin, dtype=float)
    if md.shape != mc.shape:
        raise ValueError("donor_margin and cell_margin must have the same shape")

    return tau * (np.logaddexp(md / tau, mc / tau) - np.log(2.0))


def _effective_margin(table: pd.DataFrame, margin_col: str, protection_col: str) -> np.ndarray:
    margin = pd.to_numeric(table[margin_col], errors="raise").to_numpy(dtype=float)
    protected = (
        table[protection_col].fillna(False).astype(bool).to_numpy()
        if protection_col in table.columns
        else np.zeros(len(table), dtype=bool)
    )
    return np.where(protected, -np.inf, margin)


def _drop_obsolete_internal_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Remove v0.3 confidence/debug fields from the public v0.14 table."""

    exact = {
        "donor_confidence_pre_protection",
        "donor_confidence",
        "donor_endpoint_call_from_confidence",
        "donor_packaged_risk_gate_match",
        "donor_packaged_final_call_match",
        "cell_confidence_pre_protection",
        "cell_confidence",
        "cell_endpoint_call_from_confidence",
        "cell_packaged_risk_gate_match",
        "cell_packaged_final_call_match",
    }
    return frame.drop(columns=[c for c in exact if c in frame.columns], errors="ignore")


def unified_margin_fusion(
    donor_evidence: pd.DataFrame,
    cell_evidence: pd.DataFrame,
    *,
    replication_adequacy: Optional[float] = None,
    config: Optional[UnifiedFusionConfig] = None,
) -> tuple[pd.DataFrame, dict]:
    """Fuse donor-aware and cell-level evidence on a shared signed-margin scale.

    Replication adequacy is accepted only so it can be recorded as audit metadata;
    it is not used in the mathematical fusion.
    """

    cfg = config or UnifiedFusionConfig()
    if not isinstance(cfg, UnifiedFusionConfig):
        raise TypeError("config must be UnifiedFusionConfig")

    donor = donor_soft_gate(
        donor_evidence,
        config=cfg.donor,
        soft_temperature=cfg.soft_temperature,
        boundary_epsilon=cfg.boundary_epsilon,
        scale_floor=cfg.scale_floor,
    )
    cell = cell_soft_gate(
        cell_evidence,
        config=cfg.cell,
        soft_temperature=cfg.soft_temperature,
        boundary_epsilon=cfg.boundary_epsilon,
        scale_floor=cfg.scale_floor,
    )

    donor_genes = donor["gene"].astype(str)
    cell_genes = cell["gene"].astype(str)
    if set(donor_genes) != set(cell_genes):
        only_d = sorted(set(donor_genes) - set(cell_genes))
        only_c = sorted(set(cell_genes) - set(donor_genes))
        raise ValueError(
            "Donor and cell evidence must describe the same candidate panel. "
            f"Only donor={only_d[:10]}, only cell={only_c[:10]}"
        )

    cell = cell.set_index("gene").reindex(donor_genes).reset_index()
    table = donor.merge(cell, on="gene", how="left", validate="one_to_one")
    table = _drop_obsolete_internal_columns(table)

    md_raw = pd.to_numeric(table["donor_branch_margin"], errors="raise").to_numpy(dtype=float)
    mc_raw = pd.to_numeric(table["cell_branch_margin"], errors="raise").to_numpy(dtype=float)
    md = _effective_margin(table, "donor_branch_margin", "donor_protected")
    mc = _effective_margin(table, "cell_branch_margin", "cell_protected")
    mu = centered_logmeanexp(md, mc, tau=cfg.tau)
    remove = mu > float(cfg.discovery_cutoff)

    table["donor_branch_margin_raw"] = md_raw
    table["cell_branch_margin_raw"] = mc_raw
    table["donor_branch_margin_effective"] = md
    table["cell_branch_margin_effective"] = mc
    table["unified_margin"] = mu
    table["fusion_tau"] = float(cfg.tau)
    table["discovery_cutoff"] = float(cfg.discovery_cutoff)
    table["discovery_remove"] = remove
    table["final_action"] = np.where(remove, "remove", "keep")
    table["replication_adequacy"] = (
        np.nan if replication_adequacy is None else float(replication_adequacy)
    )

    manifest = {
        "method": "dual_evidence_signed_margin_centered_logmeanexp_v0_4",
        "fusion_formula": "M_u=tau*log((exp(M_d/tau)+exp(M_c/tau))/2)",
        "tau": float(cfg.tau),
        "discovery_cutoff": float(cfg.discovery_cutoff),
        "decision_rule": "discovery_remove iff unified_margin > discovery_cutoff",
        "component_selection_used": False,
        "dataset_dependent_component_weight_used": False,
        "replication_adequacy": None if replication_adequacy is None else float(replication_adequacy),
        "replication_adequacy_role": "audit_metadata_only_not_used_in_fusion",
        "n_candidate_genes": int(len(table)),
        "n_discovery_removed": int(remove.sum()),
        "config": asdict(cfg),
    }
    return table, manifest


__all__ = [
    "UnifiedFusionConfig",
    "DualEvidenceFusionConfig",
    "DonorGateConfig",
    "CellGateConfig",
    "centered_logmeanexp",
    "unified_margin_fusion",
]
