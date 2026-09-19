"""Public design-audit API for HVGDecision."""

from dataclasses import dataclass
import math

from ._design_core import DesignAuditResult, audit_metadata


@dataclass(frozen=True)
class DesignConfig:
    """Configuration for donor-by-cell-type support auditing."""

    tau: float = 15.0

    def __post_init__(self):
        if isinstance(self.tau, bool) or not isinstance(self.tau, (int, float)):
            raise ValueError("tau must be a finite positive number")
        if not math.isfinite(self.tau) or self.tau <= 0:
            raise ValueError("tau must be a finite positive number")


def audit_design(adata, *, batch_key, label_key, reference=None, config=None):
    """Audit discovery-cohort replication adequacy and identifiability.

    The returned quantities are metadata audits only and are not used to weight
    or select either gene-level evidence component.
    """

    config = config or DesignConfig()
    if not isinstance(config, DesignConfig):
        raise TypeError("config must be DesignConfig")
    return audit_metadata(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        reference=reference,
        tau=config.tau,
    )


__all__ = ["DesignConfig", "DesignAuditResult", "audit_design"]
