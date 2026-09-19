"""HVGDecision public API."""

from .api import CountSourceResult, find_raw_counts
from .design import DesignAuditResult, DesignConfig, audit_design
from .fusion import (
    CellGateConfig,
    DonorGateConfig,
    DualEvidenceFusionConfig,
    UnifiedFusionConfig,
    centered_logmeanexp,
    unified_margin_fusion,
)
from .parameters import CellLevelConfig, DonorAwareConfig

__version__ = "0.14.0"
REFINEMENT_MODES = ("dual_evidence",)


def refine(*args, **kwargs):
    """Run the current HVGDecision refinement workflow.

    Imported lazily so the lightweight mathematical API can be inspected before
    optional scientific I/O dependencies are loaded.
    """

    from .refinement_workflow import refine as _refine

    return _refine(*args, **kwargs)


def __getattr__(name):
    if name == "RefinementResult":
        from .refinement_workflow import RefinementResult

        return RefinementResult
    raise AttributeError(name)


__all__ = [
    "refine",
    "audit_design",
    "RefinementResult",
    "DesignAuditResult",
    "DesignConfig",
    "UnifiedFusionConfig",
    "DualEvidenceFusionConfig",
    "centered_logmeanexp",
    "unified_margin_fusion",
    "CellGateConfig",
    "DonorGateConfig",
    "CellLevelConfig",
    "DonorAwareConfig",
    "CountSourceResult",
    "find_raw_counts",
    "REFINEMENT_MODES",
]
