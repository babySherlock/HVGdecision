"""HVGDecision public API."""

from .api import (
    BudgetSearchResult,
    CountSourceResult,
    HVGStudy,
    HVGRefinementResult,
    find_raw_counts,
    prepare_cross_domain_inputs,
    setup_reference_query as _legacy_setup_reference_query,
)
from .modes import normalize_mode
from .scoring import three_domain_scores

from .routing import RoutingConfig, RoutingResult, audit_design
from .parameters import CellLevelConfig, DonorAwareConfig
from .workflow import RefinementResult, refine

__version__ = "0.10.0"
# Retained for the legacy normalize_mode/setup_reference_query API only.
VALID_MODES = ("within_domain", "cross_domain")
REFINEMENT_MODES = ("auto", "cell_level", "donor_aware")


def setup_reference_query(*args, **kwargs):
    """Legacy 0.9 API; use refine/audit_design for the current manuscript workflow."""
    import warnings
    warnings.warn('setup_reference_query is the legacy 0.9 workflow and does not use '
                  'automatic risk-branch routing. Use hd.refine for the current workflow.',
                  FutureWarning, stacklevel=2)
    return _legacy_setup_reference_query(*args, **kwargs)

__all__ = [
    "refine", "audit_design", "RefinementResult", "RoutingResult",
    "RoutingConfig", "CellLevelConfig", "DonorAwareConfig", "REFINEMENT_MODES",
    "BudgetSearchResult",
    "CountSourceResult",
    "HVGStudy",
    "HVGRefinementResult",
    "VALID_MODES",
    "__version__",
    "find_raw_counts",
    "prepare_cross_domain_inputs",
    "normalize_mode",
    "setup_reference_query",
    "three_domain_scores",
]
