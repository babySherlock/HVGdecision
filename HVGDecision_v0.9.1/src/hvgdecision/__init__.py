"""HVGDecision public API."""

from .api import (
    BudgetSearchResult,
    CountSourceResult,
    HVGStudy,
    HVGRefinementResult,
    find_raw_counts,
    prepare_cross_domain_inputs,
    setup_reference_query,
)
from .modes import normalize_mode
from .scoring import three_domain_scores

__version__ = "0.9.1"
VALID_MODES = ("within_domain", "cross_domain")

__all__ = [
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
