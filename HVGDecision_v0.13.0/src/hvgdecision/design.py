"""Weighted replication adequacy; design identifiability is a separate audit."""
from dataclasses import dataclass
import math

from .routing import audit_design as _legacy_audit, RoutingConfig
from .fusion import ContinuousFusionConfig, sigmoid_weight


@dataclass(frozen=True)
class DesignConfig:
    tau: float = 15.0

    def __post_init__(self):
        if isinstance(self.tau, bool) or not math.isfinite(self.tau) or self.tau <= 0:
            raise ValueError('tau must be finite and positive')


def audit_design(adata, *, batch_key, label_key, reference=None, config=None, fusion_config=None):
    """Compute weighted A from discovery metadata only; never selects a risk branch.

    The weight column in by_celltype is v_k = min(1, n_k/(D*tau)).
    effective_donor_support = sum(v_k) / sum(v_k/(r_k+1e-12)).
    A = effective_donor_support * sqrt(donor_balance).
    ``identifiable`` additionally requires >=2 supported donors/types and a
    connected supported design. A high A alone is not an identifiability claim.
    """
    config = config or DesignConfig()
    fusion_config = fusion_config or ContinuousFusionConfig()
    if not isinstance(config, DesignConfig):
        raise TypeError('config must be DesignConfig')
    if not isinstance(fusion_config, ContinuousFusionConfig):
        raise TypeError('fusion_config must be ContinuousFusionConfig')
    result = _legacy_audit(adata, batch_key=batch_key, label_key=label_key,
                           reference=reference, config=RoutingConfig(tau=config.tau))
    summary = result.summary
    for key in ('lower_boundary', 'upper_boundary', 'donor_aware_threshold'):
        summary.pop(key, None)
    automatic = sigmoid_weight(summary['replication_adequacy'],
                               fusion_config.weight_center, fusion_config.weight_scale)
    used = automatic if fusion_config.weight_override is None else fusion_config.weight_override
    summary.update(selected_mode='continuous', route_selection_used=False,
                   reason='design_audit_only_no_hard_branch_selection',
                   replication_formula='weighted_harmonic_v1',
                   automatic_donor_aware_weight=automatic, donor_aware_weight=used,
                   cell_level_weight=1-used, weight_center=fusion_config.weight_center,
                   weight_scale=fusion_config.weight_scale,
                   weight_source=('sigmoid_replication_adequacy' if fusion_config.weight_override is None
                                  else 'explicit_endpoint_or_sensitivity_override'))
    return result
