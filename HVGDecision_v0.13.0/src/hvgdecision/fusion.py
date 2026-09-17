"""Validated public configuration around the frozen continuous-v0.3 mathematics."""
from dataclasses import dataclass, asdict
import math

from ._continuous_core import (
    ContinuousFusionConfig as _CoreConfig,
    continuous_softgate_fusion as _fuse,
    CellGateConfig, DonorGateConfig, sigmoid_weight,
)


@dataclass(frozen=True)
class ContinuousFusionConfig(_CoreConfig):
    def __post_init__(self):
        for name in ('weight_center', 'weight_scale', 'soft_temperature',
                     'unified_cutoff', 'boundary_epsilon', 'scale_floor'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be finite')
        for name in ('weight_scale', 'soft_temperature', 'boundary_epsilon', 'scale_floor'):
            if getattr(self, name) <= 0:
                raise ValueError(f'{name} must be > 0')
        if not 0 < self.unified_cutoff < 1:
            raise ValueError('unified_cutoff must lie strictly between 0 and 1')
        if self.weight_override is not None:
            if (isinstance(self.weight_override, bool) or
                not math.isfinite(self.weight_override) or not 0 <= self.weight_override <= 1):
                raise ValueError('weight_override must be finite and in [0,1]')
        for obj, cls in ((self.cell, CellGateConfig), (self.donor, DonorGateConfig)):
            if not isinstance(obj, cls):
                raise TypeError(f'Expected {cls.__name__}')
            for key, value in asdict(obj).items():
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError(f'{key} must be finite')
            for key in ('alpha', 'bootstrap_pass_fraction', 'technical_quantile', 'risk_quantile',
                        'maximum_biological_support_quantile', 'minimum_direction_agreement',
                        'minimum_single_donor_dominance', 'lodo_positive_fraction_ceiling'):
                if hasattr(obj, key) and not 0 <= getattr(obj, key) <= 1:
                    raise ValueError(f'{key} must lie in [0,1]')
        if not isinstance(self.donor.minimum_replicates, int) or self.donor.minimum_replicates < 1:
            raise ValueError('minimum_replicates must be a positive integer')


def continuous_softgate_fusion(donor_evidence, cell_evidence, *, replication_adequacy, config=None):
    config = config or ContinuousFusionConfig()
    if not isinstance(config, ContinuousFusionConfig):
        raise TypeError('config must be ContinuousFusionConfig')
    if not math.isfinite(replication_adequacy) or replication_adequacy < 0:
        raise ValueError('replication_adequacy must be finite and nonnegative')
    return _fuse(donor_evidence, cell_evidence,
                 replication_adequacy=replication_adequacy, config=config)
