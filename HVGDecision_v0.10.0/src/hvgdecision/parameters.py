"""Explicit, validated parameters for the two frozen risk branches."""
from dataclasses import asdict, dataclass
import math


def _validate(obj, integers=(), probabilities=()):
    for name, value in asdict(obj).items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f'{name} must be a finite number')
        if name in integers and (not isinstance(value, int) or value < 1):
            raise ValueError(f'{name} must be a positive integer')
        if name in probabilities and not 0 <= value <= 1:
            raise ValueError(f'{name} must lie in [0, 1]')


@dataclass(frozen=True)
class CellLevelConfig:
    max_cells_per_batch_label: int = 200
    n_permutations: int = 100
    n_bootstraps: int = 20
    alpha: float = 0.05
    bootstrap_pass_fraction: float = 0.80
    leakage_z_floor: float = 1.0
    risk_z_floor: float = 1.0
    biology_z_ceiling: float = 0.0
    marker_log_effect: float = 0.50
    marker_replication_fraction: float = 0.80
    marker_min_eligible_donors: int = 2
    marker_min_cells_per_side: int = 10

    def __post_init__(self):
        _validate(self, ('max_cells_per_batch_label', 'n_permutations', 'n_bootstraps',
                         'marker_min_eligible_donors', 'marker_min_cells_per_side'),
                        ('alpha', 'bootstrap_pass_fraction', 'marker_replication_fraction'))


@dataclass(frozen=True)
class DonorAwareConfig:
    min_group_cells: int = 15
    technical_quantile: float = 0.975
    risk_quantile: float = 0.985
    minimum_absolute_risk: float = 0.30
    maximum_biological_support_quantile: float = 0.45
    minimum_direction_agreement: float = 0.75
    minimum_single_donor_dominance: float = 0.55

    def __post_init__(self):
        _validate(self, ('min_group_cells',),
                  ('technical_quantile', 'risk_quantile', 'maximum_biological_support_quantile',
                   'minimum_direction_agreement', 'minimum_single_donor_dominance'))
