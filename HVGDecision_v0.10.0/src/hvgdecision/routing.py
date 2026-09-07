"""Metadata-only replication and supported-design identifiability audit."""
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components


@dataclass(frozen=True)
class RoutingConfig:
    tau: float = 15.0
    lower_boundary: float = 2.0
    upper_boundary: float = 4.0
    donor_aware_threshold: float = 0.75

    def __post_init__(self):
        values = np.asarray(list(asdict(self).values()), dtype=float)
        if not np.isfinite(values).all() or self.tau <= 0:
            raise ValueError('Routing parameters must be finite and tau must be positive')
        if self.upper_boundary <= self.lower_boundary:
            raise ValueError('upper_boundary must exceed lower_boundary')
        if not 0 <= self.donor_aware_threshold <= 1:
            raise ValueError('donor_aware_threshold must lie in [0, 1]')


@dataclass
class RoutingResult:
    summary: dict
    coverage: pd.DataFrame
    by_celltype: pd.DataFrame
    by_donor: pd.DataFrame

    @property
    def mode(self):
        return self.summary['selected_mode']

    @property
    def audit(self):
        return pd.DataFrame([self.summary])


def selection_mask(adata, batch_key, label_key, reference=None):
    if batch_key == label_key:
        raise ValueError('batch_key and label_key must be distinct columns')
    for key in (batch_key, label_key):
        if key not in adata.obs:
            raise KeyError(f'Missing adata.obs column: {key!r}')
        if adata.obs[key].isna().any():
            raise ValueError(f'{key!r} contains missing metadata; resolve it before refinement')
        if adata.obs[key].astype(str).str.strip().eq('').any():
            raise ValueError(f'{key!r} contains empty metadata')
    if adata.n_obs == 0:
        raise ValueError('AnnData contains no cells')
    if reference is None:
        return np.ones(adata.n_obs, dtype=bool)
    if isinstance(reference, str):
        raise TypeError('reference must be a list of donor IDs, not a string')
    reference = list(map(str, reference))
    if not reference or len(reference) != len(set(reference)):
        raise ValueError('reference must contain distinct donor IDs')
    missing = set(reference) - set(adata.obs[batch_key].astype(str))
    if missing:
        raise ValueError(f'Reference donors not found: {sorted(missing)}')
    return adata.obs[batch_key].astype(str).isin(reference).to_numpy()


def audit_design(adata, *, batch_key, label_key, reference=None, config=None):
    """Audit selection-cohort metadata without accessing expression or outcomes.

    Unsupported strata do not establish identifiability. A disconnected design,
    fewer than two supported donors or fewer than two supported cell types returns
    ``insufficient_confounded``; it is not routed to a fallback removal engine.
    """
    config = config or RoutingConfig()
    if not isinstance(config, RoutingConfig):
        raise TypeError('config must be RoutingConfig')
    mask = selection_mask(adata, batch_key, label_key, reference)
    obs = adata.obs.loc[mask, [batch_key, label_key]].astype(str)
    donors = sorted(obs[batch_key].unique())
    labels = sorted(obs[label_key].unique())
    cross = pd.crosstab(obs[batch_key], obs[label_key]).reindex(
        index=donors, columns=labels, fill_value=0)
    n = cross.to_numpy(dtype=float)
    d, k = n.shape
    support = n / (n + config.tau)
    r = support.sum(axis=0)
    v = np.minimum(1.0, n.sum(axis=0) / (d * config.tau))
    effective = float(v.sum() / np.sum(v / (r + 1e-12)))
    u = support.mean(axis=1)
    balance = float(np.clip(u.sum()**2 / (d * np.square(u).sum()), 0, 1))
    adequacy = effective * np.sqrt(balance)
    weight = float(np.clip((adequacy - config.lower_boundary) /
                          (config.upper_boundary - config.lower_boundary), 0, 1))

    edges = n >= config.tau
    active_d, active_k = edges.any(axis=1), edges.any(axis=0)
    supported = edges[np.ix_(active_d, active_k)]
    ds, ks = supported.shape
    expected = max(ds + ks - 1, 0)
    if supported.any():
        graph = sparse.bmat([[None, sparse.csr_matrix(supported)],
                             [sparse.csr_matrix(supported.T), None]], format='csr')
        components = int(connected_components(graph, directed=False, return_labels=False))
        # Exact rank of the intercept + donor/type dummy design of a bipartite graph.
        # This avoids allocating a dense design matrix for hundreds of donors/types.
        rank = ds + ks - components
    else:
        components, rank = 0, 0
    iota = float(rank / expected) if expected else 0.0
    identifiable = ds >= 2 and ks >= 2 and components == 1 and rank == expected
    if not identifiable:
        mode, reason = 'insufficient_confounded', 'insufficient_or_disconnected_supported_design'
    elif weight >= config.donor_aware_threshold:
        mode, reason = 'donor_aware', 'identifiable_and_w_at_or_above_threshold'
    else:
        mode, reason = 'cell_level', 'identifiable_but_w_below_threshold'
    summary = dict(
        selected_mode=mode, reason=reason, n_cells=len(obs), n_donors=d, n_celltypes=k,
        n_supported_donors=ds, n_supported_celltypes=ks,
        supported_cell_fraction=float(n[edges].sum() / n.sum()),
        effective_donor_support=effective, donor_balance=balance,
        replication_adequacy=float(adequacy), donor_aware_weight=weight,
        design_rank=rank, expected_design_rank=expected, identifiability=iota,
        connected_components=components, identifiable=bool(identifiable),
        **asdict(config),
    )
    coverage = cross.rename_axis(index='donor', columns='celltype').stack().rename('n_cells').reset_index()
    coverage['support'] = support.ravel()
    coverage['supported_edge'] = edges.ravel()
    by_type = pd.DataFrame(dict(celltype=labels, n_cells=n.sum(axis=0).astype(int),
                               effective_support=r, weight=v))
    by_donor = pd.DataFrame(dict(donor=donors, n_cells=n.sum(axis=1).astype(int), mean_support=u))
    return RoutingResult(summary, coverage, by_type, by_donor)
