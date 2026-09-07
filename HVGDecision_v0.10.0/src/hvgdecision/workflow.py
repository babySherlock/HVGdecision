"""Current multi-donor panel refinement API (no integration fitting or budget search)."""
from dataclasses import asdict, dataclass
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .api import CountSourceResult, find_raw_counts
from .count_validation import audit_counts
from .io import count_matrix, count_var_names, gene_hash, is_raw_source
from .parameters import CellLevelConfig, DonorAwareConfig
from .routing import RoutingConfig, RoutingResult, audit_design, selection_mask


@dataclass
class RefinementResult:
    adata: ad.AnnData
    decision_table: pd.DataFrame
    routing: RoutingResult
    counts_audit: pd.DataFrame
    run_info: dict

    @property
    def selected_mode(self):
        return self.run_info['selected_mode']

    @property
    def removed_genes(self):
        return self.decision_table.loc[self.decision_table.final_action.eq('remove'), 'gene'].tolist()

    @property
    def harmful_genes(self):
        return self.removed_genes

    @property
    def final_n_hvg(self):
        return self.adata.n_vars

    @property
    def base_n_hvg(self):
        return len(self.decision_table)

    def __repr__(self):
        return (f'RefinementResult(mode={self.selected_mode!r}, base={self.base_n_hvg}, '
                f'removed={len(self.removed_genes)}, final={self.final_n_hvg})')

    def save(self, output_dir):
        """Save a self-contained audit and raw-count AnnData to a NEW/empty folder."""
        path = Path(output_dir).expanduser()
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise FileExistsError(f'Output directory must be empty; existing results preserved: {path}')
        path.mkdir(parents=True, exist_ok=True)
        self.decision_table.to_csv(path / 'gene_decisions.csv', index=False)
        self.decision_table.loc[self.decision_table.final_action.eq('remove')].to_csv(
            path / 'removed_genes.csv', index=False)
        self.decision_table.loc[self.decision_table.final_action.eq('keep')].to_csv(
            path / 'retained_genes.csv', index=False)
        self.counts_audit.to_csv(path / 'counts_audit.csv', index=False)
        self.routing.audit.to_csv(path / 'routing_summary.csv', index=False)
        self.routing.coverage.to_csv(path / 'donor_celltype_coverage.csv', index=False)
        self.routing.by_donor.to_csv(path / 'routing_by_donor.csv', index=False)
        self.routing.by_celltype.to_csv(path / 'routing_by_celltype.csv', index=False)
        (path / 'run_manifest.json').write_text(json.dumps(self.run_info, ensure_ascii=False, indent=2))
        self.adata.write_h5ad(path / 'adata_hvg.h5ad', compression='gzip')
        return path


def _gene_list(value, name):
    if isinstance(value, (str, Path)):
        raise TypeError(f'{name} must be a gene-ID sequence, not a filename/string; read the table first')
    result = list(value)
    if any(pd.isna(g) or not str(g).strip() for g in result):
        raise ValueError(f'{name} contains missing/empty gene IDs')
    result = list(map(str, result))
    if not result or len(set(result)) != len(result):
        raise ValueError(f'{name} must contain nonempty, distinct gene IDs')
    return result


def _select_hvgs(counts, genes, obs, batch_key, n_hvg):
    import scanpy as sc
    if isinstance(n_hvg, bool) or not isinstance(n_hvg, int) or n_hvg < 2:
        raise ValueError('n_hvg must be an integer >= 2')
    if n_hvg > len(genes):
        raise ValueError(f'Requested {n_hvg} HVGs but counts have only {len(genes)} genes')
    work = ad.AnnData(X=counts.copy(), obs=obs.copy(), var=pd.DataFrame(index=genes))
    sc.pp.highly_variable_genes(work, flavor='seurat_v3', batch_key=batch_key,
                               n_top_genes=n_hvg, subset=False, inplace=True)
    selected = work.var.loc[work.var.highly_variable].copy()
    selected['gene_id_tiebreak'] = selected.index.astype(str)
    selected = selected.sort_values(
        ['highly_variable_rank', 'highly_variable_nbatches', 'gene_id_tiebreak'],
        ascending=[True, False, True], kind='mergesort')
    if len(selected) != n_hvg:
        raise RuntimeError(f'Seurat v3 returned {len(selected)} genes, expected {n_hvg}; inspect input QC')
    return selected.index.astype(str).tolist(), work.var


def refine(adata, *, batch_key, label_key, counts=None, n_hvg=2000, hvg_genes=None,
           reference=None, protected_genes=None, mode='auto', routing_config=None,
           cell_level_config=None, donor_aware_config=None, seed=20260829,
           output_dir=None, return_details=False):
    """Return a raw-count AnnData after auditable removal from a fixed base panel.

    By default all cells are the selection cohort; ``reference=[donor IDs]``
    restricts HVG selection, routing and risk evidence to that cohort. Output
    always retains ALL input cells. No held-out evaluation is performed here.
    An external ordered ``hvg_genes`` list bypasses Seurat v3 and determines
    the actual base size. Gene IDs must exactly match the chosen counts axis.
    ``mode='auto'`` is the current algorithm. Explicit branch overrides are
    sensitivity analyses and cannot bypass the identifiability safeguard.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError('adata must be an AnnData object')
    if adata.isbacked:
        raise ValueError('Use adata.to_memory() before refinement')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('AnnData cell and gene IDs must be unique; resolve duplicates explicitly')
    if mode not in ('auto', 'cell_level', 'donor_aware'):
        raise ValueError('mode must be auto, cell_level or donor_aware')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a non-negative integer')
    if output_dir is not None:
        path = Path(output_dir).expanduser()
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise FileExistsError(f'Choose a new/empty output directory: {path}')
    routing_config = routing_config or RoutingConfig()
    cell_level_config = cell_level_config or CellLevelConfig()
    donor_aware_config = donor_aware_config or DonorAwareConfig()
    for obj, cls in [(routing_config, RoutingConfig), (cell_level_config, CellLevelConfig),
                     (donor_aware_config, DonorAwareConfig)]:
        if not isinstance(obj, cls):
            raise TypeError(f'Expected {cls.__name__}, received {type(obj).__name__}')
    mask = selection_mask(adata, batch_key, label_key, reference)
    routing = audit_design(adata, batch_key=batch_key, label_key=label_key,
                           reference=reference, config=routing_config)
    source = counts if isinstance(counts, CountSourceResult) else find_raw_counts(
        adata, source='auto' if counts is None else counts)
    if not source.valid:
        raise ValueError(f'No valid raw counts: {source.error}\n{source.audit.to_string(index=False)}')
    if source.obs_names is not None and not pd.Index(source.obs_names).equals(adata.obs_names.astype(str)):
        raise ValueError('CountSourceResult belongs to a different cell axis/order; rerun find_raw_counts for this AnnData')
    matrix = source.matrix if source.source == 'external' else count_matrix(adata, source.source)
    genes = pd.Index(source.gene_names if source.source == 'external' else count_var_names(adata, source.source)).astype(str)
    if matrix.shape != (adata.n_obs, len(genes)) or not genes.is_unique:
        raise ValueError('Count-source dimensions or gene IDs do not match the current input')
    # Revalidate to catch edits made after find_raw_counts() returned.
    numeric = audit_counts(matrix, source.location)
    if not numeric['valid']:
        raise ValueError(numeric['error'])
    matrix = sparse.csr_matrix(matrix) if sparse.issparse(matrix) else np.asarray(matrix)
    if np.any(np.asarray(matrix.sum(axis=1)).ravel() <= 0):
        raise ValueError('Counts contain zero-library cells; perform QC before refinement')
    cohort_counts = matrix[mask, :]
    cohort = ad.AnnData(obs=adata.obs.loc[mask].copy())
    if hvg_genes is None:
        panel, feature_metadata = _select_hvgs(cohort_counts, genes, cohort.obs, batch_key, n_hvg)
        feature_method = 'scanpy_seurat_v3_batch_aware'
    else:
        panel = _gene_list(hvg_genes, 'hvg_genes')
        if len(panel) < 2:
            raise ValueError('At least two panel genes are required')
        feature_metadata = pd.DataFrame(index=genes)
        feature_method = 'external_ordered_panel'
    missing = set(panel) - set(genes)
    if missing:
        raise ValueError(f'{len(missing)} panel genes absent from counts: {sorted(missing)[:20]}; no silent dropping/renaming')
    protected = [] if protected_genes is None or len(protected_genes) == 0 else _gene_list(protected_genes, 'protected_genes')
    unknown_protected = sorted(set(protected) - set(genes))
    if unknown_protected:
        raise ValueError(f'Protected gene IDs absent from counts: {unknown_protected[:20]}')
    actual = routing.mode
    if actual == 'insufficient_confounded':
        warnings.warn('Supported donor-by-cell-type design is insufficient/confounded; original panel retained.', UserWarning, stacklevel=2)
        table = pd.DataFrame(dict(gene=panel, final_action='keep', harmful=False,
                                  explicitly_protected=[g in protected for g in panel],
                                  decision_reason='not_tested_insufficient_confounded_design'))
        engine_id = 'not_run'
        diagnostics = {}
    else:
        if mode != 'auto':
            actual = mode
            warnings.warn('Explicit risk-branch override is a sensitivity analysis, not automatic routing.', UserWarning, stacklevel=2)
        if actual == 'cell_level':
            from ._cell_level import find_harmful_v1_celllevel
            result = find_harmful_v1_celllevel(cohort, cohort_counts, genes, panel, batch_key,
                                             label_key, protected_genes=protected, seed=seed,
                                             **asdict(cell_level_config))
            table = result['table'].copy()
            table['final_action'] = np.where(table.harmful, 'remove', 'keep')
            gates = {
                'permutation_fdr': table.permutation_fdr <= cell_level_config.alpha,
                'donor_leakage_z': table.donor_leakage_z >= cell_level_config.leakage_z_floor,
                'risk_score_z': table.risk_score_z >= cell_level_config.risk_z_floor,
                'biology_z': table.biology_z <= cell_level_config.biology_z_ceiling,
                'bootstrap': table.bootstrap_risk_fraction >= cell_level_config.bootstrap_pass_fraction,
                'marker_protection': ~table.marker_protected,
                'explicit_protection': ~table.explicitly_protected,
            }
            engine_id = 'generic_cell_level_V1'
            diagnostics = {'sampled_cells': result['sampled_cells']}
        else:
            from ._donor_aware import find_harmful_within_protocol_donor_replicate
            result = find_harmful_within_protocol_donor_replicate(
                cohort, cohort_counts, genes, panel, batch_key, label_key,
                protected_genes=protected, seed=seed, risk_config=asdict(donor_aware_config))
            table = result['table'].copy()
            table['harmful'] = table.final_action.eq('remove')
            table['explicitly_protected'] = table.user_protected
            gates = {
                'technical_score': table.donor_nonreproducibility_score >= table.technical_threshold,
                'risk_score': table.risk_score >= table.risk_threshold,
                'biological_support': table.biological_support_score <= table.biological_support_ceiling,
                'replicate_count': table.replicate_count >= 3,
                'dominance_or_lodo': (table.single_group_effect_dominance >= donor_aware_config.minimum_single_donor_dominance) | (table.lodo_positive_fraction < 0.8),
                'direction_agreement': table.celltype_effect_direction_agreement < donor_aware_config.minimum_direction_agreement,
                'replicated_marker_protection': ~table.replicated_marker_protection,
                'rare_marker_protection': ~table.rare_marker_protection,
                'explicit_protection': ~table.user_protected,
            }
            engine_id = 'scenario_specific_v2.0_donor_replicate'
            diagnostics = {}
        for gate, passed in gates.items():
            table['pass_' + gate] = passed
        gate_frame = pd.DataFrame(gates)
        table['failed_rules'] = gate_frame.apply(lambda row: '|'.join(row.index[~row.astype(bool)]), axis=1)
        if not np.array_equal(gate_frame.all(axis=1).to_numpy(), table.harmful.to_numpy()):
            raise RuntimeError('Audit gates do not match the frozen engine removal calls')
    table['input_rank'] = table.gene.map({g: i+1 for i, g in enumerate(panel)})
    table = table.sort_values('input_rank', kind='mergesort').reset_index(drop=True)
    table['selected_mode'] = actual
    table['in_final_panel'] = table.final_action.eq('keep')
    annotation = adata.raw.var if is_raw_source(source.source) and adata.raw is not None else adata.var
    for column in ['gene_symbol', 'gene_symbols', 'feature_name', 'ensembl_id', 'gene_ids', 'feature_id']:
        if column in annotation and column not in table:
            values = annotation[column].reindex(table.gene).astype(object)
            table[column] = values.where(values.notna(), '').astype(str).to_numpy()
    retained = table.loc[table.in_final_panel, 'gene'].tolist()
    removed = table.loc[~table.in_final_panel, 'gene'].tolist()
    if not retained:
        raise RuntimeError('All panel genes would be removed; no empty AnnData exported')
    if set(removed) & set(protected):
        raise RuntimeError('Internal error: protected gene removed')
    selected_counts = matrix[:, genes.get_indexer(retained)].copy()
    if is_raw_source(source.source) and adata.raw is not None:
        var = adata.raw.var.reindex(retained).copy()
    else:
        var = adata.var.reindex(retained).copy()
    for col in feature_metadata:
        var[col] = feature_metadata.reindex(retained)[col]
    var['highly_variable'] = True
    var['hvgdecision_input_rank'] = table.set_index('gene').loc[retained, 'input_rank'].to_numpy()
    final = ad.AnnData(X=selected_counts, obs=adata.obs.copy(), var=var)
    final.layers['counts'] = selected_counts.copy()
    # Old embeddings/graphs/normalised layers are deliberately not copied.
    versions = {}
    for package in ['hvgdecision', 'scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'scikit-misc']:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = 'not_installed'
    engine_hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                     for name in ['_cell_level.py', '_donor_aware.py', 'routing.py']}
    cohort_ids = adata.obs_names[mask].astype(str).tolist()
    manifest = dict(
        schema_version='0.10.0', software_versions=versions, engine_sha256=engine_hashes,
        engine_id=engine_id, requested_mode=mode, selected_mode=actual, automatic_routing=routing.summary,
        feature_method=feature_method, base_n_hvg=len(panel), final_n_hvg=len(retained),
        removed_genes=removed, base_panel_sha256=gene_hash(panel), final_panel_sha256=gene_hash(retained),
        batch_key=batch_key, label_key=label_key, selection_donors=sorted(cohort.obs[batch_key].astype(str).unique()),
        selection_cell_ids_sha256=gene_hash(cohort_ids), n_selection_cells=len(cohort_ids),
        selection_scope='all_input_cells' if reference is None else 'specified_reference_donors',
        selection_uses_celltype_labels=True, downstream_evaluation_performed=False,
        raw_counts_location=source.location, counts_numeric_validation=numeric,
        raw_provenance_note='Integer-like counts do not establish experimental raw-count provenance.',
        protected_genes=protected, protected_genes_outside_panel=sorted(set(protected)-set(panel)),
        routing_config=asdict(routing_config), cell_level_config=asdict(cell_level_config),
        donor_aware_config=asdict(donor_aware_config), seed=seed, diagnostics=diagnostics,
    )
    final.uns['hvgdecision'] = dict(
        manifest_json=json.dumps(manifest, ensure_ascii=False),
        gene_decisions=table.copy(), routing_summary=routing.audit,
        removed_genes=np.asarray(removed, dtype=str), base_genes=np.asarray(panel, dtype=str),
        selection_cell_ids=np.asarray(cohort_ids, dtype=str), selected_mode=actual,
    )
    updated_counts_audit = source.audit.copy()
    for col, val in numeric.items():
        if col == 'source':
            continue
        selected_rows = updated_counts_audit['selected'].fillna(False).astype(bool)
        updated_counts_audit.loc[selected_rows, col] = val
    output = RefinementResult(final, table, routing, updated_counts_audit, manifest)
    if output_dir is not None:
        output.save(output_dir)
    return output if return_details else final
