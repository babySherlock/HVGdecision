"""Current multi-donor panel refinement API (no integration fitting or budget search)."""
from dataclasses import asdict, dataclass, field, replace
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
from .routing import RoutingResult, selection_mask
from .design import DesignConfig, audit_design
from .fusion import ContinuousFusionConfig, CellGateConfig, DonorGateConfig, continuous_softgate_fusion


@dataclass
class RefinementResult:
    adata: ad.AnnData
    decision_table: pd.DataFrame
    routing: RoutingResult
    counts_audit: pd.DataFrame
    run_info: dict
    reference_decision_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    hvg_audits: dict = field(default_factory=dict)
    component_evidence: dict = field(default_factory=dict)

    @property
    def reference_risk_genes(self):
        table = self.reference_decision_table
        return table.loc[table.final_action.eq('remove'), 'gene'].tolist() if not table.empty else []

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
        self.routing.audit.to_csv(path / 'design_audit.csv', index=False)
        self.routing.coverage.to_csv(path / 'donor_celltype_coverage.csv', index=False)
        self.routing.by_donor.to_csv(path / 'routing_by_donor.csv', index=False)
        self.routing.by_celltype.to_csv(path / 'routing_by_celltype.csv', index=False)
        self.reference_decision_table.to_csv(path / 'reference_gene_decisions.csv', index=False)
        self.reference_decision_table.loc[self.reference_decision_table.final_action.eq('remove')].to_csv(
            path / 'reference_risk_genes.csv', index=False)
        for role, audit in self.hvg_audits.items():
            audit.to_csv(path / f'{role}_hvg_full_audit.csv', index=False)
        for component, evidence in self.component_evidence.items():
            evidence.to_csv(path / f'{component}_evidence.csv', index=False)
        pd.DataFrame({'cell_id': self.adata.uns['hvgdecision']['selection_cell_ids']}).to_csv(
            path / 'discovery_cell_ids.csv', index=False)
        for role in ('reference', 'query'):
            panel = self.run_info.get(role + '_base_genes', [])
            if panel:
                pd.DataFrame({'gene': panel, 'rank': range(1, len(panel)+1)}).to_csv(
                    path / f'{role}_hvg_panel.csv', index=False)
        attempts = self.run_info.get('hvg_selection', {})
        pd.DataFrame([dict(role=role, **attempt) for role, info in attempts.items()
                      for attempt in info.get('attempts', [])],
                     columns=['role', 'span', 'success', 'error']).to_csv(
                         path / 'hvg_span_attempts.csv', index=False)
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


def _span_values(value):
    if isinstance(value, str) and value == 'auto':
        return [0.3, 0.5, 0.7, 1.0]
    if isinstance(value, (bool, str)) or not isinstance(value, (int, float)) or not np.isfinite(value) or not 0 < value <= 1:
        raise ValueError("hvg_span must be 'auto' or a finite number in (0, 1]")
    return [float(value)]


def _select_hvgs(counts, genes, obs, batch_key, n_hvg, hvg_span='auto'):
    import scanpy as sc
    if isinstance(n_hvg, bool) or not isinstance(n_hvg, int) or n_hvg < 2:
        raise ValueError('n_hvg must be an integer >= 2')
    if n_hvg > len(genes):
        raise ValueError(f'Requested {n_hvg} HVGs but counts have only {len(genes)} genes')
    spans = _span_values(hvg_span)
    attempts = []
    for span in spans:
        # Fresh fit per attempt; only batch metadata, never Query cell-type labels.
        work = ad.AnnData(X=counts.copy(), obs=pd.DataFrame(
            {batch_key: obs[batch_key].astype(str).to_numpy()}, index=obs.index.copy()),
            var=pd.DataFrame(index=genes))
        try:
            sc.pp.highly_variable_genes(work, flavor='seurat_v3',
                batch_key=batch_key if obs[batch_key].nunique() >= 2 else None,
                n_top_genes=n_hvg, span=span, subset=False, inplace=True)
        except (ValueError, np.linalg.LinAlgError) as error:
            attempts.append(dict(span=span, success=False, error=repr(error)))
            message = str(error).lower()
            numerical = any(term in message for term in (
                'singular', 'reciprocal condition', 'svddc', 'extrapolation not allowed',
                'zero-width neighborhood', 'near singularit'))
            if not numerical or span == spans[-1]:
                raise ValueError(f'Seurat v3 failed; hvg_span={hvg_span!r}; attempts={attempts}. '
                                 'Explicit numeric spans do not silently fall back.') from error
            warnings.warn(f'Seurat v3 span={span:g} numerical failure; retrying next span.',
                          UserWarning, stacklevel=2)
            continue
        attempts.append(dict(span=span, success=True, error=''))
        break
    # Match the original notebook: stable rank order, retaining input-axis order for ties.
    full = work.var.copy()
    full['gene'] = full.index.astype(str)
    full = full.sort_values(['highly_variable', 'highly_variable_rank'],
                            ascending=[False, True], kind='mergesort')
    selected = full.loc[full.highly_variable]
    if len(selected) != n_hvg:
        raise RuntimeError(f'Seurat v3 returned {len(selected)} genes, expected {n_hvg}; inspect input QC')
    full.attrs['hvg_selection'] = dict(requested_span=hvg_span, used_span=span,
                                      method='scanpy_seurat_v3_batch_aware', attempts=attempts)
    return selected.index.astype(str).tolist(), full


def refine(adata, *, batch_key, label_key, counts=None, n_hvg=2000, hvg_genes=None,
           reference=None, protected_genes=None, mode='continuous', design_config=None,
           cell_level_config=None, donor_aware_config=None, seed=20260829,
           output_dir=None, return_details=False, query=None, query_hvg_genes=None,
           hvg_span='auto', query_hvg_span=None, fusion_config=None, design_policy='report'):
    """Return a raw-count AnnData after auditable removal from a fixed base panel.

    By default all cells are the selection cohort; ``reference=[donor IDs]``
    restricts HVG selection, routing and risk evidence to that cohort. Output
    always retains ALL input cells. With ``query=[donor IDs]``, both lists must
    be explicit, disjoint and cover all input cells. Reference and Query HVGs
    are selected independently; only Reference risk calls intersecting the Query
    panel are removed. Query labels are not used. No evaluation is performed.
    ``hvg_span='auto'`` retries numerical LOESS failures at .3/.5/.7/1.0,
    recording all attempts. A number fixes the span with no fallback.
    ``query_hvg_span`` defaults to the same policy as ``hvg_span``.
    An external ordered ``hvg_genes`` list bypasses Seurat v3 and determines
    the actual base size. Gene IDs must exactly match the chosen counts axis.
    A is computed automatically using the weighted harmonic mean on discovery
    cells. Both evidence components are evaluated; no hard branch is selected.
    ``mode='auto'`` is a compatibility alias for continuous, not the old router.
    ``design_policy='report'`` reports/warns about unsupported identifiability;
    'error' aborts. At least two observed donors and types are always required.
    No downstream integration or evaluation is performed.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError('adata must be an AnnData object')
    if adata.isbacked:
        raise ValueError('Use adata.to_memory() before refinement')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('AnnData cell and gene IDs must be unique; resolve duplicates explicitly')
    if mode not in ('auto', 'continuous'):
        raise ValueError('Use mode="continuous". For endpoints use fusion_config.weight_override=0 or 1.')
    if design_policy not in ('report', 'error'):
        raise ValueError('design_policy must be report or error')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a non-negative integer')
    _span_values(hvg_span)
    if query_hvg_span is not None:
        _span_values(query_hvg_span)
    if query is None and (query_hvg_genes is not None or query_hvg_span is not None):
        raise ValueError('query must be specified when using query_hvg_genes/query_hvg_span')
    if query is not None and reference is None:
        raise ValueError('Specify reference explicitly together with query')
    if output_dir is not None:
        path = Path(output_dir).expanduser()
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise FileExistsError(f'Choose a new/empty output directory: {path}')
    design_config = design_config or DesignConfig()
    fusion_config = fusion_config or ContinuousFusionConfig()
    cell_level_config = cell_level_config or CellLevelConfig(marker_log_effect=.75, marker_replication_fraction=.90)
    donor_aware_config = donor_aware_config or DonorAwareConfig()
    for obj, cls in [(design_config, DesignConfig), (fusion_config, ContinuousFusionConfig), (cell_level_config, CellLevelConfig),
                     (donor_aware_config, DonorAwareConfig)]:
        if not isinstance(obj, cls):
            raise TypeError(f'Expected {cls.__name__}, received {type(obj).__name__}')
    mask = selection_mask(adata, batch_key, label_key, reference)
    query_mask = None
    if query is not None:
        query_ids = _gene_list(query, 'query')
        missing_query = set(query_ids) - set(adata.obs[batch_key].astype(str))
        if missing_query:
            raise ValueError(f'Query donors not found: {sorted(missing_query)}')
        query_mask = adata.obs[batch_key].astype(str).isin(query_ids).to_numpy()
        if np.any(mask & query_mask):
            raise ValueError('Reference and Query donors must be disjoint')
        if not np.all(mask | query_mask):
            raise ValueError('Reference/Query do not cover all input cells; subset adata explicitly first')
    routing = audit_design(adata, batch_key=batch_key, label_key=label_key,
                           reference=reference, config=design_config, fusion_config=fusion_config)
    if routing.summary['n_donors'] < 2 or routing.summary['n_celltypes'] < 2:
        raise ValueError('Continuous evidence requires at least two observed discovery donors and cell types; inspect audit_design first.')
    if not routing.summary['identifiable']:
        message = 'Supported donor-cell-type design is not identifiable; audit is separate from A, and interpretation is limited.'
        if design_policy == 'error':
            raise ValueError(message)
        warnings.warn(message, UserWarning, stacklevel=2)
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
        panel, feature_metadata = _select_hvgs(cohort_counts, genes, cohort.obs, batch_key, n_hvg, hvg_span)
        feature_method = 'scanpy_seurat_v3_batch_aware'
    else:
        panel = _gene_list(hvg_genes, 'hvg_genes')
        if len(panel) < 2:
            raise ValueError('At least two panel genes are required')
        feature_metadata = pd.DataFrame(index=genes)
        feature_metadata.attrs['hvg_selection'] = dict(method='external_ordered_panel',
                                                      requested_span=None, used_span=None, attempts=[])
        feature_method = 'external_ordered_panel'
    missing = set(panel) - set(genes)
    if missing:
        raise ValueError(f'{len(missing)} panel genes absent from counts: {sorted(missing)[:20]}; no silent dropping/renaming')
    reference_panel = list(panel)
    hvg_selection = {'reference': feature_metadata.attrs['hvg_selection']}
    hvg_audits = {'reference': feature_metadata.reset_index(drop=True)} if hvg_genes is None else {}
    query_panel = None
    if query_mask is not None:
        if query_hvg_genes is None:
            query_panel, query_metadata = _select_hvgs(
                matrix[query_mask, :], genes, adata.obs.loc[query_mask], batch_key, n_hvg,
                hvg_span if query_hvg_span is None else query_hvg_span)
            hvg_audits['query'] = query_metadata.reset_index(drop=True)
        else:
            query_panel = _gene_list(query_hvg_genes, 'query_hvg_genes')
            if len(query_panel) < 2 or set(query_panel) - set(genes):
                raise ValueError('Query panel requires >=2 genes, all present in the counts axis')
            query_metadata = pd.DataFrame(index=genes)
            query_metadata.attrs['hvg_selection'] = dict(method='external_ordered_panel',
                                                        requested_span=None, used_span=None, attempts=[])
        hvg_selection['query'] = query_metadata.attrs['hvg_selection']
    protected = [] if protected_genes is None or len(protected_genes) == 0 else _gene_list(protected_genes, 'protected_genes')
    unknown_protected = sorted(set(protected) - set(genes))
    if unknown_protected:
        raise ValueError(f'Protected gene IDs absent from counts: {unknown_protected[:20]}')
    from ._cell_level import find_harmful_v1_celllevel
    from ._donor_aware import find_harmful_within_protocol_donor_replicate
    donor_result = find_harmful_within_protocol_donor_replicate(
        cohort, cohort_counts, genes, panel, batch_key, label_key,
        protected_genes=protected, seed=seed, risk_config=asdict(donor_aware_config))
    cell_result = find_harmful_v1_celllevel(
        cohort, cohort_counts, genes, panel, batch_key, label_key,
        protected_genes=protected, seed=seed, **asdict(cell_level_config))
    component_evidence = {'donor': donor_result['table'].copy(),
                          'cell': cell_result['table'].copy()}
    cell_gates = CellGateConfig(
        alpha=cell_level_config.alpha, leakage_z_floor=cell_level_config.leakage_z_floor,
        risk_z_floor=cell_level_config.risk_z_floor,
        biology_z_ceiling=cell_level_config.biology_z_ceiling,
        bootstrap_pass_fraction=cell_level_config.bootstrap_pass_fraction)
    donor_gates = replace(fusion_config.donor,
        technical_quantile=donor_aware_config.technical_quantile,
        risk_quantile=donor_aware_config.risk_quantile,
        minimum_absolute_risk=donor_aware_config.minimum_absolute_risk,
        maximum_biological_support_quantile=donor_aware_config.maximum_biological_support_quantile,
        minimum_direction_agreement=donor_aware_config.minimum_direction_agreement,
        minimum_single_donor_dominance=donor_aware_config.minimum_single_donor_dominance)
    if fusion_config.cell != CellGateConfig() and fusion_config.cell != cell_gates:
        raise ValueError('Conflicting gates: specify cell gates consistently in cell_level_config and fusion_config')
    if fusion_config.donor != DonorGateConfig() and fusion_config.donor != donor_gates:
        raise ValueError('Conflicting gates: specify donor gates consistently in donor_aware_config and fusion_config')
    fusion_config = replace(fusion_config, cell=cell_gates, donor=donor_gates)
    table, fusion_manifest = continuous_softgate_fusion(
        component_evidence['donor'], component_evidence['cell'],
        replication_adequacy=routing.summary['replication_adequacy'], config=fusion_config)
    table['harmful'] = table.continuous_remove
    table['explicitly_protected'] = table.gene.isin(protected)
    table['decision_reason'] = np.where(table.harmful, 'removed_continuous_evidence',
                                       'kept_continuous_evidence_at_or_below_cutoff')
    table.loc[table.explicitly_protected, 'decision_reason'] = 'kept_user_protected'
    actual = 'continuous'
    engine_id = 'continuous_softgate_v0_3_weighted_A'
    diagnostics = {'sampled_cells': cell_result['sampled_cells']}
    table['input_rank'] = table.gene.map({g: i+1 for i, g in enumerate(panel)})
    table = table.sort_values('input_rank', kind='mergesort').reset_index(drop=True)
    table['selected_mode'] = actual
    table['in_final_panel'] = table.final_action.eq('keep')
    reference_table = table.copy()
    reference_risk = reference_table.loc[reference_table.final_action.eq('remove'), 'gene'].tolist()
    if query_panel is not None:
        panel = query_panel
        feature_metadata = query_metadata
        # Preserve untested evidence as NA; Query-only genes are not zero-risk genes.
        evidence = reference_table.set_index('gene')
        for col in evidence:
            if pd.api.types.is_bool_dtype(evidence[col]):
                evidence[col] = evidence[col].astype('boolean')
        table = evidence.reindex(panel).rename_axis('gene').reset_index()
        in_reference = table.gene.isin(reference_panel)
        table['reference_hvg_rank'] = table['input_rank']
        table['input_rank'] = np.arange(1, len(panel)+1)
        table['in_reference_hvg_panel'] = in_reference
        table['reference_risk_flagged'] = table.gene.isin(reference_risk)
        table['harmful'] = table['reference_risk_flagged'] & ~table.gene.isin(protected)
        table['final_action'] = np.where(table.harmful, 'remove', 'keep')
        table['in_final_panel'] = ~table.harmful
        table['explicitly_protected'] = table.gene.isin(protected)
        table['selected_mode'] = actual
        table.loc[~in_reference, 'decision_reason'] = 'not_tested_outside_reference_hvg_panel'
        table.loc[table.explicitly_protected, 'decision_reason'] = 'kept_user_protected'
        table.loc[table.harmful, 'decision_reason'] = 'removed_reference_risk_intersection_with_query_hvg'
        for col in table.select_dtypes(include=['object', 'string']).columns:
            # h5ad string columns cannot mix missing float sentinels and strings.
            table[col] = table[col].fillna('').astype(str)
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
                     for name in ['_cell_level.py', '_donor_aware.py', 'routing.py', 'design.py',
                                  '_continuous_core.py', 'fusion.py', 'continuous_workflow.py']}
    cohort_ids = adata.obs_names[mask].astype(str).tolist()
    manifest = dict(
        schema_version='0.13.0', software_versions=versions, engine_sha256=engine_hashes,
        engine_id=engine_id, requested_mode=mode, selected_mode=actual, design_audit=routing.summary,
        route_selection_used=False, replication_adequacy_source='computed_from_discovery_metadata',
        replication_formula='weighted_harmonic_v1', fusion=fusion_manifest,
        design_policy=design_policy,
        feature_method=hvg_selection['query']['method'] if query_panel is not None else feature_method,
        reference_feature_method=hvg_selection['reference']['method'],
        base_n_hvg=len(panel), final_n_hvg=len(retained),
        removed_genes=removed, base_panel_sha256=gene_hash(panel), final_panel_sha256=gene_hash(retained),
        batch_key=batch_key, label_key=label_key, selection_donors=sorted(cohort.obs[batch_key].astype(str).unique()),
        selection_cell_ids_sha256=gene_hash(cohort_ids), n_selection_cells=len(cohort_ids),
        selection_scope='all_input_cells' if reference is None else 'specified_reference_donors',
        selection_uses_celltype_labels=True, downstream_evaluation_performed=False,
        raw_counts_location=source.location, counts_numeric_validation=numeric,
        raw_provenance_note='Integer-like counts do not establish experimental raw-count provenance.',
        protected_genes=protected, protected_genes_outside_panel=sorted(set(protected)-set(panel)),
        design_config=asdict(design_config), cell_level_config=asdict(cell_level_config),
        donor_aware_config=asdict(donor_aware_config), seed=seed, diagnostics=diagnostics,
        workflow='reference_risk_to_query_panel' if query_panel is not None else 'reference_panel_refinement',
        hvg_selection=hvg_selection, reference_base_genes=reference_panel,
        query_base_genes=query_panel or [], reference_risk_genes=reference_risk,
        reference_base_panel_sha256=gene_hash(reference_panel),
        reference_risk_n=len(reference_risk),
        query_donors=sorted(adata.obs.loc[query_mask, batch_key].astype(str).unique()) if query_mask is not None else [],
        n_query_cells=int(query_mask.sum()) if query_mask is not None else 0,
        query_cell_ids_sha256=gene_hash(adata.obs_names[query_mask].astype(str).tolist()) if query_mask is not None else '',
        query_expression_used_for_hvg_selection=query_panel is not None and query_hvg_genes is None,
        query_labels_used_for_selection_or_risk=False,
        risk_evidence_scope=('reference_hvg_panel_on_reference_cells' if reference is not None
                             else 'pooled_hvg_panel_on_all_input_cells'),
        output_cells='all_input_cells',
    )
    final.uns['hvgdecision'] = dict(
        manifest_json=json.dumps(manifest, ensure_ascii=False),
        gene_decisions=table.copy(), routing_summary=routing.audit,
        removed_genes=np.asarray(removed, dtype=str), base_genes=np.asarray(panel, dtype=str),
        selection_cell_ids=np.asarray(cohort_ids, dtype=str), selected_mode=actual,
        reference_gene_decisions=reference_table.copy(),
        reference_risk_genes=np.asarray(reference_risk, dtype=str),
    )
    updated_counts_audit = source.audit.copy()
    for col, val in numeric.items():
        if col == 'source':
            continue
        selected_rows = updated_counts_audit['selected'].fillna(False).astype(bool)
        updated_counts_audit.loc[selected_rows, col] = val
    output = RefinementResult(final, table, routing, updated_counts_audit, manifest,
                              reference_table, hvg_audits, component_evidence)
    if output_dir is not None:
        output.save(output_dir)
    return output if return_details else final
