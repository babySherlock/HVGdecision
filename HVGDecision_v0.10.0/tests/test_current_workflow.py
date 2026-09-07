import json
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import hvgdecision as hd


def example(donors=6, types=3, cells=24, genes=96, sparse_counts=True):
    d = np.repeat(np.arange(donors), types*cells)
    k = np.tile(np.repeat(np.arange(types), cells), donors)
    rng = np.random.default_rng(846)
    x = rng.poisson(2.0, (len(d), genes)).astype(np.float32)
    for i in range(types):
        x[k == i, i*5:(i+1)*5] += 5
    # A donor-specific cell-type interaction, not a named real gene.
    x[:, -1] = 1
    x[k != 0, -1] += 2
    x[(d == 0) & (k == 0), -1] += 60
    obs = pd.DataFrame({'donor': [f'd{i}' for i in d], 'type': [f't{i}' for i in k]},
                       index=[f'c{i}' for i in range(len(d))])
    return ad.AnnData(sparse.csr_matrix(x) if sparse_counts else x, obs=obs,
                      var=pd.DataFrame(index=[f'g{i}' for i in range(genes)]))


def test_replication_formula_and_routing():
    a = example()
    result = hd.audit_design(a, batch_key='donor', label_key='type')
    expected_a = 6 * 24 / (24 + 15)
    assert result.mode == 'donor_aware'
    assert result.summary['replication_adequacy'] == pytest.approx(expected_a)
    assert result.summary['donor_aware_weight'] == pytest.approx((expected_a-2)/2)
    assert result.summary['design_rank'] == 8
    assert result.summary['identifiability'] == 1
    assert hd.audit_design(example(donors=3), batch_key='donor', label_key='type').mode == 'cell_level'


@pytest.mark.parametrize('same_type', [True, False])
def test_hundred_donors_one_type_each_is_not_sufficient(same_type):
    a = example(donors=100, types=1, cells=20)
    if not same_type:
        a.obs['type'] = a.obs.donor.to_numpy()
    r = hd.audit_design(a, batch_key='donor', label_key='type')
    assert r.mode == 'insufficient_confounded'
    if not same_type:
        assert r.summary['connected_components'] == 100


def test_weak_bridge_does_not_establish_identifiability():
    a = example(donors=2, types=2, cells=20)
    keep = ((a.obs.donor == 'd0') & (a.obs.type == 't0')) | ((a.obs.donor == 'd1') & (a.obs.type == 't1'))
    keep.iloc[20] = True  # one bridge cell < tau
    r = hd.audit_design(a[keep].copy(), batch_key='donor', label_key='type')
    assert r.summary['connected_components'] == 2
    assert r.mode == 'insufficient_confounded'


def test_confounded_no_deletion_and_roundtrip(tmp_path):
    a = example(donors=3, types=1)
    with pytest.warns(UserWarning, match='insufficient'):
        result = hd.refine(a, batch_key='donor', label_key='type',
                           hvg_genes=list(a.var_names), return_details=True,
                           mode='donor_aware', output_dir=tmp_path / 'run')
    assert not result.removed_genes
    assert result.selected_mode == 'insufficient_confounded'
    loaded = ad.read_h5ad(tmp_path / 'run/adata_hvg.h5ad')
    assert loaded.shape == a.shape
    assert len(loaded.uns['hvgdecision']['gene_decisions']) == a.n_vars
    assert (loaded.X != a.X).nnz == 0
    with pytest.raises(FileExistsError):
        result.save(tmp_path / 'run')


@pytest.mark.parametrize('mode', ['cell_level', 'donor_aware'])
@pytest.mark.parametrize('sparse_counts', [True, False])
def test_engines_return_auditable_adata_without_mutation(mode, sparse_counts, tmp_path):
    a = example(sparse_counts=sparse_counts)
    before = a.X.copy()
    a.obsm['X_pca'] = np.zeros((a.n_obs, 2))
    kwargs = dict(batch_key='donor', label_key='type', hvg_genes=list(a.var_names),
                  return_details=True, mode=mode, cell_level_config=hd.CellLevelConfig(n_permutations=9, n_bootstraps=3))
    with pytest.warns(UserWarning, match='override'):
        result = hd.refine(a, **kwargs)
    delta = a.X - before
    assert (delta.nnz == 0) if sparse_counts else np.array_equal(a.X, before)
    assert not result.adata.obsm
    assert result.base_n_hvg == a.n_vars
    assert set(result.adata.var_names).isdisjoint(result.removed_genes)
    assert result.final_n_hvg + len(result.removed_genes) == result.base_n_hvg
    result.save(tmp_path / mode)
    loaded = ad.read_h5ad(tmp_path / mode / 'adata_hvg.h5ad')
    pd.testing.assert_frame_equal(loaded.uns['hvgdecision']['gene_decisions'], result.decision_table, check_categorical=False)
    info = json.loads(loaded.uns['hvgdecision']['manifest_json'])
    assert info['downstream_evaluation_performed'] is False


def test_donor_positive_removal_and_protection():
    a = example()
    args = dict(batch_key='donor', label_key='type', hvg_genes=a.var_names.tolist(), return_details=True)
    r = hd.refine(a, **args)
    assert r.removed_genes, 'Synthetic positive control must exercise actual removal'
    p = hd.refine(a, **args, protected_genes=r.removed_genes)
    assert set(r.removed_genes).issubset(p.adata.var_names)
    assert not p.removed_genes


@pytest.mark.parametrize('bad', [np.nan, np.inf, -1, 0.01])
def test_counts_exhaustive_not_sampled(bad):
    a = example(sparse_counts=False)
    a.X[0, 0] = bad
    result = hd.find_raw_counts(a)
    assert not result.valid


def test_raw_axis_and_external_alignment():
    full = example()
    a = full[:, :30].copy()
    a.raw = full
    a.X = a.X.astype(float) / 7
    source = hd.find_raw_counts(a)
    assert source.source == 'raw'
    r = hd.refine(a, batch_key='donor', label_key='type', counts=source,
                  hvg_genes=full.var_names.tolist(), return_details=True)
    assert r.base_n_hvg == 96
    frame = pd.DataFrame(full.X.toarray(), index=full.obs_names, columns=full.var_names)
    external = hd.find_raw_counts(a, source=frame.iloc[::-1].T)
    assert external.valid
    assert np.array_equal(external.matrix, full.X.toarray())
    frame.index = ['wrong'+str(i) for i in range(len(frame))]
    assert not hd.find_raw_counts(a, source=frame).valid


def test_invalid_metadata_missing_genes_and_stale_counts():
    a = example()
    source = hd.find_raw_counts(a)
    a.X.data[0] = np.nan
    with pytest.raises(ValueError, match='non-finite'):
        hd.refine(a, batch_key='donor', label_key='type', counts=source, hvg_genes=a.var_names)
    a = example()
    with pytest.raises(ValueError, match='absent'):
        hd.refine(a, batch_key='donor', label_key='type', hvg_genes=['g1', 'missing'])
    a.obs.iloc[0, 0] = None
    with pytest.raises(ValueError, match='missing metadata'):
        hd.audit_design(a, batch_key='donor', label_key='type')


def test_reference_isolation():
    a = example()
    args = dict(batch_key='donor', label_key='type', reference=['d0', 'd1', 'd2'],
                hvg_genes=a.var_names.tolist(), return_details=True,
                cell_level_config=hd.CellLevelConfig(n_permutations=9, n_bootstraps=3))
    r1 = hd.refine(a, **args)
    altered = a.copy()
    altered.X = altered.X.toarray()
    mask = ~altered.obs.donor.isin(args['reference'])
    altered.X[mask] *= 9
    altered.obs.loc[mask, 'type'] = 'unseen_query_annotation'
    r2 = hd.refine(altered, **args)
    pd.testing.assert_frame_equal(r1.decision_table, r2.decision_table)
    assert r1.adata.n_obs == a.n_obs


def test_external_cached_source_rejects_reordered_cells():
    a = example()
    source = hd.find_raw_counts(a, source=(a.X, a.var_names))
    with pytest.raises(ValueError, match='cell axis/order'):
        hd.refine(a[::-1].copy(), batch_key='donor', label_key='type', counts=source,
                  hvg_genes=a.var_names.tolist())


def test_default_2000_and_downstream_preprocessing(tmp_path):
    import scanpy as sc
    a = example(genes=2500)
    r = hd.refine(a, batch_key='donor', label_key='type', return_details=True,
                  output_dir=tmp_path / 'default2000')
    assert r.base_n_hvg == 2000
    assert r.final_n_hvg == 2000 - len(r.removed_genes)
    work = ad.read_h5ad(tmp_path / 'default2000/adata_hvg.h5ad')
    before = work.layers['counts'].copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.pca(work, n_comps=20)
    sc.pp.neighbors(work, n_neighbors=10, n_pcs=20)
    assert work.obsm['X_pca'].shape == (a.n_obs, 20)
    assert (work.layers['counts'] != before).nnz == 0


def test_real_seurat_selection_and_cli(tmp_path):
    a = example(genes=500)
    r = hd.refine(a, batch_key='donor', label_key='type', n_hvg=100, return_details=True)
    assert r.base_n_hvg == 100
    assert r.run_info['feature_method'] == 'scanpy_seurat_v3_batch_aware'
    a.write_h5ad(tmp_path / 'input.h5ad')
    pd.DataFrame({'gene': list(a.var_names[:100])}).to_csv(tmp_path / 'panel.csv', index=False)
    proc = subprocess.run([sys.executable, '-m', 'hvgdecision.cli', 'refine',
                           '--input', str(tmp_path / 'input.h5ad'), '--batch-key', 'donor',
                           '--label-key', 'type', '--hvg-table', str(tmp_path / 'panel.csv'),
                           '--output', str(tmp_path / 'cli')], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / 'cli/adata_hvg.h5ad').exists()


@pytest.mark.parametrize('constructor,kwargs', [
    (hd.RoutingConfig, {'tau': 0}), (hd.RoutingConfig, {'donor_aware_threshold': 2}),
    (hd.CellLevelConfig, {'n_permutations': 0}), (hd.CellLevelConfig, {'alpha': -1}),
    (hd.DonorAwareConfig, {'min_group_cells': 0}), (hd.DonorAwareConfig, {'risk_quantile': 2}),
])
def test_parameters(constructor, kwargs):
    with pytest.raises(ValueError):
        constructor(**kwargs)
