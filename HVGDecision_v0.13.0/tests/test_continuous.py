import json
from dataclasses import replace

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
    x[:, -1] = 1
    x[k != 0, -1] += 2
    x[(d == 0) & (k == 0), -1] += 60
    obs = pd.DataFrame({'donor': [f'd{i}' for i in d], 'type': [f't{i}' for i in k]},
                       index=[f'c{i}' for i in range(len(d))])
    return ad.AnnData(sparse.csr_matrix(x) if sparse_counts else x, obs=obs,
                      var=pd.DataFrame(index=[f'g{i}' for i in range(genes)]))


def run(a, **kwargs):
    args = dict(batch_key='donor', label_key='type', hvg_genes=a.var_names.tolist(),
                cell_level_config=hd.CellLevelConfig(n_permutations=9, n_bootstraps=3),
                return_details=True)
    args.update(kwargs)
    return hd.refine(a, **args)


def test_weighted_formula_not_unweighted():
    a = example(donors=3)
    a = a[(a.obs.type != 't2') | (np.arange(a.n_obs) % 20 == 0)].copy()
    r = hd.audit_design(a, batch_key='donor', label_key='type')
    n = pd.crosstab(a.obs.donor, a.obs.type).to_numpy(float)
    s = n/(n+15); v = np.minimum(1, n.sum(0)/(len(n)*15)); rk = s.sum(0)
    u = s.mean(1); b = u.sum()**2/(len(n)*np.square(u).sum())
    weighted = v.sum()/np.sum(v/(rk+1e-12))*np.sqrt(b)
    unweighted = len(rk)/np.sum(1/(rk+1e-12))*np.sqrt(b)
    assert r.summary['replication_adequacy'] == pytest.approx(weighted)
    assert weighted != pytest.approx(unweighted)
    assert r.summary['donor_aware_weight'] == pytest.approx(1/(1+np.exp(-(weighted-3)/.5)))
    assert r.summary['replication_formula'] == 'weighted_harmonic_v1'
    assert r.mode == 'continuous'


@pytest.mark.parametrize('sparse_counts', [True, False])
def test_roundtrip_and_no_mutation(tmp_path, sparse_counts):
    a = example(sparse_counts=sparse_counts)
    before = a.X.copy()
    a.obsm['X_pca'] = np.zeros((a.n_obs, 2))
    r = run(a, output_dir=tmp_path/'run')
    delta = a.X-before
    assert delta.nnz == 0 if sparse_counts else np.array_equal(a.X, before)
    assert r.selected_mode == 'continuous' and not r.adata.obsm
    assert r.final_n_hvg + len(r.removed_genes) == a.n_vars
    assert (tmp_path/'run/donor_evidence.csv').exists()
    assert (tmp_path/'run/cell_evidence.csv').exists()
    b = ad.read_h5ad(tmp_path/'run/adata_hvg.h5ad')
    assert b.shape == r.adata.shape
    manifest = json.loads((tmp_path/'run/run_manifest.json').read_text())
    assert manifest['route_selection_used'] is False
    assert manifest['replication_adequacy_source'] == 'computed_from_discovery_metadata'
    assert manifest['cell_level_config']['marker_log_effect'] == .75
    assert manifest['cell_level_config']['marker_replication_fraction'] == .90
    with pytest.raises(FileExistsError):
        r.save(tmp_path/'run')


def test_both_engines_and_endpoints():
    a = example()
    r = run(a)
    for w, col in [(0, 'cell_hard_remove'), (1, 'donor_hard_remove')]:
        cfg = hd.ContinuousFusionConfig(weight_override=w)
        t, m = hd.continuous_softgate_fusion(r.component_evidence['donor'],
            r.component_evidence['cell'], replication_adequacy=r.routing.summary['replication_adequacy'], config=cfg)
        assert np.array_equal(t.continuous_remove, t[col])
        assert m['donor_endpoint_exact_at_w1'] and m['cell_endpoint_exact_at_w0']
    assert set(r.component_evidence) == {'donor', 'cell'}
    t, _ = hd.continuous_softgate_fusion(r.component_evidence['donor'],
        r.component_evidence['cell'], replication_adequacy=r.routing.summary['replication_adequacy'],
        config=hd.ContinuousFusionConfig())
    assert set(t.loc[t.continuous_remove, 'gene']) == set(r.removed_genes)


def test_query_labels_and_expression_not_used_for_risk(tmp_path):
    a = example(donors=8)
    args = dict(reference=[f'd{i}' for i in range(6)], query=['d6','d7'],
                hvg_genes=a.var_names[1:].tolist(), query_hvg_genes=a.var_names.tolist())
    r = run(a, **args, output_dir=tmp_path/'query')
    row = r.decision_table.set_index('gene').loc['g0']
    assert row.final_action == 'keep' and pd.isna(row.unified_harmfulness)
    assert row.decision_reason == 'not_tested_outside_reference_hvg_panel'
    q = a.obs.donor.isin(args['query'])
    a.obs.loc[q, 'type'] = None
    a.X[q.to_numpy(), :] *= 3
    r2 = run(a, **args)
    pd.testing.assert_frame_equal(r.reference_decision_table, r2.reference_decision_table)
    assert r.routing.summary == r2.routing.summary
    assert set(r.removed_genes) == set(r.reference_risk_genes) & set(args['query_hvg_genes'])
    assert r.adata.n_obs == a.n_obs


def test_explicit_protection_and_zero_removal():
    a = example()
    r = run(a, protected_genes=a.var_names.tolist())
    assert not r.removed_genes and r.final_n_hvg == a.n_vars
    assert (r.reference_decision_table.unified_harmfulness == 0).all()


def test_invalid_inputs_and_no_legacy_route():
    a = example()
    with pytest.raises(ValueError, match='mode='):
        run(a, mode='donor_aware')
    with pytest.raises(ValueError, match='disjoint'):
        run(a, reference=['d0'], query=['d0'])
    with pytest.raises(ValueError, match='absent from counts'):
        run(a, hvg_genes=['not_a_gene','g0'])
    with pytest.raises(ValueError, match='two observed'):
        run(example(types=1))
    with pytest.raises(ValueError):
        hd.ContinuousFusionConfig(weight_scale=float('nan'))
    with pytest.raises(ValueError):
        hd.ContinuousFusionConfig(weight_override=float('nan'))


def test_disconnected_audit_not_a_branch_switch():
    a = example(donors=4, types=4)
    a = a[a.obs.donor.str[1:] == a.obs.type.str[1:]].copy()
    r = hd.audit_design(a, batch_key='donor', label_key='type')
    assert not r.summary['identifiable'] and r.summary['connected_components'] == 4
    assert r.mode == 'continuous' and r.summary['route_selection_used'] is False
    with pytest.raises(ValueError, match='not identifiable'):
        run(a, design_policy='error')


def test_boundary_signs():
    from hvgdecision._continuous_core import _signed_margin
    m = _signed_margin([0.,0.], [True,False], scale=1, epsilon=1e-9)
    assert m[0] > 0 and m[1] < 0


def test_explicit_and_external_counts(tmp_path):
    a = example(sparse_counts=False)
    frame = pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names)
    frame.iloc[::-1].to_csv(tmp_path/'counts.csv')
    source = hd.find_raw_counts(a, source=tmp_path/'counts.csv')
    assert source.valid
    np.testing.assert_array_equal(source.matrix, a.X)
    a.layers['counts'] = a.X.copy()
    assert hd.find_raw_counts(a, source='counts').valid
    a.X = np.log1p(a.X)
    assert not hd.find_raw_counts(a, source='X').valid
    assert hd.find_raw_counts(a).valid
    r = run(a, counts=source)
    assert r.adata.n_obs == a.n_obs


def test_default_permutations_and_repeatability():
    a = example(donors=3, cells=12, genes=48)
    r1 = run(a, cell_level_config=hd.CellLevelConfig())
    r2 = run(a, cell_level_config=hd.CellLevelConfig())
    assert r1.run_info['cell_level_config']['n_permutations'] == 100
    assert r1.run_info['cell_level_config']['n_bootstraps'] == 20
    pd.testing.assert_frame_equal(r1.decision_table, r2.decision_table)


def test_independent_seurat_panels_and_fixed_spans():
    import scanpy as sc
    a = example(donors=8, genes=600)
    r = run(a, hvg_genes=None, n_hvg=200, hvg_span=.5, query_hvg_span=.3,
            reference=[f'd{i}' for i in range(6)], query=['d6','d7'])
    for role, donors, span in [('reference', [f'd{i}' for i in range(6)], .5), ('query',['d6','d7'],.3)]:
        b = a[a.obs.donor.isin(donors)].copy()
        sc.pp.highly_variable_genes(b, flavor='seurat_v3', n_top_genes=200, batch_key='donor', span=span)
        ordered = b.var.sort_values(['highly_variable','highly_variable_rank'], ascending=[False,True],kind='mergesort')
        assert r.run_info[role+'_base_genes'] == ordered.index[ordered.highly_variable].tolist()
        assert r.run_info['hvg_selection'][role]['used_span'] == span
