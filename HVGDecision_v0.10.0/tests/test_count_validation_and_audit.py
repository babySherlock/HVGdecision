import numpy as np
import pandas as pd
import pytest
import hvgdecision as hd
from hvgdecision.scoring import _minmax
from test_current_workflow import example


def test_missing_metric_is_not_neutral_score():
    result = _minmax(pd.Series([0.8, np.nan, 0.8]))
    assert result.iloc[0] == 0.5
    assert pd.isna(result.iloc[1])


def test_full_gate_audit_and_existing_annotation():
    a = example()
    a.var['gene_symbol'] = ['symbol_' + g for g in a.var_names]
    r = hd.refine(a, batch_key='donor', label_key='type', hvg_genes=a.var_names.tolist(), return_details=True)
    assert r.decision_table.gene_symbol.tolist() == a.var.gene_symbol.tolist()
    assert r.decision_table.loc[r.decision_table.harmful, 'failed_rules'].eq('').all()
    assert r.decision_table.loc[~r.decision_table.harmful, 'failed_rules'].str.len().gt(0).all()


def test_fractional_count_outside_old_sample_tolerance_rejected():
    a = example(sparse_counts=False)
    a.X[0, 0] = 100000.1  # relative tolerance would incorrectly accept large values
    assert not hd.find_raw_counts(a).valid


def test_no_supported_edges_and_minimum_library_qc():
    a = example(cells=1)
    design = hd.audit_design(a, batch_key='donor', label_key='type')
    assert design.mode == 'insufficient_confounded'
    assert design.summary['design_rank'] == 0
    a = example(sparse_counts=False)
    a.X[0, :] = 0
    with pytest.raises(ValueError, match='zero-library'):
        hd.refine(a, batch_key='donor', label_key='type', hvg_genes=a.var_names.tolist())
