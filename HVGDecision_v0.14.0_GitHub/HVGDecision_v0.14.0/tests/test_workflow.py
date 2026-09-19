import json

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
scanpy = pytest.importorskip("scanpy")
from scipy import sparse

import hvgdecision as hd



def _example(donors=6, types=3, cells=24, genes=96, sparse_counts=True):
    d = np.repeat(np.arange(donors), types * cells)
    k = np.tile(np.repeat(np.arange(types), cells), donors)
    rng = np.random.default_rng(846)
    x = rng.poisson(2.0, (len(d), genes)).astype(np.float32)
    for i in range(types):
        x[k == i, i * 5 : (i + 1) * 5] += 5
    x[:, -1] = 1
    x[k != 0, -1] += 2
    x[(d == 0) & (k == 0), -1] += 60
    obs = pd.DataFrame(
        {
            "donor": [f"d{i}" for i in d],
            "type": [f"t{i}" for i in k],
        },
        index=[f"c{i}" for i in range(len(d))],
    )
    return anndata.AnnData(
        sparse.csr_matrix(x) if sparse_counts else x,
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(genes)]),
    )


def _run(a, **kwargs):
    args = dict(
        batch_key="donor",
        label_key="type",
        hvg_genes=a.var_names.tolist(),
        cell_level_config=hd.CellLevelConfig(n_permutations=9, n_bootstraps=3),
        fusion_config=hd.UnifiedFusionConfig(tau=0.10),
        return_details=True,
    )
    args.update(kwargs)
    return hd.refine(a, **args)


@pytest.mark.parametrize("sparse_counts", [True, False])
def test_roundtrip_no_mutation_and_manifest(tmp_path, sparse_counts):
    a = _example(sparse_counts=sparse_counts)
    before = a.X.copy()
    a.obsm["X_pca"] = np.zeros((a.n_obs, 2))
    r = _run(a, output_dir=tmp_path / "run")

    if sparse_counts:
        assert (a.X - before).nnz == 0
    else:
        assert np.array_equal(a.X, before)

    assert r.selected_mode == "dual_evidence"
    assert not r.adata.obsm
    assert r.final_n_hvg + len(r.removed_genes) == a.n_vars
    assert (tmp_path / "run/donor_evidence.csv").exists()
    assert (tmp_path / "run/cell_evidence.csv").exists()
    assert (tmp_path / "run/design_audit.csv").exists()

    manifest = json.loads((tmp_path / "run/run_manifest.json").read_text())
    assert manifest["component_selection_used"] is False
    assert manifest["dataset_dependent_component_weight_used"] is False
    assert manifest["fusion"]["method"] == "dual_evidence_signed_margin_centered_logmeanexp_v0_4"
    assert manifest["fusion"]["tau"] == pytest.approx(0.10)

    with pytest.raises(FileExistsError):
        r.save(tmp_path / "run")


def test_query_labels_not_used_for_reference_risk():
    a = _example(donors=8)
    args = dict(
        reference=[f"d{i}" for i in range(6)],
        query=["d6", "d7"],
        hvg_genes=a.var_names[1:].tolist(),
        query_hvg_genes=a.var_names.tolist(),
    )
    r1 = _run(a, **args)

    q = a.obs.donor.isin(args["query"])
    a.obs.loc[q, "type"] = "QUERY_LABEL_CHANGED"
    a.X[q.to_numpy(), :] *= 3
    r2 = _run(a, **args)

    pd.testing.assert_frame_equal(r1.reference_decision_table, r2.reference_decision_table)
    assert r1.design.summary == r2.design.summary
    assert set(r1.removed_genes) == set(r1.reference_risk_genes) & set(args["query_hvg_genes"])


def test_explicit_protection_allows_zero_removal():
    a = _example()
    r = _run(a, protected_genes=a.var_names.tolist())
    assert not r.removed_genes
    assert r.final_n_hvg == a.n_vars
    assert not r.reference_decision_table["discovery_remove"].any()


def test_repeatability():
    a = _example(donors=3, cells=12, genes=48)
    cfg = hd.CellLevelConfig(n_permutations=9, n_bootstraps=3)
    r1 = _run(a, cell_level_config=cfg)
    r2 = _run(a, cell_level_config=cfg)
    pd.testing.assert_frame_equal(r1.reference_decision_table, r2.reference_decision_table)
