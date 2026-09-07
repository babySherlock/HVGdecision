"""Frozen statistical engine; see docs/ENGINE_PROVENANCE.json."""
import numpy as np
import pandas as pd
from scipy import sparse
import warnings

def _to_numpy_1d(x):
    return np.asarray(x).reshape(-1)

def robust_z(x, eps=1e-8):
    """
    Median/MAD robust z-score across genes.
    """
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    scale = 1.4826 * mad

    if not np.isfinite(scale) or scale < eps:
        scale = np.nanstd(x)

    if not np.isfinite(scale) or scale < eps:
        scale = 1.0

    return (x - med) / scale

def bh_fdr(p_values):
    """
    Benjamini-Hochberg FDR.
    """
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)

    valid = np.isfinite(p)
    if not valid.any():
        return q

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]

    m = len(ranked)
    raw_q = ranked * m / np.arange(1, m + 1)
    raw_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    raw_q = np.clip(raw_q, 0.0, 1.0)

    tmp = np.empty_like(raw_q)
    tmp[order] = raw_q
    q[valid] = tmp
    return q

def _group_indices(values):
    values = np.asarray(values).astype(str)
    groups = {}
    for i, value in enumerate(values):
        groups.setdefault(value, []).append(i)
    return {
        k: np.asarray(v, dtype=int)
        for k, v in groups.items()
    }

def _normalize_selected_cells(
    counts,
    gene_names,
    selected_genes,
    cell_indices=None,
    target_sum=1e4,
):
    """
    Normalize selected genes using TOTAL raw counts across all genes
    as library-size denominator. This avoids re-normalizing only within
    the 2000-HVG panel.
    """
    gene_names = pd.Index(gene_names.astype(str))
    selected_genes = [str(g) for g in selected_genes]

    gene_index = gene_names.get_indexer(selected_genes)
    if np.any(gene_index < 0):
        missing = [
            selected_genes[i]
            for i, idx in enumerate(gene_index)
            if idx < 0
        ]
        raise KeyError(
            f"{len(missing)} selected genes are absent from raw counts. "
            f"Examples: {missing[:10]}"
        )

    if cell_indices is None:
        cell_indices = np.arange(counts.shape[0], dtype=int)
    else:
        cell_indices = np.asarray(cell_indices, dtype=int)

    full = counts[cell_indices, :]
    selected = full[:, gene_index]

    total = _to_numpy_1d(full.sum(axis=1)).astype(float)
    total[~np.isfinite(total)] = 0.0
    total[total <= 0] = 1.0

    sf = target_sum / total

    if sparse.issparse(selected):
        X = selected.tocsr().astype(np.float64)
        X = sparse.diags(sf) @ X
        X.data = np.log1p(X.data)
        return X

    X = np.asarray(selected, dtype=np.float64)
    X *= sf[:, None]
    np.log1p(X, out=X)
    return X

def _biology_eta2(X, labels):
    """
    Fraction of cell-level expression variance explained by biology label.
    """
    labels = np.asarray(labels).astype(str)
    Xd = X.toarray() if sparse.issparse(X) else np.asarray(X)

    global_mean = Xd.mean(axis=0)
    total_ss = np.sum(
        (Xd - global_mean) ** 2,
        axis=0,
    )

    between_ss = np.zeros(Xd.shape[1], dtype=float)

    for _, idx in _group_indices(labels).items():
        if len(idx) == 0:
            continue
        mu = Xd[idx].mean(axis=0)
        between_ss += len(idx) * (mu - global_mean) ** 2

    return between_ss / (total_ss + 1e-12)

def _residualize_by_label(X, labels):
    labels = np.asarray(labels).astype(str)
    Xd = X.toarray() if sparse.issparse(X) else np.asarray(X)
    R = Xd.copy()

    for _, idx in _group_indices(labels).items():
        R[idx] -= Xd[idx].mean(axis=0)

    return R

def _donor_leakage_from_label_residual(R, batches):
    """
    Fraction of label-residual variance explained by donor identity.
    """
    batches = np.asarray(batches).astype(str)

    total_ss = np.sum(R ** 2, axis=0)
    donor_ss = np.zeros(R.shape[1], dtype=float)

    for _, idx in _group_indices(batches).items():
        if len(idx) == 0:
            continue
        mu = R[idx].mean(axis=0)
        donor_ss += len(idx) * (mu ** 2)

    return donor_ss / (total_ss + 1e-12)

def _celllevel_interaction_range(X, labels, batches):
    """
    Baseline interaction metric:
      within each biology label,
      donor mean max-min divided by within-label SD,
      then median across eligible labels.
    This intentionally represents the older / less robust baseline.
    """
    labels = np.asarray(labels).astype(str)
    batches = np.asarray(batches).astype(str)
    Xd = X.toarray() if sparse.issparse(X) else np.asarray(X)

    metrics = []

    for label, label_idx in _group_indices(labels).items():
        donor_values = batches[label_idx]
        donor_groups = _group_indices(donor_values)

        if len(donor_groups) < 2:
            continue

        donor_means = []
        for _, local_idx in donor_groups.items():
            idx = label_idx[local_idx]
            donor_means.append(
                Xd[idx].mean(axis=0)
            )

        donor_means = np.vstack(donor_means)
        donor_range = (
            np.max(donor_means, axis=0)
            - np.min(donor_means, axis=0)
        )

        scale = np.std(
            Xd[label_idx],
            axis=0,
            ddof=1,
        )
        scale[~np.isfinite(scale)] = 0.0

        metrics.append(
            donor_range / (scale + 1e-6)
        )

    if not metrics:
        return np.zeros(Xd.shape[1], dtype=float)

    return np.nanmedian(
        np.vstack(metrics),
        axis=0,
    )

def _stratified_sample_indices(
    obs,
    batch_key,
    label_key,
    max_cells_per_stratum,
    seed,
):
    rng = np.random.default_rng(seed)

    batch = obs[batch_key].astype(str).to_numpy()
    label = obs[label_key].astype(str).to_numpy()

    strata = {}
    for i, (b, l) in enumerate(zip(batch, label)):
        strata.setdefault((b, l), []).append(i)

    selected = []

    for key in sorted(strata):
        idx = np.asarray(strata[key], dtype=int)

        if len(idx) > max_cells_per_stratum:
            idx = rng.choice(
                idx,
                size=max_cells_per_stratum,
                replace=False,
            )

        selected.extend(idx.tolist())

    return np.asarray(
        sorted(selected),
        dtype=int,
    )

def _marker_protection_from_groups(
    group_expr,
    group_batches,
    group_labels,
    group_n_cells,
    genes,
    log_effect=0.50,
    replication_fraction=0.80,
    min_eligible_donors=2,
    min_cells_per_side=10,
):
    """
    Replicated positive-marker protection.

    Denominator = ELIGIBLE donors for each target label, not all donors.
    A donor is eligible when:
      target-label cells >= min_cells_per_side
      other-label cells in same donor >= min_cells_per_side
    """
    genes = np.asarray(genes).astype(str)
    group_batches = np.asarray(group_batches).astype(str)
    group_labels = np.asarray(group_labels).astype(str)
    group_n_cells = np.asarray(group_n_cells, dtype=int)
    Y = np.asarray(group_expr, dtype=float)

    protected = np.zeros(len(genes), dtype=bool)
    best_label = np.array([""] * len(genes), dtype=object)
    best_fraction = np.zeros(len(genes), dtype=float)
    best_n_eligible = np.zeros(len(genes), dtype=int)

    batches = sorted(pd.unique(group_batches))
    labels = sorted(pd.unique(group_labels))

    for label in labels:
        effects = []
        eligible_donors = []

        for batch in batches:
            target_idx = np.where(
                (group_batches == batch)
                & (group_labels == label)
            )[0]

            other_idx = np.where(
                (group_batches == batch)
                & (group_labels != label)
            )[0]

            if len(target_idx) == 0 or len(other_idx) == 0:
                continue

            target_cells = int(
                group_n_cells[target_idx].sum()
            )
            other_cells = int(
                group_n_cells[other_idx].sum()
            )

            if (
                target_cells < min_cells_per_side
                or other_cells < min_cells_per_side
            ):
                continue

            # Usually exactly one target group per donor×label.
            target_mean = np.average(
                Y[target_idx],
                axis=0,
                weights=group_n_cells[target_idx],
            )

            other_mean = np.average(
                Y[other_idx],
                axis=0,
                weights=group_n_cells[other_idx],
            )

            effects.append(target_mean - other_mean)
            eligible_donors.append(batch)

        n_eligible = len(eligible_donors)

        if n_eligible < min_eligible_donors:
            continue

        effects = np.vstack(effects)
        pass_fraction = np.mean(
            effects >= log_effect,
            axis=0,
        )

        positive = (
            pass_fraction >= replication_fraction
        )

        update = positive & (
            (pass_fraction > best_fraction)
            | (
                (pass_fraction == best_fraction)
                & (n_eligible > best_n_eligible)
            )
        )

        protected |= positive
        best_label[update] = label
        best_fraction[update] = pass_fraction[update]
        best_n_eligible[update] = n_eligible

    return pd.DataFrame(
        {
            "gene": genes,
            "marker_protected": protected,
            "marker_protected_label": best_label,
            "marker_replication_fraction": best_fraction,
            "marker_n_eligible_donors": best_n_eligible,
        }
    )

def build_donor_label_aggregate(
    adata,
    counts,
    gene_names,
    selected_genes,
    batch_key,
    label_key,
):
    """
    Build donor×biology pseudoaggregate expression from raw counts.

    For each donor×label group:
      - raw counts for selected genes are SUMMED
      - library-size denominator uses SUM of full-transcriptome raw counts
      - normalized to 1e6 and log1p transformed

    This makes donor×label the statistical unit rather than individual cells.
    """
    batches = (
        adata.obs[batch_key]
        .astype(str)
        .to_numpy()
    )
    labels = (
        adata.obs[label_key]
        .astype(str)
        .to_numpy()
    )

    gene_names = pd.Index(gene_names.astype(str))
    selected_genes = [str(g) for g in selected_genes]
    gene_index = gene_names.get_indexer(selected_genes)

    if np.any(gene_index < 0):
        missing = [
            selected_genes[i]
            for i, x in enumerate(gene_index)
            if x < 0
        ]
        raise KeyError(
            f"Missing genes in raw counts: {missing[:10]}"
        )

    full_totals = _to_numpy_1d(
        counts.sum(axis=1)
    ).astype(float)

    selected_counts = counts[:, gene_index]

    strata = {}
    for i, (b, l) in enumerate(zip(batches, labels)):
        strata.setdefault((b, l), []).append(i)

    rows = []
    batch_out = []
    label_out = []
    n_cells_out = []

    for (batch, label), idx_list in sorted(strata.items()):
        idx = np.asarray(idx_list, dtype=int)

        if sparse.issparse(selected_counts):
            summed = _to_numpy_1d(
                selected_counts[idx].sum(axis=0)
            ).astype(float)
        else:
            summed = np.asarray(
                selected_counts[idx].sum(axis=0)
            ).reshape(-1).astype(float)

        library_total = float(
            np.sum(full_totals[idx])
        )

        if not np.isfinite(library_total) or library_total <= 0:
            library_total = 1.0

        expr = np.log1p(
            summed / library_total * 1e6
        )

        rows.append(expr)
        batch_out.append(batch)
        label_out.append(label)
        n_cells_out.append(len(idx))

    return (
        np.vstack(rows),
        np.asarray(batch_out, dtype=str),
        np.asarray(label_out, dtype=str),
        np.asarray(n_cells_out, dtype=int),
    )

def find_harmful_v1_celllevel(
    adata,
    counts,
    gene_names,
    hvg_genes,
    batch_key,
    label_key,
    *,
    max_cells_per_batch_label=200,
    n_permutations=100,
    n_bootstraps=20,
    alpha=0.05,
    bootstrap_pass_fraction=0.80,
    leakage_z_floor=1.0,
    risk_z_floor=1.0,
    biology_z_ceiling=0.0,
    marker_log_effect=0.50,
    marker_replication_fraction=0.80,
    marker_min_eligible_donors=2,
    marker_min_cells_per_side=10,
    protected_genes=None,
    seed=20260829,
):
    """
    V1 = intentionally retained CELL-LEVEL baseline.

    Statistical unit:
        sampled individual cells.

    Core evidence:
        donor leakage after biology-label residualization
        + donor×biology max-min instability
        - biological signal

    Significance:
        conditional donor-label permutation within biology label.

    Stability:
        cell-level stratified bootstrap/subsampling recurrence.

    Biology protection:
        replicated marker protection across ELIGIBLE donors
        + explicit protected_genes.

    IMPORTANT:
      This is the baseline comparator. It still contains the weaknesses
      we want V2 to address: cell-level pseudoreplication risk and
      max-min donor instability.
    """
    protected_genes = set(
        map(str, protected_genes or [])
    )
    hvg_genes = list(map(str, hvg_genes))
    n_genes = len(hvg_genes)

    # --------------------------------------------------------
    # Stratified cell sampling
    # --------------------------------------------------------
    sample_idx = _stratified_sample_indices(
        adata.obs,
        batch_key,
        label_key,
        max_cells_per_batch_label,
        seed,
    )

    X = _normalize_selected_cells(
        counts,
        gene_names,
        hvg_genes,
        cell_indices=sample_idx,
    )

    obs = adata.obs.iloc[sample_idx].copy()
    batches = (
        obs[batch_key]
        .astype(str)
        .to_numpy()
    )
    labels = (
        obs[label_key]
        .astype(str)
        .to_numpy()
    )

    # --------------------------------------------------------
    # Observed metrics
    # --------------------------------------------------------
    biology = _biology_eta2(
        X,
        labels,
    )

    R = _residualize_by_label(
        X,
        labels,
    )

    leakage = _donor_leakage_from_label_residual(
        R,
        batches,
    )

    interaction = _celllevel_interaction_range(
        X,
        labels,
        batches,
    )

    leakage_z = robust_z(leakage)
    interaction_z = robust_z(interaction)
    biology_z = robust_z(biology)

    # Baseline-style score. Used for ranking AND V1 gate.
    raw_risk = (
        leakage_z
        + 0.75 * interaction_z
        - biology_z
    )
    risk_z = robust_z(raw_risk)

    # --------------------------------------------------------
    # Conditional permutation:
    # shuffle donor labels WITHIN biology labels
    # --------------------------------------------------------
    rng = np.random.default_rng(seed)
    perm_exceed = np.zeros(
        n_genes,
        dtype=int,
    )
    perm_max = np.zeros(
        n_permutations,
        dtype=float,
    )

    label_groups = _group_indices(labels)

    for p in range(n_permutations):
        perm_batch = batches.copy()

        for _, idx in label_groups.items():
            perm_batch[idx] = rng.permutation(
                perm_batch[idx]
            )

        perm_leakage = (
            _donor_leakage_from_label_residual(
                R,
                perm_batch,
            )
        )

        perm_exceed += (
            perm_leakage >= leakage
        )
        perm_max[p] = np.nanmax(
            perm_leakage
        )

    perm_p = (
        1 + perm_exceed
    ) / (
        1 + n_permutations
    )

    perm_fdr = bh_fdr(perm_p)

    maxT_p = np.asarray(
        [
            (
                1
                + np.sum(
                    perm_max >= obs_value
                )
            )
            / (1 + n_permutations)
            for obs_value in leakage
        ],
        dtype=float,
    )

    # --------------------------------------------------------
    # Cell-level bootstrap/subsampling recurrence
    # --------------------------------------------------------
    bootstrap_pass = np.zeros(
        n_genes,
        dtype=int,
    )

    strata = {}
    for i, (b, l) in enumerate(
        zip(batches, labels)
    ):
        strata.setdefault((b, l), []).append(i)

    for b in range(n_bootstraps):
        boot_idx = []

        for key in sorted(strata):
            idx = np.asarray(
                strata[key],
                dtype=int,
            )

            # Bootstrap within donor×biology strata.
            draw = rng.choice(
                idx,
                size=len(idx),
                replace=True,
            )
            boot_idx.extend(draw.tolist())

        boot_idx = np.asarray(
            boot_idx,
            dtype=int,
        )

        Xb = (
            X[boot_idx]
            if sparse.issparse(X)
            else np.asarray(X)[boot_idx]
        )

        labels_b = labels[boot_idx]
        batches_b = batches[boot_idx]

        biology_b = _biology_eta2(
            Xb,
            labels_b,
        )
        Rb = _residualize_by_label(
            Xb,
            labels_b,
        )
        leakage_b = (
            _donor_leakage_from_label_residual(
                Rb,
                batches_b,
            )
        )
        interaction_b = (
            _celllevel_interaction_range(
                Xb,
                labels_b,
                batches_b,
            )
        )

        leak_z_b = robust_z(
            leakage_b
        )
        int_z_b = robust_z(
            interaction_b
        )
        bio_z_b = robust_z(
            biology_b
        )

        risk_b = (
            leak_z_b
            + 0.75 * int_z_b
            - bio_z_b
        )
        risk_z_b = robust_z(
            risk_b
        )

        bootstrap_pass += (
            (leak_z_b >= leakage_z_floor)
            & (risk_z_b >= risk_z_floor)
            & (bio_z_b <= biology_z_ceiling)
        )

    bootstrap_fraction = (
        bootstrap_pass
        / max(1, n_bootstraps)
    )

    # --------------------------------------------------------
    # Marker protection is calculated from full data aggregates.
    # --------------------------------------------------------
    (
        group_expr,
        group_batches,
        group_labels,
        group_n_cells,
    ) = build_donor_label_aggregate(
        adata,
        counts,
        gene_names,
        hvg_genes,
        batch_key,
        label_key,
    )

    marker = _marker_protection_from_groups(
        group_expr,
        group_batches,
        group_labels,
        group_n_cells,
        hvg_genes,
        log_effect=marker_log_effect,
        replication_fraction=marker_replication_fraction,
        min_eligible_donors=marker_min_eligible_donors,
        min_cells_per_side=marker_min_cells_per_side,
    )

    marker_map = marker.set_index("gene")

    marker_protected = np.asarray(
        [
            bool(
                marker_map.loc[g, "marker_protected"]
            )
            for g in hvg_genes
        ],
        dtype=bool,
    )

    explicit_protection = np.asarray(
        [
            g in protected_genes
            for g in hvg_genes
        ],
        dtype=bool,
    )

    # --------------------------------------------------------
    # V1 final rule
    # --------------------------------------------------------
    significant = (
        perm_fdr <= alpha
    )

    risk_flagged = (
        significant
        & (leakage_z >= leakage_z_floor)
        & (risk_z >= risk_z_floor)
        & (biology_z <= biology_z_ceiling)
        & (
            bootstrap_fraction
            >= bootstrap_pass_fraction
        )
    )

    harmful = (
        risk_flagged
        & (~marker_protected)
        & (~explicit_protection)
    )

    table = pd.DataFrame(
        {
            "gene": hvg_genes,
            "donor_leakage": leakage,
            "donor_leakage_z": leakage_z,
            "interaction_range": interaction,
            "interaction_range_z": interaction_z,
            "biology_eta2": biology,
            "biology_z": biology_z,
            "risk_score_raw": raw_risk,
            "risk_score_z": risk_z,
            "permutation_p": perm_p,
            "permutation_fdr": perm_fdr,
            "permutation_maxT_p": maxT_p,
            "bootstrap_risk_fraction": bootstrap_fraction,
            "risk_flagged_before_protection": risk_flagged,
            "marker_protected": marker_protected,
            "explicitly_protected": explicit_protection,
            "harmful": harmful,
        }
    )

    table = table.merge(
        marker,
        on=[
            "gene",
            "marker_protected",
        ],
        how="left",
    )

    table["decision_reason"] = np.select(
        [
            table["explicitly_protected"],
            table["marker_protected"],
            table["harmful"],
            table["risk_flagged_before_protection"],
        ],
        [
            "explicit_protection",
            "replicated_marker_protection",
            "harmful_celllevel_v1",
            "risk_but_not_removed",
        ],
        default="retain",
    )

    table = table.sort_values(
        [
            "harmful",
            "risk_score_z",
        ],
        ascending=[
            False,
            False,
        ],
        kind="mergesort",
    ).reset_index(drop=True)

    harmful_genes = (
        table.loc[
            table["harmful"],
            "gene",
        ]
        .astype(str)
        .tolist()
    )

    return {
        "method": "V1_celllevel_baseline",
        "harmful_genes": harmful_genes,
        "table": table,
        "sampled_cells": int(len(sample_idx)),
        "n_hvg": int(len(hvg_genes)),
    }
