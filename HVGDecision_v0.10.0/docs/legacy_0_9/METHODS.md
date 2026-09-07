# HVGDecision 0.9.1 — formal method definition

## Explicit mode selection

Every run must select one of two modes.

```text
within_domain
cross_domain
```

The modes share raw-count validation, Reference/Query bookkeeping, marker protection and auditable final-panel output, but use different risk definitions.

---

## Mode 1: within_domain

### Biological explained variance

For gene \(g\), biological support is the fraction of total expression variance explained by the declared biology strata:

$$
B_g=
\frac{\sum_c n_c(\bar{x}_{cg}-\bar{x}_g)^2}
{\sum_i(x_{ig}-\bar{x}_g)^2+\epsilon}.
$$

### Cell-type/biology residualization and donor leakage

$$
r_{ig}=x_{ig}-\bar{x}_{c_i g}
$$

$$
L_g=
\frac{\sum_d n_d\bar r_{dg}^{2}}
{\sum_i r_{ig}^{2}+\epsilon}.
$$

### Donor × biology interaction instability

For each eligible biology stratum \(c\):

$$
I_{gc}=\frac{\max_d\bar{x}_{cdg}-\min_d\bar{x}_{cdg}}
{s_{cg}+\epsilon},
$$

and

$$
I_g=\operatorname{median}_c(I_{gc}).
$$

### Robust normalization and raw risk

Each component is robustly standardized across genes using median/MAD, with an SD fallback when MAD is degenerate.

$$
\boxed{R_g^{within}=Z(L_g)+0.75Z(I_g)-Z(B_g)}
$$

A second robust standardization gives \(Z_R(g)=Z(R_g^{within})\).

### Conditional permutation

Donor labels are permuted within biological strata. Per-gene leakage significance is

$$
p_g=\frac{1+\sum_{p=1}^{P}\mathbf{1}(L_g^{(p)}\ge L_g)}{P+1},
$$

followed by Benjamini-Hochberg FDR correction.

Defaults:

```text
P = 100
alpha = 0.05
```

### Bootstrap recurrence

Cells are resampled within donor × biology strata. A bootstrap replicate passes when

$$
Z(L_g^{(b)})\ge1,
\quad Z_R^{(b)}(g)\ge1,
\quad Z(B_g^{(b)})\le0.
$$

The recurrence fraction is

$$
F_g^{boot}=\frac1{N_{boot}}\sum_b\mathbf1[\mathrm{pass}_b].
$$

Defaults:

```text
N_boot = 20
minimum recurrence = 0.80
```

### Final within-domain harmful gate

$$
\boxed{
H_g^{within}=\mathbf1[
q_g\le0.05\land
Z(L_g)\ge1\land
Z_R(g)\ge1\land
Z(B_g)\le0\land
F_g^{boot}\ge0.80\land
\neg P_g]
}
$$

where \(P_g\) is replicated-marker or explicit biological protection.

The final panel is

$$
\mathcal H_{final}=\mathcal H_{base}\setminus\{g:H_g^{within}=1\}.
$$

Zero deletion is valid.

### Within-domain HVG-budget recommendation

The existing Reference-only v0.8 minimum-sufficient budget search is retained for `study.run()` / `find_best_hvg()` in within-domain mode. It uses donor-held-out biological sufficiency, rare-cell support, within-biology donor-instability guardrails and HVG-stability quality floors. Query expression and Query labels are excluded from budget selection.

---

## Mode 2: cross_domain

Cross-domain mode refines a **Query HVG panel**. A same-size Reference HVG panel is fitted independently and scored by the Within-domain V1 evidence model above.

### Reference Rule percentile

The Reference raw risk is converted to a percentile:

$$
R_g\in[0,1].
$$

### Label-free Reference→Query shift

Expression is library-size normalized and log transformed. The standardized mean shift is

$$
E_g=
\frac{|\mu_{Q,g}-\mu_{R,g}|}
{\sqrt{(Var_{R,g}+Var_{Q,g})/2+0.05^2}}.
$$

Detection-rate shift is

$$
D_g=|\pi_{Q,g}-\pi_{R,g}|.
$$

After percentile ranking,

$$
\boxed{S_g=0.70E_g^*+0.30D_g^*}.
$$

Query true labels are not used.

### Reference biological protection

$$
\boxed{
B_g=0.60B_g^{celltype}+0.25M_g+0.15P_g
}
$$

where \(B_g^{celltype}\) is the percentile of Reference biology eta-squared, \(M_g\) is marker replication fraction and \(P_g\) is a hard protection flag. Hard-protected genes are excluded from deletion.

### Dual technical consensus

$$
T_g=\min(R_g,S_g)
$$

$$
A_g=1-|R_g-S_g|.
$$

### Final Cross-domain Rule V3

$$
\boxed{
R_g^{cross}
=
T_g(0.75+0.25A_g)(1-0.75B_g)
}
$$

or equivalently

$$
\boxed{
R_g^{cross}
=
\min(R_g,S_g)
[0.75+0.25(1-|R_g-S_g|)]
(1-0.75B_g)
}.
$$

The user explicitly selects a deletion budget \(k\). Primary audited budgets are 5, 10 and 20; the software accepts any non-negative integer not exceeding the eligible candidate count.

$$
\mathcal H_{final}^{Q}=\mathcal H_{base}^{Q}\setminus\operatorname{TopK}(R_g^{cross}).
$$

---

## Marker protection

For every declared `biology_key`, donor-specific biological enrichment is evaluated. A marker can receive hard protection when the same enriched group is reproduced across a sufficient fraction of eligible Reference donors. Defaults are:

```text
marker_fdr = 0.01
marker_min_log_effect = 0.50
marker_replication_fraction = 0.80
marker_min_eligible_donors = 2
marker_min_cells_per_side = 10
```

`protected_genes` adds project-specific explicit hard protection.

---

## Query leakage discipline

### within_domain

Feature-selection risk uses Reference expression, Reference batch IDs and Reference biological labels. Query expression and Query labels are excluded from the risk decision.

### cross_domain

Feature selection uses:

- Reference counts, batch IDs and biological labels;
- Query counts and Query batch IDs for Query HVG selection / shift estimation.

Query **true biological labels are never used** by the Cross-domain Rule.

---

## Optional three-domain benchmark aggregation

This is an evaluation helper and is not part of feature selection.

Within each fixed `dataset × integration_method` task, each metric is min-max normalized across feature panels.

Transfer:

$$
S_T=mean(\widetilde{MacroF1},\widetilde{BA},\widetilde{RareF1}).
$$

Batch correction:

$$
S_B=mean(\widetilde{DatasetMixing},\widetilde{DonorMixingWithinCellType}).
$$

Biological fidelity:

$$
S_F=mean(\widetilde{ARI},\widetilde{NMI},\widetilde{CellTypeSilhouette}).
$$

Overall:

$$
\boxed{S_{overall}=0.30S_T+0.40S_B+0.30S_F}.
$$

This should be described as a **scIB-inspired three-domain composite score**, not an exact scIB score.
