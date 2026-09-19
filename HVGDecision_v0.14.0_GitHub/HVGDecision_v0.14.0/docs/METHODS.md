# HVGDecision v0.14.0 — mathematical specification

## Overview

HVGDecision refines an existing HVG panel before downstream integration. For every candidate gene, it computes two complementary evidence components and combines them on a common standardized signed-margin scale.

The two evidence components are:

- **Cell-level conditional donor-leakage evidence**
- **Donor-aware evidence**

Replication adequacy and design identifiability are reported separately as dataset-level audits and do not enter the gene-level fusion equation.

## 1. Replication-adequacy audit

Let `n_dk` be the number of cells from donor `d` assigned to cell type `k`.

Donor–cell-type support is

\[
s_{dk}=\frac{n_{dk}}{n_{dk}+\tau_s},\qquad \tau_s=15.
\]

Effective support for cell type `k` is

\[
r_k=\sum_d s_{dk}.
\]

The cell-type audit weight is

\[
v_k=\min\left(1,\frac{\sum_d n_{dk}}{D\tau_s}\right).
\]

The weighted harmonic effective donor support is

\[
D_{eff}=\frac{\sum_k v_k}{\sum_k v_k/(r_k+\epsilon)}.
\]

For donor `d`,

\[
u_d=\frac{1}{K}\sum_k s_{dk},
\]

and donor balance is

\[
B_D=\frac{(\sum_d u_d)^2}{D\sum_d u_d^2+\epsilon}.
\]

Replication adequacy is

\[
A=D_{eff}\sqrt{B_D}.
\]

`A` is an **audit value only** in v0.14.0.

## 2. Design-identifiability audit

A donor–cell-type combination is treated as supported when

\[
e_{dk}=\mathbf{1}(n_{dk}\ge\tau_s).
\]

The supported bipartite design is identifiable when the corresponding intercept + donor + cell-type design has full rank, equivalently one connected supported donor–cell-type component with at least two supported donors and two supported cell types.

The reported identifiability ratio is

\[
\iota=\frac{\operatorname{rank}(X_{design})}{D_s+K_s-1}.
\]

## 3. Cell-level conditional donor-leakage evidence

The cell-level component is evaluated from individual-cell expression after conditioning on cell type. It includes:

- conditional donor leakage;
- within-cell-type donor instability;
- cell-type biological support;
- within-cell-type donor-label permutation FDR;
- within-stratum bootstrap stability.

The component reproduces the manuscript criteria and produces criterion-specific standardized signed margins.

## 4. Donor-aware evidence

The donor-aware component includes:

- within-cell-type donor-associated nuisance;
- single-donor effect dominance;
- cross-donor heterogeneity of the cell-type effect;
- leave-one-donor-out instability;
- cell-type-associated biological support;
- directional agreement across donors;
- effective donor count.

The technical non-replication score is

\[
T_g = 0.35\,PR(N_g)+0.25\,PR(U_g)+0.20\,PR(H_g)+0.20\,PR(J_g).
\]

Cross-donor biological support is

\[
V_g=0.55\,PR(C_g^{(d)})+0.45a_g,
\]

and donor-aware net risk is

\[
R_g^{(d)}=T_g-V_g.
\]

## 5. Standardized signed margins

For criterion value `x_gj`, threshold `t_j` and robust scale `s_j`:

\[
m_{gj}=\frac{x_{gj}-t_j}{s_j}
\]

when higher values favor risk, and

\[
m_{gj}=\frac{t_j-x_{gj}}{s_j}
\]

when lower values favor risk.

Positive margins satisfy the corresponding risk condition; negative margins fail it. Exact threshold equalities are sign-preserved with a small numerical epsilon.

Conjunctive criteria use the minimum margin. Alternative criteria use the maximum margin.

Thus,

\[
M_g^{(c)}=\min_j m_{gj}^{(c)},
\]

and the donor-aware component has the structure

\[
M_g^{(d)}=
\min\left\{
 m_{T,g},
 m_{R_d,g},
 m_{V,g},
 m_{m,g},
 \max(m_{U,g},m_{LODO,g}),
 m_{a,g}
\right\}.
\]

## 6. Component-specific biological protection

For component `b` in `{c,d}`, let `pi_g^(b)` be the corresponding protection indicator.

The effective component margin is

\[
\widetilde M_g^{(b)}=
\begin{cases}
M_g^{(b)}, & \pi_g^{(b)}=0,\\
-\infty, & \pi_g^{(b)}=1.
\end{cases}
\]

Protection disables harmful evidence from the corresponding component without imposing a global veto on the other component.

## 7. Dual-evidence centered log-mean-exp fusion

The final unified signed margin is

\[
M_g^{(u)}
=\tau_f\log\left[
\frac{
\exp(\widetilde M_g^{(d)}/\tau_f)
+\exp(\widetilde M_g^{(c)}/\tau_f)
}{2}
\right].
\]

Default:

\[
\tau_f=0.10.
\]

The final discovery-risk indicator is

\[
h_g=\mathbf{1}(M_g^{(u)}>0).
\]

### Mathematical properties

If the two effective component margins are equal to `m`, then

\[
M_g^{(u)}=m.
\]

When one component strongly dominates the other,

\[
M_g^{(u)}\approx
\max(\widetilde M_g^{(d)},\widetilde M_g^{(c)})
-\tau_f\log2.
\]

As `tau_f -> 0+`, the unified margin approaches the larger effective component margin. Increasing `tau_f` increases the disagreement penalty.

There is no dataset-dependent donor-versus-cell fusion weight.

## 8. Reference → Query transport

For a Reference discovery panel,

\[
R_{ref}=\{g\in G_{ref}:M_g^{(u)}>0\}.
\]

If a Query HVG panel is supplied, the final removed set is

\[
G_{removed}=G_{query}\cap R_{ref},
\]

and

\[
G_{final}=G_{query}\setminus G_{removed}.
\]

Thus discovery-risk counts and final removal counts can differ.

## 9. Pooled analysis

For pooled analysis, the discovery panel and final base panel are the same pooled HVG panel, so every discovery-risk gene in that panel is removed.

## 10. Default parameters

Current release defaults relevant to the final fusion:

```text
fusion temperature tau_f = 0.10
discovery cutoff = 0
support half-saturation tau_s = 15
boundary epsilon = 1e-9
scale floor = 1e-8
```

The remaining donor-aware and cell-level gate parameters are exposed through `DonorAwareConfig`, `CellLevelConfig`, `DonorGateConfig` and `CellGateConfig` and are recorded in the run manifest.
