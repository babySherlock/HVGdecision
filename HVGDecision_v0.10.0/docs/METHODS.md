# Methods and default parameters

The current entry point is `hvgdecision.refine`. It does not run pseudobulk
differential-expression discovery, fit an integration model, or optimize a
downstream score. It uses the full validated count axis for library normalization
and assesses genes only in the supplied or independently fitted base panel.

## Selection cohort and base panel

All input cells form the selection cohort unless `reference` donor IDs are
provided. Metadata must be nonmissing and cell/gene IDs unique. The built-in
selector calls `scanpy.pp.highly_variable_genes(flavor='seurat_v3',
batch_key=batch_key, n_top_genes=n_hvg)` on raw selection-cohort counts.
The default base size is 2000. Each requested size is an independent fit.
Selected genes are ordered by `highly_variable_rank` ascending, then batch count
descending, then gene ID to resolve remaining ties. External panel order is
preserved and its length overrides `n_hvg`; no automatic rank reinterpretation
is performed. The algorithm never replenishes removed genes.

## Replication adequacy

Let n(d,k) be selection-cohort cell counts for donor d and cell type k, with D
donors and K types. The implementation uses:

```text
s(d,k) = n(d,k) / (n(d,k) + tau)
r(k)   = sum_d s(d,k)
v(k)   = min(1, sum_d n(d,k) / (D*tau))
D_eff  = sum_k v(k) / sum_k [v(k) / (r(k) + 1e-12)]
u(d)   = mean_k s(d,k)
B_D    = [sum_d u(d)]^2 / [D * sum_d u(d)^2]
A      = D_eff * sqrt(B_D)
w      = clip((A - 2) / (4 - 2), 0, 1)
```

Default `tau=15`. The harmonic mean downweights poorly supported types; v(k)
reduces the influence of types represented by very few cells. D is the number
of observed selection donors, not the number originally collected.

## Identifiability and routing

An edge connects donor d and type k when n(d,k) >= tau. Remove isolated vertices
and let D_s, K_s be the numbers of supported donors and types, and C the number
of connected components. The rank of the intercept + donor/type dummy design
on supported strata is D_s + K_s - C. The expected identifiable rank is
D_s + K_s - 1. Their ratio is iota. The code computes this exact graph identity
instead of allocating a dense design matrix.

If D_s < 2, K_s < 2, or C != 1, the algorithm returns
`insufficient_confounded` and makes no removal calls. Otherwise `w >= 0.75`
selects `donor_aware`; lower w selects `cell_level`. w does not blend the two
risk scores. It is a fixed design heuristic, not a calibrated probability.
An explicit branch override is recorded and cannot bypass this safeguard.

Identifiability is checked on supported strata, not a guarantee that every
donor/type pair is adequately sampled. The output includes unsupported pairs
and the fraction of cells belonging to supported strata.

## Cell-level branch: Generic V1

Sample at most 200 cells without replacement in each donor-by-type stratum.
Normalize selected genes using each cell's full count-library total (target
10,000) and log1p. Define B_g as the fraction of expression variance explained
by cell type. Center expression within each type; L_g is the fraction of this
residual sum of squares explained by donor means. For each type represented
by >=2 donors, compute donor-mean range divided by within-type standard
deviation (+1e-6); I_g is its median across eligible types.

```text
R_g = Z( Z(L_g) + 0.75*Z(I_g) - Z(B_g) )
```

Z uses the across-panel median and 1.4826*MAD. A near-zero MAD falls back to
standard deviation, then to 1. The permutation statistic is L_g, NOT R_g.
Permute donor identities within cell types 100 times, using
`p_g=(1 + count(L_perm >= L_observed))/(1 + N_perm)` and BH correction across
the entire input panel. maxT-adjusted p values are reported only as diagnostics.

Run 20 bootstrap repetitions by sampling with replacement, at the original
sampled stratum size, within donor-by-type strata. Recompute the robust scores
for each resample. A bootstrap pass requires Z(L)>=1, R>=1 and Z(B)<=0.
The final pre-protection risk call requires all of:

| Parameter | Default criterion |
| --- | --- |
| permutation FDR | <= 0.05 |
| donor leakage Z | >= 1.0 |
| risk score Z | >= 1.0 |
| biology Z | <= 0.0 |
| bootstrap pass fraction | >= 0.80 |

Replicated marker protection uses donor-by-type sums of counts from the full
selection cohort, with full-gene library denominators (log1p CPM, target 1e6).
For each target type, an eligible donor has >=10 target cells and >=10 other
cells. A gene is protected when target-minus-other log-expression is >=0.50
in >=80% of eligible donors, with at least two eligible donors. This aggregation
is used for protection only; it is not the risk-discovery statistic or the
Figure 1 illustrative pseudobulk analysis. Explicit protected genes are also
retained. The cell-level engine has no separate rare-marker gate.

This is cell-level conditional evidence, not donor-level hypothesis testing.
Permutation resolution is 1/101 by default; BH correction and bootstrap
stability can legitimately yield zero removals. Do not replace BH FDR with the
maxT diagnostic simply because the latter permits a desired gene.

## Donor-aware branch: donor-replicate v2.0

Log-normalize counts to 10,000 using full libraries. Compute B_g (cell-type eta2)
and N_g (donor-by-type mean deviations from the type mean, weighted by stratum
cell counts and divided by total variance). Strata require >=15 cells.

For each gene choose the type with highest pooled mean. In each donor with
adequate target and other-type cells, calculate target-minus-other expression
effect e(d,g). Let A_g be the fraction of finite effects >0, D_g the largest
absolute donor effect divided by the sum of absolute effects, and H_g the
standard deviation of donor effects divided by their mean absolute value.
J_g is the corresponding instability of leave-one-donor-out mean effects.
H_g and J_g are clipped to [0,5]. Let P denote average-tie percentile rank
across the input panel (nonfinite ranks map to zero):

```text
T_g = 0.35*P(N_g) + 0.25*P(D_g) + 0.20*P(H_g) + 0.20*P(J_g)
S_g = 0.55*P(B_g) + 0.45*A_g
Q_g = T_g - S_g
```

Here A_g is effect direction agreement, distinct from design-level adequacy A.
It is not a bootstrap frequency. No permutation FDR is computed in this branch.
All these risk conditions must hold:

| Condition | Default |
| --- | --- |
| Technical evidence T_g | >= panel 97.5th percentile |
| Risk Q_g | >= max(0.30, panel 98.5th percentile) |
| Biological support S_g | <= panel 45th percentile |
| Finite donor effects | >=3 |
| Dominance / LODO evidence | D_g >=0.55 OR LODO positive fraction <0.80 |
| Direction agreement A_g | <0.75 |

Replicated markers require B_g >= its panel 90th percentile, direction >=0.80,
at least max(2, ceil(D/2)) finite donor effects and LODO positive fraction >=0.80.
Rare-marker protection requires a highest-mean type with pooled prevalence <=5%,
B_g >= its panel 75th percentile, direction >=0.75 and the same minimum donor
replication. These gene-protection gates are distinct from a downstream rare-cell
F1 definition. The package does not calculate rare-cell F1 or implement the
evaluation criterion of 50 cells and 25% donor coverage here.

All replicated markers, rare markers and explicit protected genes are retained.
There is no fixed removal count and no built-in real gene name list. Remaining
hard-coded constants belong to the frozen engine and are documented above;
they are not data-trained parameters.

## Interpretation and reproducibility

The two branches use different statistics and protections. They share the
panel-removal interface and identifiability safeguard; their numeric risk scores
are not interchangeable. Do not report donor-aware rows with fictitious FDR or
bootstrap values. Software schema columns unavailable in a branch remain absent
or missing, never zero-filled to imply testing.

Risk scores are relative to the tested gene panel. Identical counts, cell order,
gene order, selection cohort, annotations, source implementation, parameters
and random seed are needed for reproduction. The release freezes formulas but
does not assert reproduction of historical PBMC or Lung deletion counts without
the original inputs. Engine hashes and version records are saved for every run.

The cell-level path can use substantial dense memory after stratified sampling.
Use the metadata audit before starting a large run; reducing the sampling cap
changes the analysis and is recorded in the manifest.
