# Migration: HVGDecision v0.13 → v0.14

## What changed

Version 0.13 used the following final decision layer:

```text
component signed margins
→ sigmoid component confidences
→ replication-adequacy-dependent weighted mean
→ threshold at 0.5
```

Formally, the v0.13 fused score was of the form

\[
F_g=w(A)\psi_g^{(d)}+[1-w(A)]\psi_g^{(c)}.
\]

Version 0.14 replaces that final layer with direct fusion on the signed-margin scale:

```text
component signed margins
→ component-specific protection
→ centered log-mean-exp
→ threshold at 0
```

\[
M_g^{(u)}=\tau_f\log\left[
\frac{
\exp(\widetilde M_g^{(d)}/\tau_f)
+\exp(\widetilde M_g^{(c)}/\tau_f)
}{2}
\right].
\]

The default is `tau_f = 0.10` and the discovery call is `M_u > 0`.

## What did not change conceptually

The two evidence engines remain the basis of the method:

- cell-level conditional donor-leakage evidence;
- donor-aware cross-donor evidence.

The signed gate reconstruction also remains: individual risk criteria are converted to standardized signed margins, AND conditions use a minimum and OR conditions use a maximum.

## Replication adequacy

Replication adequacy is still computed because it is useful for characterizing donor × cell-type support, but it no longer changes the gene-level fusion result.

In v0.14:

```text
replication adequacy = audit metadata only
```

## API migration

### Old

```python
cfg = hd.ContinuousFusionConfig(...)
```

### New

```python
cfg = hd.UnifiedFusionConfig(tau=0.10)
```

Then:

```python
result = hd.refine(
    adata,
    batch_key='donor',
    label_key='cell_type',
    fusion_config=cfg,
    return_details=True,
)
```

The main high-level call remains `hd.refine(...)`.

## Result object

The primary design audit is now:

```python
result.design
```

For transition from v0.13, `result.routing` remains as a deprecated compatibility alias and emits a `FutureWarning`.

## Important output-column changes

Current decision tables use:

```text
donor_branch_margin_raw
cell_branch_margin_raw
donor_branch_margin_effective
cell_branch_margin_effective
unified_margin
fusion_tau
discovery_cutoff
discovery_remove
final_action
```

The v0.13 `unified_harmfulness`, donor/cell confidence weighting, and dataset-dependent fusion-weight fields are not part of the v0.14 public decision rule.

## Reproducibility

When reproducing a manuscript benchmark generated with a frozen HVG panel, pass the same ordered Reference / pooled and Query HVG lists to v0.14. Matching final panel sizes is not sufficient; exact gene-set identity should be verified.
