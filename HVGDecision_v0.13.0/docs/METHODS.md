# Continuous evidence fusion

## Replication adequacy

Let \(n_{dk}\) be the discovery-cell count for donor \(d\) and cell type \(k\), with \(D\) observed donors and \(K\) observed types. With \(\tau=15\),

\[
s_{dk}=\frac{n_{dk}}{n_{dk}+\tau},\quad r_k=\sum_d s_{dk},
\quad v_k=\min\left(1,\frac{\sum_d n_{dk}}{D\tau}\right).
\]

The effective support is a weighted harmonic mean:

\[
D_{\mathrm{eff}}=\frac{\sum_k v_k}{\sum_k v_k/(r_k+10^{-12})}.
\]

Writing \(u_d=K^{-1}\sum_k s_{dk}\),

\[
B_{\mathrm{bal}}=\frac{(\sum_d u_d)^2}{D\sum_d u_d^2},\qquad
A=D_{\mathrm{eff}}\sqrt{B_{\mathrm{bal}}}.
\]

Use \(B_{\mathrm{bal}}\) for donor balance to distinguish it from gene-level biological evidence. All quantities use the discovery cohort only.

## Fusion and decisions

Both evidence components are calculated. Their donor-component weight is

\[
w(A)=\frac{1}{1+\exp[-(A-3)/0.5]}.
\]

Supported-design identifiability is audited separately; a high A does not establish identifiability.

Each component retains its evidence statistics and thresholds. Gate distances are oriented toward passing and divided by a robust scale (MAD, with standard-deviation and unit fallbacks). At exact boundaries, signs follow the original comparison operators, using epsilon \(10^{-9}\). Nonfinite distances fail.

Conjunctions use the minimum signed margin. The donor dominance-or-LODO condition uses the maximum of those two margins within the conjunction. The sigmoid of the resulting margin, with temperature 1, gives component confidence. Component-specific protected genes receive confidence zero.

\[
U_g=w(A)H_{g,\mathrm{donor}}+[1-w(A)]H_{g,\mathrm{cell}}.
\]

The final discovery call is strictly \(U_g>0.5\). This is not a posterior probability, and no FDR control is asserted for the fused set. Explicit user protection is applied to both components. For Reference/Query analysis, intersect discovery calls with the independently selected Query panel.

All executable thresholds and operators are in parameters.py, fusion.py and _continuous_core.py. Saved runs include their configurations and component evidence.

