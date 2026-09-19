# HVGDecision

HVGDecision 用于多 donor 单细胞 RNA-seq 数据整合前的 HVG 面板精炼。当前 **v0.14.0** 对每个候选基因同时评估两类证据：

1. **Cell-level conditional donor-leakage evidence**：控制 cell type 后，从单细胞表达中检测 donor-associated variation；
2. **Donor-aware evidence**：评估跨 donor 的非重复性、single-donor dominance、heterogeneity、LODO instability、direction agreement 和 biological support。

两个 evidence component 都会被转成 standardized signed margin，并在各自 biological protection 后，通过 centered log-mean-exp 融合：

\[
M_g^{(u)} = \tau_f\log\left[\frac{\exp(\widetilde M_g^{(d)}/\tau_f)+\exp(\widetilde M_g^{(c)}/\tau_f)}{2}\right],
\qquad \tau_f=0.10.
\]

最终 discovery-risk 判定：

\[
M_g^{(u)}>0.
\]

**Replication adequacy 和 design identifiability 仅作为 dataset-level audit，不进入 gene-level fusion。**

[英文 README](README.md) · [数学说明](docs/METHODS.md) · [v0.13 → v0.14 迁移](docs/MIGRATION.md) · [Counts](docs/COUNTS.md)

## 安装

推荐 Python 3.10 及以上，并安装到你实际运行 notebook 的同一个环境。

直接安装 wheel：

```bash
python -m pip install ./dist/hvgdecision-0.14.0-py3-none-any.whl
```

或者从源码目录安装：

```bash
python -m pip install .
```

安装后：

```python
import hvgdecision as hd
print(hd.__version__)
```

应输出：

```text
0.14.0
```

## Reference / Query 用法

```python
from datetime import datetime
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad('/path/to/filtered_counts.h5ad')

counts = hd.find_raw_counts(adata)
if not counts.valid:
    raise ValueError(counts.audit.to_string(index=False))

reference = [
    'patient_101', 'patient_1015', 'patient_1016',
    'patient_1039', 'patient_107', 'patient_1244',
]
query = ['patient_1256', 'patient_1488']

result = hd.refine(
    adata,
    batch_key='replicate',
    label_key='cell_type',
    counts=counts,
    reference=reference,
    query=query,
    n_hvg=2000,
    hvg_span=0.5,
    query_hvg_span=0.3,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    output_dir='./HVGDecision_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'),
    return_details=True,
)

print(result)
print('Reference discovery-risk genes:', result.reference_risk_genes)
print('Final removed genes:', result.removed_genes)
print(result.design.audit.to_string(index=False))

adata_hvg = result.adata
```

Reference 和 Query 必须互斥，并覆盖当前输入对象中的全部 donor。Reference 用于生成风险证据；Query HVG 独立计算，但 Query cell-type label 不参与 feature refinement。

最终删除集合为：

\[
R_{final}=R_{Reference\ discovery}\cap HVG_{Query}.
\]

因此，**Reference discovery-risk gene 数量大于最终 Query 删除数量是正常设计结果**。

## Pooled 用法

如果是同一组织、同一分析场景下的多 donor pooled 数据，不传 `reference` / `query`：

```python
result = hd.refine(
    adata,
    batch_key='donor_id',
    label_key='cell_type',
    counts=hd.find_raw_counts(adata),
    n_hvg=2000,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    return_details=True,
)
```

此时 discovery panel 和最终 base panel 是同一个 pooled HVG panel。

## 使用已经冻结的 HVG panel

如果你已经有 Reference / Query HVG2000：

```python
import pandas as pd

reference_genes = pd.read_csv('reference_hvg_panel.csv')['gene'].astype(str).tolist()
query_genes = pd.read_csv('query_hvg_panel.csv')['gene'].astype(str).tolist()

result = hd.refine(
    adata,
    batch_key='replicate',
    label_key='cell_type',
    counts=counts,
    reference=reference,
    query=query,
    hvg_genes=reference_genes,
    query_hvg_genes=query_genes,
    fusion_config=hd.UnifiedFusionConfig(tau=0.10),
    return_details=True,
)
```

外部 panel 会跳过 Seurat v3 HVG fitting。不存在于 counts gene axis 的基因会直接报错，不会静默删除或自动改名。

## 当前数学逻辑

### 1. Cell-level conditional donor-leakage evidence

它不是“只看细胞、不看 donor”。它是在控制 cell type 后，从 individual-cell expression 中检查：

- donor leakage；
- within-cell-type donor instability；
- cell-type biological support；
- permutation FDR；
- bootstrap stability。

### 2. Donor-aware evidence

主要检查：

- technical non-replication；
- single-donor dominance；
- cross-donor heterogeneity；
- leave-one-donor-out instability；
- direction agreement；
- biological support。

### 3. Signed margin

每一个 gate 都被转成 signed standardized margin：

- 正值：满足 risk criterion；
- 负值：不满足 risk criterion。

AND 条件取 minimum，OR 条件取 maximum，得到：

\[
M_g^{(c)},\qquad M_g^{(d)}.
\]

### 4. Biological protection

\[
\widetilde M_g^{(b)}=
\begin{cases}
M_g^{(b)}, & \pi_g^{(b)}=0,\\
-\infty, & \pi_g^{(b)}=1.
\end{cases}
\]

被某个 component 保护后，该 component 不再提供 harmful evidence，但不会全局否决另一个 component。

### 5. Centered log-mean-exp fusion

\[
M_g^{(u)} = \tau_f\log\left[\frac{\exp(\widetilde M_g^{(d)}/\tau_f)+\exp(\widetilde M_g^{(c)}/\tau_f)}{2}\right].
\]

默认：

```python
hd.UnifiedFusionConfig(tau=0.10)
```

性质：

- 两个 component 数学上对称；
- 如果两个 margin 都等于 `m`，则 unified margin 也严格等于 `m`；
- 一个 component 单独占优时会承担约 `tau_f * log(2)` 的 disagreement penalty；
- `tau_f` 越小，越接近较大的 component margin；
- `tau_f` 越大，对两个 component 不一致的惩罚越强；
- 不存在 dataset-dependent donor/cell weight。

## Replication adequacy 现在做什么

仍然会计算：

\[
A=D_{eff}\sqrt{B_D}
\]

但它只用于描述数据集 donor × cell-type replication structure，不进入最终 gene-level fusion。

```python
report = hd.audit_design(
    adata,
    batch_key='donor',
    label_key='cell_type',
)
print(report.audit)
```

## 输出

`return_details=True` 时，`result` 主要包含：

- `result.adata`：最终 retained HVG 的 raw-count AnnData；
- `result.removed_genes`：最终 base/query panel 中被删除的基因；
- `result.reference_risk_genes`：Reference 或 pooled discovery-risk genes；
- `result.decision_table`：最终 target panel 决策；
- `result.reference_decision_table`：discovery panel 的完整 gene-level evidence 和决策；
- `result.component_evidence['donor']`；
- `result.component_evidence['cell']`；
- `result.design`：replication adequacy + design identifiability audit；
- `result.run_info`：参数、版本、hash 和 fusion manifest。

如果指定 `output_dir`，会输出：

```text
adata_hvg.h5ad
gene_decisions.csv
reference_gene_decisions.csv
reference_risk_genes.csv
removed_genes.csv
retained_genes.csv
donor_evidence.csv
cell_evidence.csv
design_audit.csv
donor_celltype_coverage.csv
design_by_donor.csv
design_by_celltype.csv
run_manifest.json
...
```

## 命令行

仅审计设计结构：

```bash
hvgdecision audit \
  --input data.h5ad \
  --batch-key donor \
  --label-key cell_type \
  --output audit_out
```

运行 refinement：

```bash
hvgdecision refine \
  --input data.h5ad \
  --batch-key donor \
  --label-key cell_type \
  --n-hvg 2000 \
  --tau-f 0.10 \
  --output refinement_out
```

## v0.13 和 v0.14 最大区别

v0.13 的最终层是：

```text
signed margin
→ sigmoid confidence
→ replication-adequacy-dependent weighted mean
→ threshold 0.5
```

v0.14 改成：

```text
signed margin
→ component-specific protection
→ centered log-mean-exp
→ threshold 0
```

也就是说，**前面的 donor-aware / cell-level evidence engine 基本保留，真正重写的是最终 fusion 数学逻辑。**

详细说明见 [docs/MIGRATION.md](docs/MIGRATION.md)。

## 测试

```bash
python -m pip install '.[dev]'
python -m pytest -q
```

## License

MIT License。仅用于科研软件，不用于临床解释。
