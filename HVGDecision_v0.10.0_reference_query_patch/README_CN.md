# HVGDecision v0.10.0 Reference→Query 补丁（core-preserving）

这个补丁**不重写 HVGDecision v0.10.0 的 donor-aware 风险算法**。Reference 风险学习始终调用当前已安装的：

```python
import hvgdecision as hd
hd.refine(...)
```

补丁只补上两个当前教程/工作流缺失的部分：

1. 每一次 Seurat-v3 HVG 选择独立执行 `span=0.3 → 0.5 → 0.7 → 1.0` fallback；
2. Reference 风险名单转移到**独立选择的 Query HVG2000**，最终只删除交集。

因此 Reference 如果在 `span=0.3` 失败、`0.5` 成功，不会使 Query 直接从 `0.5` 开始；Query 会重新从 `0.3` 开始。

## 安装

先保留并安装你的 HVGDecision v0.10.0，然后安装本补丁：

```bash
python -m pip install ./hvgdecision_reference_query_patch-0.10.1-py3-none-any.whl --no-deps
```

检查：

```python
import hvgdecision as hd
import hvgdecision_reference_query_patch as hrq

print(hd.__version__)
print(hrq.__version__)
```

## Reference 参数

Reference 中的值对应 `batch_key`：

```python
reference_donors = [
    "patient_101",
    "patient_1015",
    "patient_1016",
    "patient_1039",
    "patient_107",
    "patient_1244",
]

result = hd.refine(
    adata,
    batch_key="replicate",
    label_key="cell_type",
    counts=raw_counts,
    reference=reference_donors,
    n_hvg=2000,
    return_details=True,
)
```

这里 `result.removed_genes` 是 **Reference-side risk list**。在 Reference→Query 设计中，它不等于最终 Query 删除列表。

## 推荐：完整 Reference→Query 工作流

```python
import hvgdecision_reference_query_patch as hrq

reference_donors = [
    "patient_101",
    "patient_1015",
    "patient_1016",
    "patient_1039",
    "patient_107",
    "patient_1244",
]

query_donors = [
    "patient_1256",
    "patient_1488",
]

result = hrq.refine_reference_query(
    adata,
    batch_key="replicate",
    label_key="cell_type",
    counts=raw_counts,
    reference=reference_donors,
    query=query_donors,
    n_hvg=2000,
    output_dir="./HVGDecision_Kang_reference6_to_query2",
)

print(result)
print("Reference risk:", result.reference_risk_genes)
print("Final removed from Query:", result.removed_genes)
print("Final n:", result.final_n_hvg)

adata_hvg = result.adata
```

逻辑为：

```text
Reference donors
  ↓
独立 Seurat-v3 HVG2000
每次从 0.3 开始，失败才 0.5 → 0.7 → 1.0
  ↓
原 HVGDecision v0.10.0 hd.refine()
  ↓
Reference risk list

Query donors
  ↓
第二次独立 Seurat-v3 HVG2000
重新从 0.3 开始
  ↓
Reference risk ∩ Query HVG2000
  ↓
最终 Query refined panel
```

数学上：

```text
H_final = H_query \ (R_reference ∩ H_query)
```

## Kang2018 已复现案例

你当前实验实际得到：

```text
Reference: 8176 cells / 6 donors
Reference span: 0.3 fail → 0.5 success
Reference HVG: 2000
Reference risk: 19

Query: 4139 cells / 2 donors
Query span: 0.3 success
Query HVG: 2000

Reference risk ∩ Query HVG:
CSF2
GJB2
AC147651.3

Final Query panel: 1997
```

这三个基因不是硬编码在补丁里；它们来自当前安装的 HVGDecision v0.10.0 风险结果与 Query-specific HVG2000 的交集。

## 输出审计

工作目录会保存：

- `REFERENCE_span_fallback_audit.csv`
- `REFERENCE_hvg_2000.csv`
- `REFERENCE_risk_genes.csv`
- `QUERY_span_fallback_audit.csv`
- `QUERY_hvg_2000.csv`
- `REFERENCE_to_QUERY_transfer_audit.csv`
- `QUERY_removed_by_REFERENCE_risk.csv`
- `QUERY_refined_hvg_panel.csv`
- `seuratv3_span_fallback_audit.csv`
- `experiment_summary.csv`

## 为什么这是补丁而不是替换版主 wheel

当前会话没有 HVGDecision v0.10.0 的完整源码/原始 wheel。为了不把 v0.10.0 的 `donor_aware` 核心误替换成旧版算法，本发布物有意把风险学习委托给你已经安装的 v0.10.0。

如果要发布正式的 `hvgdecision-0.10.1-py3-none-any.whl` 主 wheel，应在 v0.10.0 原源码上把同样的 span fallback 和 Reference→Query wrapper 合并进去后再构建。
