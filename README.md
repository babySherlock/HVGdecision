# HVGDecision 0.9.0

HVGDecision 是单细胞 RNA-seq integration 之前的 **integration-aware HVG refinement** 工具。

0.9.0 开始必须显式选择算法 `mode`：

```text
within_domain
    同一组织 / 相近实验体系，多 donor / batch
    → Reference-only Within-domain V1

cross_domain
    不同 dataset / technology 的 Reference → Query
    → Cross-domain Rule V3
```

HVGDecision 本身不是 Harmony、CCA、BBKNN、Scanorama 或 scVI 的替代品。它先决定用于这些下游方法的 feature panel。

---

## 1. 安装

```bash
pip install -e .
```

或：

```bash
pip install .
```

检查版本：

```python
import hvgdecision as hd
print(hd.__version__)
print(hd.VALID_MODES)
```

---

# 2. 必须先选 mode

## A. within_domain

适合：

- 同一数据来源或近似 protocol；
- 多 donor / batch；
- 目标是去除在 biology 条件下仍表现出 donor/batch 风险的 HVG。

最终 Reference-only raw risk：

\[
R_g^{within}=Z(L_g)+0.75Z(I_g)-Z(B_g)
\]

最终 harmful gate：

\[
q_g\le0.05,
\quad Z(L_g)\ge1,
\quad Z_R(g)\ge1,
\quad Z(B_g)\le0,
\quad F_g^{boot}\ge0.80,
\quad \neg P_g
\]

其中 permutation 在 biology strata 内打乱 donor label，使用 per-gene leakage p-value + BH-FDR；marker/显式 biology protection 可以阻止删除。

## B. cross_domain

适合：

- cross-dataset；
- cross-technology；
- Reference 和 Query 的表达分布明显不同。

Query 的真实 cell-type label 不参与 feature selection。Query expression 只用于无标签的 Reference→Query shift。

\[
S_g=0.70E_g^*+0.30D_g^*
\]

\[
B_g=0.60B_g^{celltype}+0.25M_g+0.15P_g
\]

\[
R_g^{cross}
=
\min(R_g,S_g)
\left[0.75+0.25(1-|R_g-S_g|)\right]
(1-0.75B_g)
\]

然后从 Query HVG panel 中删除 top-k unprotected genes。

软件暴露：

```text
cross_domain_delete_budget = 5 / 10 / 20 / other explicit integer
```

论文主分析目前可使用保守的 `k=5`；参数仍然显式保留，避免把删除深度藏在代码里。

---

# 3. Python API：within_domain

```python
import scanpy as sc
import hvgdecision as hd

adata = sc.read_h5ad("your_data.h5ad")

study = hd.setup_reference_query(
    adata,
    mode="within_domain",          # 必须显式选择
    batch_key="donor",
    label_key="cell_type",
    reference=["D1", "D2", "D3", "D4", "D5", "D6"],
    query=["D7", "D8"],
    counts_layer="counts",
    output_dir="HVGDecision_results/PBMC",
)

result = study.run(
    hvg_method="seurat_v3",
    return_details=True,
)

print(result.recommended_n_hvg)
print(result.harmful_genes)
print(result.final_n_hvg)
```

`within_domain` 的 `study.run()` 会先做 Reference-only minimum-sufficient HVG budget，再做 harmful-gene refinement。

也可继续直接调用：

```python
result = study.find_best_hvg(return_details=True)
```

但 `find_best_hvg()` 只允许 `within_domain`。

---

# 4. Python API：cross_domain

例如 Skin 10X → Smart-seq2：

```python
study = hd.setup_reference_query(
    adata,
    mode="cross_domain",           # 必须显式选择
    batch_key="donor",
    label_key="cell_type",
    split_key="method",
    reference=["10X"],
    query=["smartseq2"],
    counts_layer="counts",
    output_dir="HVGDecision_results/Skin_cross_technology",
)

result = study.run(
    n_hvg=2000,
    hvg_method="seurat_v3",
    cross_domain_delete_budget=5,
    return_details=True,
)

print(result.harmful_genes)
print(result.final_n_hvg)
```

Cross-domain 流程自动执行：

```text
Query HVG panel
      +
Reference independent HVG panel
      ↓
Reference Within-domain V1 raw risk percentile
      +
Reference→Query label-free expression shift
      −
Reference biology / marker protection
      ↓
Cross-domain Rule V3 ranking
      ↓
remove top-k from Query HVG panel
```

Query label 可以存在于 `adata.obs` 供最终外部 benchmark 使用，但 HVGDecision 不读取它做 Cross-domain risk。

---

# 5. 已有外部 HVG 列表

如果 HVG 是 Seurat/R 或其它 Python 工具先选好的：

```python
result = study.refine_hvg(
    ["GENE1", "GENE2", "GENE3", ...],
    method_name="Seurat_v3",
    initial_n_hvg=2000,
    return_details=True,
)
```

`within_domain` 会把它当作冻结的 Reference base panel；`cross_domain` 会把它当作冻结的 Query base panel。

---

# 6. 主要输出

Python API 的最终结果：

```python
adata_final = result.adata

final_hvgs = adata_final.var_names[
    adata_final.var["highly_variable"]
].tolist()

harmful = adata_final.var_names[
    adata_final.var["hvgdecision_harmful"]
].tolist()
```

对象不会裁成只剩 HVG；所有 count-source genes 仍保留。

核心 `var` 字段：

```text
highly_variable
hvgdecision_candidate
hvgdecision_risk_score
hvgdecision_risk_flagged
hvgdecision_marker_protected
hvgdecision_harmful
hvgdecision_removed
hvgdecision_final
hvgdecision_reason
```

核心文件：

```text
01_refine/base_hvg_ranking.csv
01_refine/gene_risk_evidence.csv
01_refine/selected_risk_genes.csv
01_refine/refinement_audit.json
```

Cross-domain 另外输出：

```text
01_refine/QUERY_hvg_ranking.csv
01_refine/REFERENCE_hvg_ranking.csv
01_refine/REFERENCE_RULE_V1_full_audit_for_crossdomain.csv
01_refine/CROSSDOMAIN_reference_query_shift_audit.csv
01_refine/CROSSDOMAIN_RULE_V3_ranking.csv
01_refine/CROSSDOMAIN_RULE_V3_budget_gene_membership.csv
CROSSDOMAIN_RULE_V3_selection_manifest.csv
CROSSDOMAIN_RULE_V3_removal_audit.csv
```

---

# 7. CLI：mode 也必须明确

Within-domain：

```bash
hvgdecision refine \
  --config configs/within_domain_example.yml \
  --mode within_domain
```

Cross-domain：

```bash
hvgdecision refine \
  --config configs/cross_domain_example.yml \
  --mode cross_domain \
  --delete-budget 5
```

也可以把 `mode:` 写在 YAML 里；如果 YAML 和 CLI 同时给出，以 CLI 为准。

---

# 8. 论文 benchmark 的三维综合分数（可选）

0.9.0 附带 `three_domain_scores()`，只用于论文评价，不参与 feature selection。

三个维度：

```text
Transfer
    Macro-F1
    Balanced accuracy
    Rare-cell Macro-F1

Batch correction
    Dataset mixing
    Donor mixing within cell type

Biological fidelity
    Leiden ARI
    Leiden NMI
    Cell-type silhouette
```

固定 `dataset × integration_method` 内先对不同 feature panel 做 min-max，再计算：

\[
Overall=0.30Transfer+0.40Batch+0.30Biological\ fidelity
\]

示例：

```python
scores = hd.three_domain_scores(benchmark_table)
```

这应称为 **scIB-inspired three-domain composite score**，不是原版 scIB score。

---

# 9. 关键原则

```text
Biological informativeness != integration suitability
```

Within-domain 与 Cross-domain 是两个不同的问题，因此 0.9.0 不再让软件隐式猜测算法，而要求用户明确选择 `mode`。
